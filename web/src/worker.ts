// Web Worker: loads vikngs-core.wasm, streams the uploaded files into MEMFS,
// runs runVikNGS, and posts a completion message back to the UI thread.

import type { RunRequest, SimRunRequest, UiToWorker, WorkerMessage } from "./types";

// The Emscripten-generated module is UMD/CJS and not module-import-friendly
// under Vite's ESM pipeline. Load it via importScripts in a classic worker,
// or dynamic import of the public URL. We expose it on globalThis for
// simplicity.
declare global {
    var createVikNGS: ((opts?: object) => Promise<VikNGSModule>) | undefined;
}

interface VikNGSModule {
    FS: {
        mkdir: (path: string) => void;
        writeFile: (path: string, data: Uint8Array) => void;
        analyzePath?: (path: string) => unknown;
    };
    print?: (t: string) => void;
    printErr?: (t: string) => void;
    runVikNGS: (req: unknown) => {
        rows: { size: () => number; get: (i: number) => {
            chrom: string; pos: number; ref: string; alt: string;
            pvalue: number; testDesc: string;
        }};
        evaluationTime: number;
        variantsParsed: number;
        errorMessage: string;
    };
    VectorSimGroup: new () => {
        push_back: (g: unknown) => void;
    };
    runSimulation: (req: unknown) => {
        rows: { size: () => number; get: (i: number) => {
            stepIdx: number; sampleSize: number; testIdx: number;
            statName: string; genotypeSource: string;
            variantIdx: number; pvalue: number;
        }};
        processingTime: number;
        evaluationTime: number;
        variantsParsed: number;
        errorMessage: string;
    };
    hello: (name: string) => string;
}

function post(msg: WorkerMessage) {
    (self as unknown as Worker).postMessage(msg);
}

async function loadModule(): Promise<VikNGSModule> {
    post({ kind: "log", level: "info", text: "Loading vikngs-core.wasm…" });
    // Load the UMD module via importScripts (classic-worker context) or
    // via a fetch+eval shim. Since Vite's worker format: "es" is ESM-only,
    // we fetch the file and evaluate via `new Function` to expose the
    // factory as a local variable.
    const cjsUrl = new URL("../public/vikngs-core.cjs", import.meta.url).href;
    const wasmUrl = new URL("../public/vikngs-core.wasm", import.meta.url).href;

    const cjsText = await (await fetch(cjsUrl)).text();
    const wasmBytes = new Uint8Array(await (await fetch(wasmUrl)).arrayBuffer());

    // The UMD wrapper assigns `module.exports = createVikNGS`. Simulate that.
    const moduleHarness: { exports: unknown } = { exports: {} };
    const fn = new Function("module", "exports", cjsText + "\nreturn module.exports;");
    const createVikNGS = fn(moduleHarness, moduleHarness.exports) as
        (opts?: object) => Promise<VikNGSModule>;

    // Forward the core's stdout/stderr (printInfo/printWarning/printError
    // in src/cmd/Log.cpp) to the UI log pane so users see live progress.
    return await createVikNGS({
        wasmBinary: wasmBytes,
        print:    (t: string) => post({ kind: "log", level: "info",  text: t }),
        printErr: (t: string) => {
            // Emscripten sends its own "warning: …" diagnostics through
            // printErr too — those aren't errors from our core.
            const level = t.startsWith("[ERROR]") ? "error" : "info";
            post({ kind: "log", level, text: t });
        },
    });
}

async function runSim(Module: VikNGSModule, req: SimRunRequest) {
    post({ kind: "log", level: "info", text: `Simulating (stat=${req.statistic}, nsnp=${req.nsnp}, steps=${req.steps})…` });

    const groupVec = new Module.VectorSimGroup();
    for (const g of req.groups) {
        groupVec.push_back({
            n: g.n,
            nIncrement: g.nIncrement,
            isCase: g.isCase,
            normalMean: g.normalMean,
            normalSd: g.normalSd,
            meanDepth: g.meanDepth,
            sdDepth: g.sdDepth,
            errorRate: g.errorRate,
            readDepth: g.readDepth,
        });
    }

    const cppReq = {
        nsnp: req.nsnp,
        effectSize: req.effectSize,
        mafMin: req.mafMin,
        mafMax: req.mafMax,
        steps: req.steps,
        family: req.family,
        statistic: req.statistic,
        collapse: req.collapse,
        nboot: req.nboot,
        stopEarly: req.stopEarly,
        covariate: req.covariate,
        corX: req.corX,
        groups: groupVec,
        seed: req.seed,
    };

    const result = Module.runSimulation(cppReq);
    if (result.errorMessage) {
        post({ kind: "error", message: result.errorMessage });
        return;
    }

    const rows = [];
    for (let i = 0; i < result.rows.size(); i++) {
        const r = result.rows.get(i);
        rows.push({
            stepIdx: r.stepIdx, sampleSize: r.sampleSize, testIdx: r.testIdx,
            statName: r.statName, genotypeSource: r.genotypeSource,
            variantIdx: r.variantIdx, pvalue: r.pvalue,
        });
    }

    post({
        kind: "sim-done",
        rows,
        steps: req.steps,
        processingTime: result.processingTime,
        evaluationTime: result.evaluationTime,
        variantsParsed: result.variantsParsed,
    });
}

self.onmessage = async (ev: MessageEvent<UiToWorker>) => {
    const req = ev.data;
    try {
        const Module = await loadModule();
        post({ kind: "log", level: "info", text: "Module ready. Staging inputs…" });

        if ((req as SimRunRequest).kind === "sim") {
            await runSim(Module, req as SimRunRequest);
            return;
        }

        const aReq = req as RunRequest;

        // Ensure /work exists in MEMFS.
        try { Module.FS.mkdir("/work"); } catch { /* already exists */ }

        Module.FS.writeFile("/work/input.vcf", new Uint8Array(await aReq.vcf.arrayBuffer()));
        Module.FS.writeFile("/work/sample.txt", new Uint8Array(await aReq.sample.arrayBuffer()));
        let bedPath = "";
        if (aReq.bed) {
            Module.FS.writeFile("/work/regions.bed", new Uint8Array(await aReq.bed.arrayBuffer()));
            bedPath = "/work/regions.bed";
        }

        post({ kind: "log", level: "info", text: `Running analysis (stat=${aReq.statistic}, gt=${aReq.genotype})…` });

        const cppReq = {
            vcfPath: "/work/input.vcf",
            samplePath: "/work/sample.txt",
            bedPath,
            outputDir: "/work",
            maf: aReq.maf,
            depth: aReq.depth,
            missing: aReq.missing,
            mustPass: aReq.mustPass,
            chrFilter: aReq.chrFilter,
            fromPos: aReq.fromPos,
            toPos: aReq.toPos,
            statistic: aReq.statistic,
            genotype: aReq.genotype,
            nboot: aReq.nboot,
            stopEarly: aReq.stopEarly,
            collapseMode: aReq.collapseMode,
            collapseK: aReq.collapseK,
            batchSize: aReq.batchSize,
            threads: 1,
        };

        const result = Module.runVikNGS(cppReq);
        if (result.errorMessage) {
            post({ kind: "error", message: result.errorMessage });
            return;
        }

        const rows = [];
        for (let i = 0; i < result.rows.size(); i++) {
            const r = result.rows.get(i);
            rows.push({
                chrom: r.chrom, pos: r.pos, ref: r.ref, alt: r.alt,
                pvalue: r.pvalue, testDesc: r.testDesc,
            });
        }

        post({
            kind: "done",
            rows,
            variantsParsed: result.variantsParsed,
            evaluationTime: result.evaluationTime,
        });
    } catch (e) {
        post({ kind: "error", message: e instanceof Error ? e.message : String(e) });
    }
};
