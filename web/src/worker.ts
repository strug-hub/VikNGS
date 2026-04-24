// Web Worker: loads vikngs-core.wasm, streams the uploaded files into MEMFS,
// runs runVikNGS, and posts a completion message back to the UI thread.

import type { RunRequest, WorkerMessage } from "./types";

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
    runVikNGS: (req: unknown) => {
        rows: { size: () => number; get: (i: number) => {
            chrom: string; pos: number; ref: string; alt: string;
            pvalue: number; testDesc: string;
        }};
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

    return await createVikNGS({ wasmBinary: wasmBytes });
}

self.onmessage = async (ev: MessageEvent<RunRequest>) => {
    const req = ev.data;
    try {
        const Module = await loadModule();
        post({ kind: "log", level: "info", text: "Module ready. Staging inputs…" });

        // Ensure /work exists in MEMFS.
        try { Module.FS.mkdir("/work"); } catch { /* already exists */ }

        Module.FS.writeFile("/work/input.vcf", new Uint8Array(await req.vcf.arrayBuffer()));
        Module.FS.writeFile("/work/sample.txt", new Uint8Array(await req.sample.arrayBuffer()));
        let bedPath = "";
        if (req.bed) {
            Module.FS.writeFile("/work/regions.bed", new Uint8Array(await req.bed.arrayBuffer()));
            bedPath = "/work/regions.bed";
        }

        post({ kind: "log", level: "info", text: `Running analysis (stat=${req.statistic}, gt=${req.genotype})…` });

        const cppReq = {
            vcfPath: "/work/input.vcf",
            samplePath: "/work/sample.txt",
            bedPath,
            outputDir: "/work",
            maf: req.maf,
            depth: req.depth,
            missing: req.missing,
            mustPass: req.mustPass,
            chrFilter: req.chrFilter,
            fromPos: req.fromPos,
            toPos: req.toPos,
            statistic: req.statistic,
            genotype: req.genotype,
            nboot: req.nboot,
            stopEarly: req.stopEarly,
            collapseMode: req.collapseMode,
            collapseK: req.collapseK,
            batchSize: req.batchSize,
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
