// Node CJS smoke test. The Emscripten module is UMD; we use require().
const createVikNGS = require("./public/vikngs-core.cjs");
const { readFileSync } = require("node:fs");
const { join, dirname } = require("node:path");

const repoRoot = dirname(__dirname);
const vcfBytes   = readFileSync(join(repoRoot, "example/example.vcf"));
const sampleText = readFileSync(join(repoRoot, "example/example_info.txt"));

(async () => {
    // Node 22's global fetch confuses Emscripten 3.1.6's Node detection; it
    // tries to fetch() the .wasm URL. Hand it the bytes directly.
    const wasmBinary = readFileSync(join(__dirname, "public/vikngs-core.wasm"));
    const Module = await createVikNGS({ wasmBinary });
    console.log("module loaded, hello() says:", Module.hello("test"));

    Module.FS.mkdir("/work");
    Module.FS.writeFile("/work/input.vcf",  vcfBytes);
    Module.FS.writeFile("/work/sample.txt", sampleText);

    const req = {
        vcfPath: "/work/input.vcf",
        samplePath: "/work/sample.txt",
        bedPath: "",
        outputDir: "/work",
        maf: 0.05, depth: 30, missing: 0.1, mustPass: false,
        chrFilter: "", fromPos: -1, toPos: -1,
        statistic: "common", genotype: "expected",
        nboot: 1, stopEarly: false,
        collapseMode: "", collapseK: -1,
        batchSize: 1000, threads: 1,
    };

    const result = Module.runVikNGS(req);
    if (result.errorMessage) {
        console.error("runVikNGS failed:", result.errorMessage);
        process.exit(1);
    }

    const rows = [];
    for (let i = 0; i < result.rows.size(); i++) rows.push(result.rows.get(i));
    console.log(`rows=${rows.length} variantsParsed=${result.variantsParsed} eval=${result.evaluationTime.toFixed(3)}s`);
    console.log("first 3 rows:");
    for (const r of rows.slice(0, 3)) console.log(" ", r);

    const golden = readFileSync(
        join(repoRoot, "tests/e2e/golden/common_default.pvalues.txt"), "utf8"
    ).trim().split("\n");
    if (rows.length !== golden.length) {
        console.error(`row-count mismatch: wasm=${rows.length} golden=${golden.length}`);
        process.exit(1);
    }
    // Golden p-values are formatted to 6 decimals by Variant::toString (via
    // std::to_string), so rtol=5e-6 accepts formatted-golden rounding drift
    // while still catching any real numerical divergence.
    const rtol = 5e-6;
    let drift = 0;
    for (let i = 0; i < golden.length; i++) {
        const [gchr, gpos, , , gp] = golden[i].split("\t");
        const r = rows[i];
        const dp = Math.abs(parseFloat(gp) - r.pvalue) / Math.max(Math.abs(parseFloat(gp)), 1e-30);
        if (r.chrom !== gchr)             { console.error(`row ${i}: chr ${r.chrom} vs ${gchr}`); drift++; }
        if (r.pos !== parseInt(gpos, 10)) { console.error(`row ${i}: pos ${r.pos} vs ${gpos}`); drift++; }
        if (dp > rtol)                    { console.error(`row ${i}: p ${r.pvalue} vs ${gp} (drift ${dp})`); drift++; }
    }
    if (drift === 0) {
        console.log(`OK: all ${rows.length} rows match golden within rtol=${rtol}`);
    } else {
        console.error(`FAIL: ${drift} discrepancies`);
        process.exit(1);
    }
})();
