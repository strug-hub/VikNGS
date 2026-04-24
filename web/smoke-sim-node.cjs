// Sim smoke test. Cross-toolchain RNG divergence (libc++ under Emscripten vs
// libstdc++ under native gcc) means we can't compare p-values bit-for-bit
// against sim-cli — std::normal_distribution etc. differ by stdlib. Instead:
//   1. WASM runSimulation with a fixed seed produces identical output on
//      repeat runs (determinism gate).
//   2. Under the null (effect=1.0) p-values live in [0,1] and pass a loose
//      uniformity sanity check (KS on empirical CDF).
// Any numerical regression vs the native sim is already covered by the
// existing stats-test suite driven through sim-cli.
const createVikNGS = require("./public/vikngs-core.cjs");
const { readFileSync } = require("node:fs");
const { join, dirname } = require("node:path");

const SEED = 42;
const NSNP = 200;

function newGroups(Module) {
    const g = new Module.VectorSimGroup();
    g.push_back({ n: 200, nIncrement: 0, isCase: true,  meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" });
    g.push_back({ n: 200, nIncrement: 0, isCase: false, meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" });
    return g;
}

function buildReq(groups) {
    return {
        nsnp: NSNP,
        effectSize: 1.0,
        mafMin: 0.05,
        mafMax: 0.5,
        steps: 1,
        family: "binomial",
        statistic: "common",
        collapse: 1,
        nboot: 1,
        stopEarly: false,
        groups,
        seed: SEED,
    };
}

function unpack(result) {
    const rows = [];
    for (let i = 0; i < result.rows.size(); i++) rows.push(result.rows.get(i));
    return rows;
}

(async () => {
    const wasmBinary = readFileSync(join(__dirname, "public/vikngs-core.wasm"));
    const Module = await createVikNGS({ wasmBinary });

    // --- Gate 1: determinism ---
    const a = unpack(Module.runSimulation(buildReq(newGroups(Module))));
    const b = unpack(Module.runSimulation(buildReq(newGroups(Module))));
    if (a.length !== b.length) { console.error(`length ${a.length} vs ${b.length}`); process.exit(1); }
    let det = 0;
    for (let i = 0; i < a.length; i++) {
        if (a[i].pvalue !== b[i].pvalue) det++;
    }
    if (det) { console.error(`FAIL: ${det}/${a.length} differ on reseeded rerun`); process.exit(1); }
    console.log(`OK: ${a.length} rows reproduce under seed ${SEED}`);

    // --- Gate 2: under-null sanity ---
    for (const src of ["true", "expected", "call"]) {
        const ps = a.filter(r => r.genotypeSource === src).map(r => r.pvalue).sort((x, y) => x - y);
        if (ps.length === 0) continue;
        if (ps.some(p => p < 0 || p > 1 || Number.isNaN(p))) {
            console.error(`FAIL: ${src} p-values outside [0,1]`);
            process.exit(1);
        }
        // KS statistic against Uniform(0,1).
        let d = 0;
        for (let i = 0; i < ps.length; i++) {
            d = Math.max(d, Math.abs((i + 1) / ps.length - ps[i]), Math.abs(i / ps.length - ps[i]));
        }
        // Loose: 1.63/sqrt(n) = 99% critical under H0. We want the test to
        // flag if something's *very* wrong (e.g., all p-values near 0.5).
        const crit = 1.63 / Math.sqrt(ps.length);
        const pass = d < crit;
        console.log(`  ${src}: n=${ps.length} KS=${d.toFixed(4)} crit(99%)=${crit.toFixed(4)} ${pass ? "OK" : "WARN"}`);
        if (!pass) { console.error(`FAIL: ${src} p-values reject uniformity`); process.exit(1); }
    }
    console.log("sim smoke OK");
})();
