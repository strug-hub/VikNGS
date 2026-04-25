// Emscripten bindings for VikNGS.
//
// JS passes a Request value-struct and paths into MEMFS. This binding
// translates to the existing Request object, calls startVikNGS, and packages
// the resulting p-values into a JS-accessible array of plain-object rows.
//
// MEMFS model: JS is expected to FS.writeFile() the VCF / sample / BED
// contents to a path (conventionally /work/input.vcf etc.) before calling
// runVikNGS. The output directory is /work and is expected to exist.
//
// This is MVP-scope: one test at a time, no simulation, up to ~2GB input.
#include "../vikNGS.h"
#include "../Request.h"
#include "../Enum/Statistic.h"
#include "../Enum/GenotypeSource.h"
#include "../Enum/Variance.h"
#include "../Enum/CollapseType.h"
#include "../Enum/Family.h"
#include "../Enum/Depth.h"
#include "../Math/Math.h"
#include "../Simulation/Simulation.h"

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <emscripten/console.h>

#include <cmath>
#include <string>
#include <vector>

#define TRACE(msg) emscripten_console_log("[bindings] " msg)

using emscripten::class_;
using emscripten::function;
using emscripten::register_vector;
using emscripten::value_object;

// JS-facing request shape. Mirrors the flags that src/cmd/vikNGScmd.cpp sets.
struct JsRequest {
    std::string vcfPath;
    std::string samplePath;
    std::string bedPath;      // empty string = no BED
    std::string outputDir;    // required writable dir (use /work under MEMFS)

    // filters
    double maf = 0.05;
    int    depth = 30;
    double missing = 0.1;
    bool   mustPass = false;
    std::string chrFilter;    // empty = no chromosome filter
    int    fromPos = -1;
    int    toPos = -1;

    // test config (single test for M1 — matches CLI's single-test behavior)
    std::string statistic = "common";  // "common" | "cast" | "skat"
    std::string genotype  = "expected"; // "expected" | "call" | "vcf"

    // bootstrapping (common-variant tests default to 1 = asymptotic)
    int    nboot = 1;
    bool   stopEarly = false;

    // BED collapsing
    std::string collapseMode;  // "", "gene", "exon", "k"
    int    collapseK = -1;

    // performance
    int    batchSize = 1000;
    int    threads = 1;
};

// One row per variant per test (matches CLI pvalues_*.txt output).
struct JsResultRow {
    std::string chrom;
    int         pos;
    std::string ref;
    std::string alt;
    double      pvalue;
    std::string testDesc;  // short description: "No vRVS SKAT" etc.
};

struct JsResults {
    std::vector<JsResultRow> rows;
    double                    evaluationTime;
    int                       variantsParsed;
    std::string               errorMessage;  // empty on success
};

// One row per sample for the drill-down table behind a Manhattan point.
struct JsSampleGenotype {
    int    sampleIdx;
    int    group;
    double phenotype;
    // Genotype-source dosage vectors. Each entry is a 2-element [source-name,
    // dosage] pair packed into a JsSampleDosage so embind can serialise it.
    double trueDosage     = std::nan("");
    double expectedDosage = std::nan("");
    double callDosage     = std::nan("");
    double vcfDosage      = std::nan("");
};

struct JsAnalysisDetail {
    std::string chrom;
    int         pos;
    std::string ref;
    std::string alt;
    std::vector<JsSampleGenotype> samples;
    std::string errorMessage;
};

// ---------------------------------------------------------------------------
// Helpers that map JS strings to enums. Throws std::runtime_error on
// unknown values; the emscripten exception propagation turns that into a
// JS-side reject/throw.
// ---------------------------------------------------------------------------
static GenotypeSource parseGenotypeSource(const std::string& s) {
    if (s == "expected") return GenotypeSource::EXPECTED;
    if (s == "call")     return GenotypeSource::CALL;
    if (s == "vcf")      return GenotypeSource::VCF_CALL;
    throw std::runtime_error("unknown genotype source: " + s);
}

// Mirror the quirky statistic/genotype → (Statistic, Variance) mapping that
// src/cmd/vikNGScmd.cpp applies. "common" + expected → SKAT, "common" + vcf →
// SKAT+RVS, "common" + call → COMMON. Rare tests use their named statistic.
static TestSettings buildTestSettings(const std::string& statName, GenotypeSource gt) {
    if (statName == "cast") {
        if (gt == GenotypeSource::CALL)     return TestSettings(GenotypeSource::CALL,     Statistic::CAST, Variance::REGULAR);
        if (gt == GenotypeSource::VCF_CALL) return TestSettings(GenotypeSource::VCF_CALL, Statistic::CAST, Variance::RVS);
        return                                    TestSettings(GenotypeSource::EXPECTED, Statistic::CAST, Variance::REGULAR);
    }
    if (statName == "skat") {
        if (gt == GenotypeSource::CALL)     return TestSettings(GenotypeSource::CALL,     Statistic::SKAT, Variance::REGULAR);
        if (gt == GenotypeSource::VCF_CALL) return TestSettings(GenotypeSource::VCF_CALL, Statistic::SKAT, Variance::RVS);
        return                                    TestSettings(GenotypeSource::EXPECTED, Statistic::SKAT, Variance::REGULAR);
    }
    // "common" — CLI default path
    if (gt == GenotypeSource::CALL)     return TestSettings(GenotypeSource::CALL,     Statistic::COMMON, Variance::REGULAR);
    if (gt == GenotypeSource::VCF_CALL) return TestSettings(GenotypeSource::VCF_CALL, Statistic::SKAT,   Variance::RVS);
    return                                    TestSettings(GenotypeSource::EXPECTED, Statistic::SKAT,   Variance::REGULAR);
}

// ---------------------------------------------------------------------------
// Cached state from the most recent runVikNGS, used by getAnalysisDetail
// for the Manhattan drill-down. Keyed by the flat row index emitted into
// JsResults.rows.
// ---------------------------------------------------------------------------
static Data g_lastAnalysis;
static bool g_hasLastAnalysis = false;
static std::vector<std::pair<int,int>> g_detailIndex;  // (vsIdx, variantInVsIdx) per result row

// ---------------------------------------------------------------------------
// Main entry: build a Request, call startVikNGS, unpack results.
// ---------------------------------------------------------------------------
static JsResults runVikNGS(JsRequest js) {
    JsResults out;
    try {
        Request req = getDefaultRequest();
        req.setInputFiles(js.vcfPath, js.samplePath);
        req.setOutputDir(js.outputDir);

        req.setMafCutOff(js.maf);
        req.setHighLowCutOff(js.depth);
        req.setMissingThreshold(js.missing);
        req.setMustPASS(js.mustPass);
        if (!js.chrFilter.empty()) req.setChromosomeFilter(js.chrFilter);
        if (js.fromPos >= 0) req.setMinPos(js.fromPos);
        if (js.toPos   >= 0) req.setMaxPos(js.toPos);

        if (js.nboot > 1) {
            req.setBootstrap(js.nboot);
            req.setStopEarly(js.stopEarly);
        }

        if (!js.bedPath.empty()) {
            req.setCollapseFile(js.bedPath);
            if (js.collapseMode == "gene")      req.setCollapseGene();
            else if (js.collapseMode == "exon") req.setCollapseExon();
        }
        if (js.collapseMode == "k" && js.collapseK >= 2) req.setCollapse(js.collapseK);

        GenotypeSource gt = parseGenotypeSource(js.genotype);
        req.addTest(buildTestSettings(js.statistic, gt));

        req.setBatchSize(js.batchSize);
        req.setNumberThreads(js.threads);
        // Keep per-sample genotypes available for the drill-down table
        // exposed via getAnalysisDetail.
        req.setRetainGenotypes(true);

        TRACE("calling startVikNGS");
        Data data = startVikNGS(req);
        TRACE("startVikNGS returned");

        out.evaluationTime = data.evaluationTime;
        out.variantsParsed = static_cast<int>(data.variantsParsed);

        auto tests = req.getTests();
        std::string testDesc = tests.empty() ? "" : tests[0].toShortString();

        // Reset detail index alongside the row vector so they stay in sync.
        g_detailIndex.clear();

        // Emit one row per underlying Variant inside each VariantSet with
        // pvals — matches the CLI's output format (one row per variant even
        // when collapsed into multi-variant sets).
        for (size_t vsIdx = 0; vsIdx < data.variants.size(); vsIdx++) {
            auto& vs = data.variants[vsIdx];
            if (vs.nPvals() == 0) continue;
            double p = vs.getPval(0);
            auto* variants = vs.getVariants();
            for (size_t vIdx = 0; vIdx < variants->size(); vIdx++) {
                auto& v = (*variants)[vIdx];
                if (!v.isValid()) continue;
                JsResultRow row;
                row.chrom    = v.getChromosome();
                row.pos      = v.getPosition();
                row.ref      = v.getRef();
                row.alt      = v.getAlt();
                row.pvalue   = p;
                row.testDesc = testDesc;
                out.rows.push_back(std::move(row));
                g_detailIndex.emplace_back(static_cast<int>(vsIdx), static_cast<int>(vIdx));
            }
        }
        // Cache the full Data so getAnalysisDetail can return per-sample
        // genotypes without rerunning the analysis.
        g_lastAnalysis = std::move(data);
        g_hasLastAnalysis = true;
        TRACE("result packing done");
    } catch (const std::exception& e) {
        emscripten_console_error(e.what());
        out.errorMessage = e.what();
    } catch (...) {
        emscripten_console_error("[bindings] unknown (non-std) exception");
        out.errorMessage = "unknown error";
    }
    return out;
}

// ---------------------------------------------------------------------------
// Per-row drill-down. Pulls per-sample dosages from the cached `Data`
// for the variant referenced by `rowIdx` (the same flat index used by
// JsResults.rows). Cheap to call repeatedly — no recomputation.
// ---------------------------------------------------------------------------
static JsAnalysisDetail getAnalysisDetail(int rowIdx) {
    JsAnalysisDetail out;
    if (!g_hasLastAnalysis) {
        out.errorMessage = "no analysis result cached; run an analysis first";
        return out;
    }
    if (rowIdx < 0 || static_cast<size_t>(rowIdx) >= g_detailIndex.size()) {
        out.errorMessage = "rowIdx out of range";
        return out;
    }
    auto [vsIdx, vIdx] = g_detailIndex[rowIdx];
    auto& vs = g_lastAnalysis.variants[vsIdx];
    auto* variants = vs.getVariants();
    auto& v = (*variants)[vIdx];

    out.chrom = v.getChromosome();
    out.pos   = v.getPosition();
    out.ref   = v.getRef();
    out.alt   = v.getAlt();

    // Pull per-sample vectors for whichever sources are populated.
    auto sources = v.getAllGenotypes();
    VectorXd* trueGT = nullptr;
    VectorXd* expGT  = nullptr;
    VectorXd* callGT = nullptr;
    VectorXd* vcfGT  = nullptr;
    for (auto src : sources) {
        VectorXd* g = v.getGenotype(src);
        if      (src == GenotypeSource::TRUEGT)   trueGT = g;
        else if (src == GenotypeSource::EXPECTED) expGT  = g;
        else if (src == GenotypeSource::CALL)     callGT = g;
        else if (src == GenotypeSource::VCF_CALL) vcfGT  = g;
    }

    VectorXd y = g_lastAnalysis.sampleInfo.getY();
    VectorXi g = g_lastAnalysis.sampleInfo.getG();

    // Pick whichever genotype vector exists to size the table.
    int n = 0;
    for (VectorXd* p : {trueGT, expGT, callGT, vcfGT}) if (p) { n = p->size(); break; }

    out.samples.reserve(n);
    for (int i = 0; i < n; i++) {
        JsSampleGenotype row;
        row.sampleIdx = i;
        row.group     = (i < g.size()) ? g[i] : -1;
        row.phenotype = (i < y.size()) ? y[i] : std::nan("");
        if (trueGT) row.trueDosage     = (*trueGT)[i];
        if (expGT)  row.expectedDosage = (*expGT)[i];
        if (callGT) row.callDosage     = (*callGT)[i];
        if (vcfGT)  row.vcfDosage      = (*vcfGT)[i];
        out.samples.push_back(std::move(row));
    }
    return out;
}

// ---------------------------------------------------------------------------
// Simulation bindings (MVP: binomial case-control, multi-step supported).
// ---------------------------------------------------------------------------
struct JsSimGroup {
    int    n = 500;
    int    nIncrement = 0;           // per-step sample-size increment
    bool   isCase = true;            // ignored when family=normal
    double normalMean = 0.0;         // NORMAL family per-group phenotype mean
    double normalSd = 1.0;           // NORMAL family per-group phenotype SD
    double meanDepth = 20.0;
    double sdDepth = 2.0;
    double errorRate = 0.01;
    std::string readDepth = "high";  // "high" | "low"
};

struct JsSimRequest {
    int    nsnp = 100;
    double effectSize = 1.0;         // OR (binomial) or R^2 (normal)
    double mafMin = 0.05;
    double mafMax = 0.5;
    int    steps = 1;

    std::string family = "binomial"; // "binomial" | "normal"
    std::string statistic = "common";
    int    collapse = 1;
    int    nboot = 1;
    bool   stopEarly = false;

    // Covariate simulation. covariate < 0 disables; corX=true correlates
    // covariate with genotype X, false with phenotype Y.
    double covariate = -1.0;
    bool   corX = true;

    std::vector<JsSimGroup> groups;

    int    seed = 0;                 // 0 = nondeterministic
};

struct JsSimPvalRow {
    int         stepIdx;
    int         sampleSize;
    int         testIdx;
    std::string statName;
    std::string genotypeSource;      // "true" | "expected" | "call"
    int         variantIdx;
    double      pvalue;
};

struct JsSimResult {
    std::vector<JsSimPvalRow> rows;
    double                    processingTime = 0.0;
    double                    evaluationTime = 0.0;
    int                       variantsParsed = 0;
    std::string               errorMessage;
};

static const char* gtSourceName(GenotypeSource g) {
    switch (g) {
        case GenotypeSource::EXPECTED: return "expected";
        case GenotypeSource::TRUEGT:   return "true";
        case GenotypeSource::CALL:     return "call";
        case GenotypeSource::VCF_CALL: return "vcf";
        default:                       return "none";
    }
}

static const char* statShortName(Statistic s) {
    switch (s) {
        case Statistic::COMMON: return "common";
        case Statistic::CAST:   return "cast";
        case Statistic::SKAT:   return "skat";
        case Statistic::CALPHA: return "calpha";
        default:                return "none";
    }
}

static JsSimResult runSimulation(JsSimRequest js) {
    JsSimResult out;
    try {
        if (js.seed != 0) setRandomSeed(static_cast<uint64_t>(js.seed));

        SimulationRequest req;
        req.nsnp         = js.nsnp;
        req.effectSize   = js.effectSize;
        req.mafMin       = js.mafMin;
        req.mafMax       = js.mafMax;
        req.steps        = js.steps < 1 ? 1 : js.steps;
        req.collapse     = js.collapse;
        req.nboot        = js.nboot;
        req.useBootstrap = js.nboot > 1;
        req.stopEarly    = js.stopEarly;
        req.nthreads     = 1;          // WASM MVP: single-threaded
        req.covariate    = js.covariate;
        req.corX         = js.corX;

        if (js.family == "binomial")    req.family = Family::BINOMIAL;
        else if (js.family == "normal") req.family = Family::NORMAL;
        else throw std::runtime_error("unknown family: " + js.family);

        if (js.statistic == "common")      req.testStatistic = Statistic::COMMON;
        else if (js.statistic == "cast")   req.testStatistic = Statistic::CAST;
        else if (js.statistic == "skat")   req.testStatistic = Statistic::SKAT;
        else if (js.statistic == "calpha") req.testStatistic = Statistic::CALPHA;
        else throw std::runtime_error("unknown statistic: " + js.statistic);

        for (size_t i = 0; i < js.groups.size(); i++) {
            const JsSimGroup& jg = js.groups[i];
            SimulationRequestGroup g;
            g.index       = static_cast<int>(i);
            g.n           = jg.n;
            g.n_increment = jg.nIncrement;
            g.family      = req.family;
            g.isCase      = jg.isCase;
            g.normalMean  = jg.normalMean;
            g.normalSd    = jg.normalSd;
            g.meanDepth   = jg.meanDepth;
            g.sdDepth     = jg.sdDepth;
            g.errorRate   = jg.errorRate;
            g.readDepth   = (jg.readDepth == "low") ? Depth::LOW : Depth::HIGH;
            req.groups.push_back(g);
        }

        req.validate();

        TRACE("calling startSimulation");
        Data data = startSimulation(req);
        TRACE("startSimulation returned");

        out.processingTime = data.processingTime;
        out.evaluationTime = data.evaluationTime;
        out.variantsParsed = static_cast<int>(data.variants.size());

        int totalTests = static_cast<int>(data.tests.size());
        int testsPerStep = req.steps > 0 ? (totalTests / req.steps) : totalTests;
        if (testsPerStep <= 0) testsPerStep = 1;

        for (size_t i = 0; i < data.variants.size(); i++) {
            int np = data.variants[i].nPvals();
            for (int k = 0; k < np; k++) {
                TestSettings& ts = data.tests[k];
                int stepIdx = k / testsPerStep;
                JsSimPvalRow row;
                row.stepIdx        = stepIdx;
                row.sampleSize     = req.nsamp(stepIdx);
                row.testIdx        = k % testsPerStep;
                row.statName       = statShortName(ts.getStatistic());
                row.genotypeSource = gtSourceName(ts.getGenotype());
                row.variantIdx     = static_cast<int>(i);
                row.pvalue         = data.variants[i].getPval(k);
                out.rows.push_back(std::move(row));
            }
        }
        TRACE("sim result packing done");
    } catch (const std::exception& e) {
        emscripten_console_error(e.what());
        out.errorMessage = e.what();
    } catch (...) {
        emscripten_console_error("[bindings] unknown (non-std) exception in sim");
        out.errorMessage = "unknown error";
    }
    return out;
}

// ---------------------------------------------------------------------------
// Streaming variant of runVikNGS. Same JsRequest fields, plus a JS reader
// object whose `readNextChunk()` method returns
//   { done: bool, bytes: Uint8Array }
// synchronously. In a Web Worker, JS can build that on top of FileReaderSync
// to pump bytes from a Blob without staging the whole VCF into MEMFS — so
// multi-GB inputs work without hitting the 2GB MEMFS ceiling.
//
// `vcfPath` should still point to a small header-only file in MEMFS so
// SampleParser can map sample IDs.
// ---------------------------------------------------------------------------
static JsResults runVikNGSStreaming(JsRequest js, emscripten::val jsReader) {
    JsResults out;
    try {
        Request req = getDefaultRequest();
        req.setInputFiles(js.vcfPath, js.samplePath);
        req.setOutputDir(js.outputDir);

        req.setMafCutOff(js.maf);
        req.setHighLowCutOff(js.depth);
        req.setMissingThreshold(js.missing);
        req.setMustPASS(js.mustPass);
        if (!js.chrFilter.empty()) req.setChromosomeFilter(js.chrFilter);
        if (js.fromPos >= 0) req.setMinPos(js.fromPos);
        if (js.toPos   >= 0) req.setMaxPos(js.toPos);

        if (js.nboot > 1) {
            req.setBootstrap(js.nboot);
            req.setStopEarly(js.stopEarly);
        }

        if (!js.bedPath.empty()) {
            req.setCollapseFile(js.bedPath);
            if (js.collapseMode == "gene")      req.setCollapseGene();
            else if (js.collapseMode == "exon") req.setCollapseExon();
        }
        if (js.collapseMode == "k" && js.collapseK >= 2) req.setCollapse(js.collapseK);

        GenotypeSource gt = parseGenotypeSource(js.genotype);
        req.addTest(buildTestSettings(js.statistic, gt));

        req.setBatchSize(js.batchSize);
        req.setNumberThreads(js.threads);
        req.setRetainGenotypes(true);

        // Wrap the JS reader as a synchronous ChunkCallback.
        // The Reader returns { done: bool, bytes: Uint8Array }.
        emscripten::val reader = jsReader;
        req.setVcfStreamSource([reader](std::string& out) -> bool {
            emscripten::val chunk = reader.call<emscripten::val>("readNextChunk");
            if (chunk["done"].as<bool>()) return false;
            emscripten::val bytes = chunk["bytes"];
            int len = bytes["length"].as<int>();
            if (len <= 0) return false;
            // Resize once and copy via the typed-memory helper.
            out.resize(static_cast<size_t>(len));
            emscripten::val view{ emscripten::typed_memory_view(
                static_cast<size_t>(len),
                reinterpret_cast<uint8_t*>(out.data())) };
            view.call<void>("set", bytes);
            return true;
        });

        TRACE("calling startVikNGS (streaming)");
        Data data = startVikNGS(req);
        TRACE("startVikNGS returned");

        out.evaluationTime = data.evaluationTime;
        out.variantsParsed = static_cast<int>(data.variantsParsed);

        auto tests = req.getTests();
        std::string testDesc = tests.empty() ? "" : tests[0].toShortString();

        g_detailIndex.clear();
        for (size_t vsIdx = 0; vsIdx < data.variants.size(); vsIdx++) {
            auto& vs = data.variants[vsIdx];
            if (vs.nPvals() == 0) continue;
            double p = vs.getPval(0);
            auto* variants = vs.getVariants();
            for (size_t vIdx = 0; vIdx < variants->size(); vIdx++) {
                auto& v = (*variants)[vIdx];
                if (!v.isValid()) continue;
                JsResultRow row;
                row.chrom    = v.getChromosome();
                row.pos      = v.getPosition();
                row.ref      = v.getRef();
                row.alt      = v.getAlt();
                row.pvalue   = p;
                row.testDesc = testDesc;
                out.rows.push_back(std::move(row));
                g_detailIndex.emplace_back(static_cast<int>(vsIdx), static_cast<int>(vIdx));
            }
        }
        g_lastAnalysis = std::move(data);
        g_hasLastAnalysis = true;
    } catch (const std::exception& e) {
        emscripten_console_error(e.what());
        out.errorMessage = e.what();
    } catch (...) {
        emscripten_console_error("[bindings] unknown (non-std) exception in streaming sim");
        out.errorMessage = "unknown error";
    }
    return out;
}

// Backwards-compat smoke ping used by the original Phase-C hello test.
static std::string hello(const std::string& name) {
    return "hello, " + name + ", from vikngs-core wasm";
}

EMSCRIPTEN_BINDINGS(vikngs_core) {
    value_object<JsRequest>("Request")
        .field("vcfPath",      &JsRequest::vcfPath)
        .field("samplePath",   &JsRequest::samplePath)
        .field("bedPath",      &JsRequest::bedPath)
        .field("outputDir",    &JsRequest::outputDir)
        .field("maf",          &JsRequest::maf)
        .field("depth",        &JsRequest::depth)
        .field("missing",      &JsRequest::missing)
        .field("mustPass",     &JsRequest::mustPass)
        .field("chrFilter",    &JsRequest::chrFilter)
        .field("fromPos",      &JsRequest::fromPos)
        .field("toPos",        &JsRequest::toPos)
        .field("statistic",    &JsRequest::statistic)
        .field("genotype",     &JsRequest::genotype)
        .field("nboot",        &JsRequest::nboot)
        .field("stopEarly",    &JsRequest::stopEarly)
        .field("collapseMode", &JsRequest::collapseMode)
        .field("collapseK",    &JsRequest::collapseK)
        .field("batchSize",    &JsRequest::batchSize)
        .field("threads",      &JsRequest::threads);

    value_object<JsResultRow>("ResultRow")
        .field("chrom",    &JsResultRow::chrom)
        .field("pos",      &JsResultRow::pos)
        .field("ref",      &JsResultRow::ref)
        .field("alt",      &JsResultRow::alt)
        .field("pvalue",   &JsResultRow::pvalue)
        .field("testDesc", &JsResultRow::testDesc);

    register_vector<JsResultRow>("VectorResultRow");

    value_object<JsResults>("Results")
        .field("rows",           &JsResults::rows)
        .field("evaluationTime", &JsResults::evaluationTime)
        .field("variantsParsed", &JsResults::variantsParsed)
        .field("errorMessage",   &JsResults::errorMessage);

    // -----------------------------------------------------------------------
    // Simulation surface
    // -----------------------------------------------------------------------
    value_object<JsSimGroup>("SimGroup")
        .field("n",          &JsSimGroup::n)
        .field("nIncrement", &JsSimGroup::nIncrement)
        .field("isCase",     &JsSimGroup::isCase)
        .field("normalMean", &JsSimGroup::normalMean)
        .field("normalSd",   &JsSimGroup::normalSd)
        .field("meanDepth",  &JsSimGroup::meanDepth)
        .field("sdDepth",    &JsSimGroup::sdDepth)
        .field("errorRate",  &JsSimGroup::errorRate)
        .field("readDepth",  &JsSimGroup::readDepth);

    register_vector<JsSimGroup>("VectorSimGroup");

    value_object<JsSimRequest>("SimRequest")
        .field("nsnp",       &JsSimRequest::nsnp)
        .field("effectSize", &JsSimRequest::effectSize)
        .field("mafMin",     &JsSimRequest::mafMin)
        .field("mafMax",     &JsSimRequest::mafMax)
        .field("steps",      &JsSimRequest::steps)
        .field("family",     &JsSimRequest::family)
        .field("statistic",  &JsSimRequest::statistic)
        .field("collapse",   &JsSimRequest::collapse)
        .field("nboot",      &JsSimRequest::nboot)
        .field("stopEarly",  &JsSimRequest::stopEarly)
        .field("covariate",  &JsSimRequest::covariate)
        .field("corX",       &JsSimRequest::corX)
        .field("groups",     &JsSimRequest::groups)
        .field("seed",       &JsSimRequest::seed);

    value_object<JsSimPvalRow>("SimPvalRow")
        .field("stepIdx",        &JsSimPvalRow::stepIdx)
        .field("sampleSize",     &JsSimPvalRow::sampleSize)
        .field("testIdx",        &JsSimPvalRow::testIdx)
        .field("statName",       &JsSimPvalRow::statName)
        .field("genotypeSource", &JsSimPvalRow::genotypeSource)
        .field("variantIdx",     &JsSimPvalRow::variantIdx)
        .field("pvalue",         &JsSimPvalRow::pvalue);

    register_vector<JsSimPvalRow>("VectorSimPvalRow");

    value_object<JsSimResult>("SimResult")
        .field("rows",           &JsSimResult::rows)
        .field("processingTime", &JsSimResult::processingTime)
        .field("evaluationTime", &JsSimResult::evaluationTime)
        .field("variantsParsed", &JsSimResult::variantsParsed)
        .field("errorMessage",   &JsSimResult::errorMessage);

    function("runSimulation", &runSimulation);

    function("runVikNGS", &runVikNGS);

    // -----------------------------------------------------------------------
    // Drill-down: per-sample table for one variant.
    // -----------------------------------------------------------------------
    value_object<JsSampleGenotype>("SampleGenotype")
        .field("sampleIdx",      &JsSampleGenotype::sampleIdx)
        .field("group",          &JsSampleGenotype::group)
        .field("phenotype",      &JsSampleGenotype::phenotype)
        .field("trueDosage",     &JsSampleGenotype::trueDosage)
        .field("expectedDosage", &JsSampleGenotype::expectedDosage)
        .field("callDosage",     &JsSampleGenotype::callDosage)
        .field("vcfDosage",      &JsSampleGenotype::vcfDosage);

    register_vector<JsSampleGenotype>("VectorSampleGenotype");

    value_object<JsAnalysisDetail>("AnalysisDetail")
        .field("chrom",        &JsAnalysisDetail::chrom)
        .field("pos",          &JsAnalysisDetail::pos)
        .field("ref",          &JsAnalysisDetail::ref)
        .field("alt",          &JsAnalysisDetail::alt)
        .field("samples",      &JsAnalysisDetail::samples)
        .field("errorMessage", &JsAnalysisDetail::errorMessage);

    function("getAnalysisDetail", &getAnalysisDetail);
    function("runVikNGSStreaming", &runVikNGSStreaming);

    function("hello",     &hello);
}
