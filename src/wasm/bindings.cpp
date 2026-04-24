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

        TRACE("calling startVikNGS");
        Data data = startVikNGS(req);
        TRACE("startVikNGS returned");

        out.evaluationTime = data.evaluationTime;
        out.variantsParsed = static_cast<int>(data.variantsParsed);

        auto tests = req.getTests();
        std::string testDesc = tests.empty() ? "" : tests[0].toShortString();

        // Emit one row per underlying Variant inside each VariantSet with
        // pvals — matches the CLI's output format (one row per variant even
        // when collapsed into multi-variant sets).
        for (auto& vs : data.variants) {
            if (vs.nPvals() == 0) continue;
            double p = vs.getPval(0);
            auto* variants = vs.getVariants();
            for (auto& v : *variants) {
                if (!v.isValid()) continue;
                JsResultRow row;
                row.chrom    = v.getChromosome();
                row.pos      = v.getPosition();
                row.ref      = v.getRef();
                row.alt      = v.getAlt();
                row.pvalue   = p;
                row.testDesc = testDesc;
                out.rows.push_back(std::move(row));
            }
        }
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

    function("runVikNGS", &runVikNGS);
    function("hello",     &hello);
}
