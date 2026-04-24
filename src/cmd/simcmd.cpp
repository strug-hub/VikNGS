// sim-cli: headless driver for the simulation engine.
// Builds a SimulationRequest from flags, runs startSimulation, emits
// one TSV row per (variant, test) with pvalue.
#include "../vikNGS.h"
#include "../Log.h"
#include "../Enum/Family.h"
#include "../Enum/Statistic.h"
#include "../Enum/GenotypeSource.h"
#include "../Enum/Variance.h"
#include "../Enum/Depth.h"
#include "../Math/Math.h"
#include "../gui/src/simulation/Simulation.h"
#include "CLI11.h"

#include <iostream>
#include <iomanip>
#include <string>

static const char* gtSourceName(GenotypeSource g) {
    switch (g) {
        case GenotypeSource::EXPECTED: return "expected";
        case GenotypeSource::TRUEGT:   return "true";
        case GenotypeSource::CALL:     return "call";
        case GenotypeSource::VCF_CALL: return "vcf";
        default:                       return "none";
    }
}

static const char* statName(Statistic s) {
    switch (s) {
        case Statistic::COMMON: return "common";
        case Statistic::CAST:   return "cast";
        case Statistic::SKAT:   return "skat";
        case Statistic::CALPHA: return "calpha";
        default:                return "none";
    }
}

int main(int argc, char* argv[]) {
    CLI::App app{ "sim-cli: headless VikNGS simulation driver" };

    std::string family = "binomial";      // binomial | normal
    std::string stat   = "common";        // common | cast | skat
    int    nsnp        = 100;
    double effect      = 1.0;             // OR for binomial, R^2 for normal
    double mafMin      = 0.05;
    double mafMax      = 0.5;
    int    collapse    = 1;
    int    nboot       = 1;
    int    ncase       = 500;
    int    ncontrol    = 500;
    double caseDepth   = 20.0;
    double caseDepthSd = 2.0;
    double caseErr     = 0.01;
    double ctrlDepth   = 20.0;
    double ctrlDepthSd = 2.0;
    double ctrlErr     = 0.01;
    int    caseHigh    = 1;               // 1 = Depth::HIGH, 0 = Depth::LOW
    int    ctrlHigh    = 1;
    bool   header      = false;
    uint64_t seed      = 0;               // 0 = nondeterministic (use random_device init)

    app.add_option("--family", family,       "binomial (default) or normal");
    app.add_option("--stat",   stat,         "common | cast | skat (default common)");
    app.add_option("--nsnp",   nsnp,         "Number of variants (independent replicates under null)");
    app.add_option("--effect", effect,       "OR (binomial) or R^2 (normal). 1.0/0.0 = null.");
    app.add_option("--maf-min", mafMin,      "Minimum MAF");
    app.add_option("--maf-max", mafMax,      "Maximum MAF");
    app.add_option("--collapse", collapse,   "Collapse size for rare tests");
    app.add_option("--nboot",  nboot,        "Bootstrap iterations (1 = asymptotic)");
    app.add_option("--ncase",  ncase,        "Case sample size");
    app.add_option("--ncontrol", ncontrol,   "Control sample size");
    app.add_option("--case-depth", caseDepth,"Case mean read depth");
    app.add_option("--case-depth-sd", caseDepthSd, "Case read depth SD");
    app.add_option("--case-error", caseErr,  "Case error rate");
    app.add_option("--control-depth", ctrlDepth, "Control mean read depth");
    app.add_option("--control-depth-sd", ctrlDepthSd, "Control read depth SD");
    app.add_option("--control-error", ctrlErr, "Control error rate");
    app.add_option("--case-high", caseHigh,  "1 = case is high-depth cohort (default), 0 = low");
    app.add_option("--control-high", ctrlHigh, "1 = control is high-depth cohort (default), 0 = low");
    app.add_flag("--header", header,         "Emit a TSV header line");
    app.add_option("--seed", seed,           "RNG seed (nonzero = deterministic)");

    CLI11_PARSE(app, argc, argv);

    if (seed != 0) setRandomSeed(seed);

    SimulationRequest req;
    req.nsnp       = nsnp;
    req.effectSize = effect;
    req.mafMin     = mafMin;
    req.mafMax     = mafMax;
    req.steps      = 1;
    req.collapse   = collapse;
    req.nboot      = nboot;
    req.useBootstrap = nboot > 1;
    req.stopEarly  = false;
    req.nthreads   = 1;

    if (family == "binomial")    req.family = Family::BINOMIAL;
    else if (family == "normal") req.family = Family::NORMAL;
    else { std::cerr << "unknown --family: " << family << "\n"; return 2; }

    if (stat == "common")   req.testStatistic = Statistic::COMMON;
    else if (stat == "cast") req.testStatistic = Statistic::CAST;
    else if (stat == "skat") req.testStatistic = Statistic::SKAT;
    else { std::cerr << "unknown --stat: " << stat << "\n"; return 2; }

    if (req.family == Family::BINOMIAL) {
        SimulationRequestGroup caseG;
        caseG.index = 0;
        caseG.n = ncase;
        caseG.family = Family::BINOMIAL;
        caseG.isCase = true;
        caseG.meanDepth = caseDepth;
        caseG.sdDepth = caseDepthSd;
        caseG.errorRate = caseErr;
        caseG.readDepth = caseHigh ? Depth::HIGH : Depth::LOW;
        req.groups.push_back(caseG);

        SimulationRequestGroup ctrlG;
        ctrlG.index = 1;
        ctrlG.n = ncontrol;
        ctrlG.family = Family::BINOMIAL;
        ctrlG.isCase = false;
        ctrlG.meanDepth = ctrlDepth;
        ctrlG.sdDepth = ctrlDepthSd;
        ctrlG.errorRate = ctrlErr;
        ctrlG.readDepth = ctrlHigh ? Depth::HIGH : Depth::LOW;
        req.groups.push_back(ctrlG);
    } else {
        SimulationRequestGroup g;
        g.index = 0;
        g.n = ncase + ncontrol;
        g.family = Family::NORMAL;
        g.isCase = false;
        g.normalMean = 0.0;
        g.normalSd = 1.0;
        g.meanDepth = caseDepth;
        g.sdDepth = caseDepthSd;
        g.errorRate = caseErr;
        g.readDepth = caseHigh ? Depth::HIGH : Depth::LOW;
        req.groups.push_back(g);
    }

    req.validate();

    Data result = startSimulation(req);

    if (header) {
        std::cout << "variant_idx\ttest_source\ttest_stat\tpvalue\n";
    }

    std::cout << std::setprecision(10);
    for (size_t i = 0; i < result.variants.size(); i++) {
        int np = result.variants[i].nPvals();
        for (int k = 0; k < np; k++) {
            // tests were pushed in order: trueGT, expectedGT, [calledGT]
            TestSettings& ts = result.tests[k];
            double p = result.variants[i].getPval(k);
            std::cout << i << "\t"
                      << gtSourceName(ts.getGenotype()) << "\t"
                      << statName(ts.getStatistic()) << "\t"
                      << p << "\n";
        }
    }
    return 0;
}
