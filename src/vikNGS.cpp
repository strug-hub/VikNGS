#include "../src/vikNGS.h"
#include "Parser/Parser.h"
#include "Test/Test.h"
#include "Output/OutputHandler.h"
#include "Log.h"

#include <chrono>

Data startVikNGS(Request req) {

    printInfo("Starting vikNGS...");

    initializeOutputFiles(req.getOutputDir(), req.getRequestName());

    printInfo("Parsing files...");

    Data result;

    result.tests = req.getTests();
    result.sampleInfo = parseSampleInfo(req);

    if(req.shouldCollapseBed()){
        result.intervals = parseBEDLines(req.getBEDDir(), req.getCollapseType());
        req.setIntervals(&result.intervals);
    }

    auto startTime = std::chrono::high_resolution_clock::now();

    result.variants = processVCF(req, result.sampleInfo, result.variantsParsed);

    auto finishTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finishTime - startTime;
    result.evaluationTime = elapsed.count();
    result.processingTime = 0;

    outputPvals(result.variants, req.getOutputDir(), result.tests, req.getRequestName());
    printInfo("Results written to " + req.getOutputDir());
    return result;
}




