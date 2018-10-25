#include "../src/vikNGS.h"
#include "Parser/Parser.h"
#include "Test/Test.h"
#include "Output/OutputHandler.h"

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>


Data startVikNGS(Request req) {

    printInfo("Starting vikNGS...");

    initializeOutputFiles(req.getOutputDir());

    printInfo("Parsing files...");

    Data result;

    result.tests = req.getTests();
    result.sampleInfo = parseSampleInfo(req);

    if(req.shouldCollapseBed()){
        result.intervals = parseBEDLines(req.getBEDDir(), req.getCollapseType());
        req.setIntervals(&result.intervals);
    }

    result.variants = processVCF(req, result.sampleInfo);
    outputPvals(result.variants, req.getOutputDir(), result.tests);

    return result;
}




