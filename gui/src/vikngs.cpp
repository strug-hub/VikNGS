#include "../src/RVS.h"
#include "simulation/simulation.h"
#include "global.h"

#include <iostream>  
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

std::vector<Variant> startVikNGS(Request req) {

    createFile(req.outputDir);

    printInfo("Parsing files...");
    TestInput input = parseInfo(req);

    if(input.hasCovariates())
        printInfo(std::to_string(input.countCovariates()) + " covariates parsed");

    std::vector<Variant> variants = processVCF(input, req);
    return variants;
}

std::vector<Variant> runTest(TestInput &input, Request &req){

    //common

    if (req.useCommon()) {
        std::vector<Variant> variants = runCommonTest(req, input);
        return variants;
    }

    //rare

    std::vector<Variant> variants = runRareTest(req, input);
    return variants;
}

std::vector<std::vector<Variant>> startSimulation(SimulationRequest& simReq) {

    std::vector<TestInput> inputs = simulate(simReq);
    std::vector<std::vector<Variant>> results;

    printInfo("Starting tests...");

    for(int i = 0; i< simReq.steps; i++){
        printInfo("Running step " + std::to_string(i+1) + " of " + std::to_string(simReq.steps) + ".");

        TestInput input = inputs[i];

        initializeRequest();

        if(simReq.test=="calpha")
          useRareTest(simReq.test);
        else if(simReq.test=="cast")
          useRareTest(simReq.test);
        else
          useCommonTest();

        if(simReq.nboot > 1){
            useBootstrap(simReq.nboot);
            setStopEarly(simReq.stopEarly);
        }

        //todo
        setNumberThreads(1);

        Request req = getRequest();

        req.regularTest = true;
        req.useTrueGenotypes = true;
        req.rvs = false;

        if (req.useCommon())
            input.variants = runCommonTest(req, input);
        else
            input.variants = runRareTest(req, input);

        req.useTrueGenotypes = false;

        if (req.useCommon())
            input.variants = runCommonTest(req, input);
        else
            input.variants = runRareTest(req, input);

        req.regularTest = false;
        req.rvs = true;

        if (req.useCommon())
            input.variants = runCommonTest(req, input);
        else
            input.variants = runRareTest(req, input);

        if(STOP_RUNNING_THREAD)
            return results;

        results.push_back(input.variants);
    }

    return results;
}

