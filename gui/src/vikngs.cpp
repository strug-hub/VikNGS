#include "../src/RVS.h"
#include "simulation/simulation.h"

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

    if (req.useCommon()) {

        std::vector<Variant> variants = runCommonTest(req, input);

        for (size_t i = 0; i < variants.size(); i++) {
            std::cout << variants[i].toString() << "\t" << variants[i].pvalue;
            std::cout << '\n';
        }

        return variants;
    }

  //if (!req.useCommon())

    std::vector<Variant> variants = runRareTest(req, input);

    for (size_t i = 0; i < variants.size(); i++) {
      std::cout << variants[i].toString() << "\t" << variants[i].pvalue;
      std::cout << '\n';
    }

    return variants;
}

std::vector<std::vector<Variant>> startSimulation(std::vector<SimulationRequest> simReqs) {

    std::vector<TestInput> inputs = simulate(simReqs);
    std::vector<std::vector<Variant>> results;
    printInfo("Starting tests...");

    for(int i = 0; i< simReqs.size(); i++){
        printInfo("Running step " + std::to_string(i) + " of " + std::to_string(simReqs.size()) + ".");

        SimulationRequest simReq = simReqs[i];
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

        setRegularTest(simReq.regular);
        setRVS(simReq.rvs);

        //todo
        setNumberThreads(1);

        Request req = getRequest();

        if (req.useCommon())
            results.push_back(runCommonTest(req, input));
        else
            results.push_back(runRareTest(req, input));
    }

    return results;
}

