#include "../src/vikNGS.h"
#include "simulation/simulation.h"
#include "Parser/InputParser.h"

#include <iostream>  
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

Data startVikNGS(Request req) {

    printInfo("Starting vikNGS...");

    initializeOutputFiles(req.outputDir);

    printInfo("Parsing files...");

    Data result;
    TestInput input = parseInfo(req);
    if(input.hasCovariates())
        printInfo(std::to_string(input.ncovariates()) + " covariates parsed");

    result.variants = processVCF(input, req);
    result.input = input;

    return result;
}

std::vector<Variant> runTest(TestInput &input, Request &req){

    if (req.useCommon()) {
        if(req.retainVariants){

        req.regularTest=false;
        req.rvs=true;
        input.variants = runCommonTest(req, input);

        req.regularTest=true;
        req.rvs=false;
        }
         return runCommonTest(req, input);

    }



    //rare
   if(req.retainVariants){
    req.regularTest=false;
    req.rvs=true;
    input.variants = runRareTest(req, input);
    req.regularTest=true;
   req.rvs=false;
   }
     return runRareTest(req, input);

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

