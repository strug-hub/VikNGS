#include "../src/vikNGS.h"
#include "./src/simulation/Simulation.h"
#include "Parser/Parser.h"
#include "Test/Test.h"

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <future>
#include <thread>

Data startVikNGS(Request req) {

    printInfo("Starting vikNGS...");

   // initializeOutputFiles(req.getOutputDir());

    printInfo("Parsing files...");

    Data result;

    if(req.shouldCollapseBed()){
        result.intervals = parseBEDLines(req.getBEDDir(), req.getCollapseType());
        req.setIntervals(& result.intervals);
    }

    result.tests = req.getTests();
    result.sampleInfo = parseSampleInfo(req);
    processVCF(req, result.sampleInfo, &result.variants);

    return result;
}

bool testBatch(SampleInfo* sampleInfo, std::vector<VariantSet*>& variants, Test& test, int nboot){
    if(variants.size() <= 0)
        return true;

    for(VariantSet* vs : variants){
        double pval = runTest(sampleInfo, vs, test, nboot, (nboot > 1));
        vs->addPval(pval);
    }

    return true;
}

class ParallelTest
{
private:
    SampleInfo* sampleInfo;
    Test test;

    bool running;
    std::future<bool> testingDone;
    std::vector<VariantSet*> pointers;
    int nboot;

public:
    ParallelTest(SampleInfo *si, Test& t) : sampleInfo(si), test(t) { running = false;}

    void beginAssociationTest(std::vector<VariantSet*>& variants, Test& test, int nboot){
        running = true;
        this->pointers = variants;
        this->test = test;
        this->nboot = nboot;
        testingDone = std::async(std::launch::async,
                [this] { return testBatch(sampleInfo, pointers, this->test, this->nboot); }) ;
    }
    inline bool isTestingDone(){
        return running && testingDone.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }
    inline void setDone(){
        running = false;
        pointers.clear();
    }
    inline bool isRunning(){ return running; }

};


Data startSimulation(SimulationRequest& simReq) {

    printInfo("Simulating data.");

    Data result;
    result.sampleInfo = simulateSampleInfo(simReq);
    result.variants = simulateVariants(simReq);

    printInfo("Starting tests...");

    int nboot = 0;
    if(simReq.useBootstrap)
        nboot = simReq.nboot;
    if(!simReq.useBootstrap)
        nboot = 1;

    std::vector<Test> tests;
    Test trueGT(Genotype::TRUE, simReq.testStatistic);
    Test expectedGT(Genotype::EXPECTED, simReq.testStatistic);
    Test calledGT(Genotype::CALL, simReq.testStatistic);

    tests.push_back(trueGT);
    tests.push_back(calledGT);
    tests.push_back(expectedGT);

    size_t nthreads = simReq.nthreads;
    std::vector<ParallelTest> threads;
    Test null(Genotype::NONE, Statistic::NONE);
    for(size_t i = 0; i < nthreads; i++)
        threads.emplace_back(&result.sampleInfo, null);

    for(int i = 0; i < simReq.steps; i++){

        printInfo("Running step " + std::to_string(i+1) + " of " + std::to_string(simReq.steps) + ".");

        for(size_t j = 0; j < tests.size(); j++){
            tests[j].setSampleSize(simReq.nsamp(i));
            result.tests.push_back(tests[j]);

            if(nthreads <= 1){
                for(size_t k = 0; k < result.variants.size(); k++){
                    double pval = runTest(&result.sampleInfo, &result.variants[k], tests[j], nboot, simReq.stopEarly);
                    result.variants[k].addPval(pval);
                }
            }
            //multithreaded
            else{
                size_t batchSize = result.variants.size()/nthreads;
                size_t counter = 0;
                size_t threadsRunning = 0;
                std::vector<VariantSet*> vs;
                while(counter < result.variants.size()){

                    vs.push_back(&result.variants[counter]);

                    if(vs.size() >= batchSize && threadsRunning < (nthreads-1)){
                        threads[threadsRunning].beginAssociationTest(vs, tests[j], nboot);
                        threadsRunning++;
                        vs.clear();
                    }
                    counter++;
                }

                if(vs.size() > 0)
                    threads[threadsRunning].beginAssociationTest(vs, tests[j], nboot);

                while(true){

                    bool stop = true;
                    for(size_t t = 0; t < nthreads; t++){
                        if(threads[t].isRunning()){
                            if(!threads[t].isTestingDone())
                                stop = false;
                            else
                                threads[t].setDone();
                        }
                    }

                    if(stop)
                        break;
                    else
                        std::this_thread::sleep_for(std::chrono::milliseconds(100));
                }
            }
        }

    }

    return result;
}

