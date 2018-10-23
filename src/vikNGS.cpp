#include "../src/vikNGS.h"
#include "./src/simulation/Simulation.h"
#include "Parser/Parser.h"
#include "Test/Test.h"
#include "Output/OutputHandler.h"

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <future>
#include <thread>

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

bool testBatch(SampleInfo* sampleInfo, std::vector<VariantSet*>& variants, std::vector<Test>& tests, int nboot){
    if(variants.size() <= 0)
        return true;

    for(VariantSet* vs : variants){
        for(size_t i=0; i < tests.size(); i++){
            double pval = runTest(sampleInfo, vs, tests[i], nboot, (nboot > 1));
            vs->addPval(pval);
        }
    }

    return true;
}

class ParallelTest
{
private:
    SampleInfo sampleInfo;
    std::vector<Test> tests;

    bool running;
    std::future<bool> testingDone;
    std::vector<VariantSet*> pointers;
    int nboot;

public:
    ParallelTest(SampleInfo si) : sampleInfo(si) { running = false;}
    void updateSampleInfoY(VectorXd Y){ this->sampleInfo.setY(Y); }
    void beginAssociationTest(std::vector<VariantSet*>& variants, Test& test, int nboot){
        std::vector<Test> v; v.push_back(test);
        beginAssociationTest(variants, v, nboot);
    }
    void beginAssociationTest(std::vector<VariantSet*>& variants, std::vector<Test>& tests, int nboot){
        running = true;
        this->pointers = variants;
        this->tests = tests;
        this->nboot = nboot;
        testingDone = std::async(std::launch::async,
                [this] { return testBatch(&sampleInfo, pointers, this->tests, this->nboot); }) ;
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

    if(simReq.family == Family::NORMAL)
        return startQuantitativeSimulation(simReq);

    printInfo("Simulating case-control data.");

    Data result;
    result.sampleInfo = simulateSampleInfo(simReq);

    if(STOP_RUNNING_THREAD)
        return result;

    result.variants = simulateVariant(simReq);

    if(STOP_RUNNING_THREAD)
        return result;

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
    tests.push_back(expectedGT);
    if(simReq.underNull())
        tests.push_back(calledGT);

    size_t nthreads = simReq.nthreads;
    std::vector<ParallelTest> threads;
    Test null(Genotype::NONE, Statistic::NONE);
    for(size_t i = 0; i < nthreads; i++)
        threads.emplace_back(result.sampleInfo);

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

                    if(STOP_RUNNING_THREAD)
                        break;

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

Data startQuantitativeSimulation(SimulationRequest& simReq) {

    printInfo("Simulating quantitative data.");

    Data result;
    result.sampleInfo = simulateSampleInfo(simReq);

    if(STOP_RUNNING_THREAD)
        return result;

    result.variants = simulateVariants(simReq);

    if(STOP_RUNNING_THREAD)
        return result;

    MatrixXd Y;
    if(simReq.family == Family::NORMAL)
       Y = addEffectOnY(simReq, result.variants);

    //VectorXd Y;
    //if(simReq.family == Family::NORMAL)
    //    Y = addEffectOnY(simReq, result.variants.at(0));

    if(STOP_RUNNING_THREAD)
        return result;

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
    Test rvsFalse(Genotype::EXPECTED, simReq.testStatistic);
    rvsFalse.setRVSFalse();

    tests.push_back(trueGT);
    tests.push_back(expectedGT);
    //if(simReq.underNull())
        tests.push_back(calledGT);
    //tests.push_back(rvsFalse);

    size_t nthreads = simReq.nthreads;
    std::vector<ParallelTest> threads;
    Test null(Genotype::NONE, Statistic::NONE);
    for(size_t i = 0; i < nthreads; i++)
        threads.emplace_back(result.sampleInfo);

    for(int i = 0; i < simReq.steps; i++){

        printInfo("Running step " + std::to_string(i+1) + " of " + std::to_string(simReq.steps) + ".");

        for(size_t j = 0; j < tests.size(); j++){
            tests[j].setSampleSize(simReq.nsamp(i));
            result.tests.push_back(tests[j]);
        }

        if(nthreads <= 1){
            for(size_t k = 0; k < result.variants.size(); k++){

                if(STOP_RUNNING_THREAD)
                    return result;

                result.sampleInfo.setY(Y.col(k));
                result.sampleInfo.setFamily(Family::NORMAL);

                for(size_t j = 0; j < tests.size(); j++){
                    double pval = runTest(&result.sampleInfo, &result.variants[k], tests[j], nboot, simReq.stopEarly);
                    result.variants[k].addPval(pval);
                }
            }
        }

        //multithreaded
        else{
            size_t counter = 0;
            std::vector<VariantSet*> vs;

            while(counter < result.variants.size()){

                if(STOP_RUNNING_THREAD)
                    break;

                for(size_t k = 0; k < nthreads; k++){
                    if(threads[k].isTestingDone())
                        threads[k].setDone();

                    if(!threads[k].isRunning()){
                        threads[k].updateSampleInfoY(Y.col(counter));
                        vs.push_back(&result.variants[counter]);
                        threads[k].beginAssociationTest(vs, tests, nboot);
                        vs.clear();
                        counter++;
                        if(counter >= result.variants.size())
                            break;
                    }
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }

            while(true){

                bool isDone = true;
                for(size_t k = 0; k < nthreads; k++){
                    if(threads[k].isTestingDone())
                        threads[k].setDone();
                    if(threads[k].isRunning())
                        isDone = false;
                }

                if(isDone)
                    break;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        }
    }

    return result;
}

