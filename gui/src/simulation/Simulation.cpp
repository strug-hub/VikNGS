#include "../src/Test/Test.h"
#include "Simulation.h"
#include <future>
#include <thread>
#include <chrono>

bool testBatch(SampleInfo* sampleInfo, MatrixXd& Y, std::vector<VariantSet*>& variants, Test& test, int nboot){
    if(variants.size() <= 0)
        return true;

    sampleInfo->setY(Y.col(0));
    int counter = 0;

    for(VariantSet* vs : variants){

        if(Y.cols() > 1){
            sampleInfo->setY(Y.col(counter));
            counter++;
        }

        double pval = runTest(sampleInfo, vs, test, nboot, (nboot > 1));
        vs->addPval(pval);
    }

    return true;
}


class ParallelTest
{
private:
    SampleInfo sampleInfo;
    std::vector<Test> tests;
    MatrixXd* Y;
    MatrixXd Ysubset;
    Test t = Test(Genotype::NONE, Statistic::NONE);

    bool running;
    std::future<bool> testingDone;
    std::vector<VariantSet*> pointers;
    int nboot;

public:
    ParallelTest(SampleInfo si, MatrixXd* phenotypes) : sampleInfo(si), Y(phenotypes) { running = false;}
    void updateSampleInfoY(VectorXd Y){ this->sampleInfo.setY(Y); }
    void beginAssociationTest(int startIdx, std::vector<VariantSet*>& variants, Test& test, int nboot){
        running = true;
        this->pointers = variants;
        this->nboot = nboot;
        t = test;

        MatrixXd tempY = *Y;
        if(Y->cols() > 1)
            Ysubset = Y->block(0, startIdx, Y->rows(), variants.size());
        else
            Ysubset = Y->col(0);

        testingDone = std::async(std::launch::async,
                [this] { return testBatch(&sampleInfo, Ysubset, pointers, t, this->nboot); }) ;
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
    // Record start time
    auto startTime = std::chrono::high_resolution_clock::now();

    Data result;
    Family fam = simReq.family;

    if(fam == Family::NORMAL)
        printInfo("Simulating quantitative data.");
    else if (fam == Family::BINOMIAL)
        printInfo("Simulating case-control data.");

    if(STOP_RUNNING_THREAD)
        return result;

    SampleInfo info;

    info.setFamily(fam);
    info.setG(simReq.getGroupVector());
    info.setGroupDepthMap(simReq.getReadDepthMap());

    if(STOP_RUNNING_THREAD)
        return result;

    printInfo("Generating genotype data.");

    VectorXd mafs = generateMafs(simReq);

    MatrixXd X;
    if(simReq.family == Family::BINOMIAL)
        X = simulateXCaseControl(simReq, mafs);
    else if(simReq.family == Family::NORMAL)
        X = simulateXNormal(simReq, mafs);

    printInfo("Simulating sequencing experiment.");

    std::vector<std::vector<Vector3d>> likelihoods = simulateSequencing(simReq, X);

    for(int i = 0; i < X.cols(); i += simReq.collapse){

        if(STOP_RUNNING_THREAD){
            Data empty;
            return empty;
        }

        VariantSet vs;
        result.variants.push_back(vs);

        for (int j = 0; j < simReq.collapse; j++){

            Variant v = randomVariant();
            VectorXd x = X.col(i+j);
            v.setTrueGenotypes(x);
            v.setExpectedGenotypes(likelihoods[i+j]);
            v.setCallGenotypes(likelihoods[i+j]);
            result.variants.back().addVariant(v);
        }
    }

    printInfo("Simulating phenotypes.");

    MatrixXd Y;
    if(simReq.family == Family::BINOMIAL)
        Y = simulateYCaseControl(simReq);
    else if(simReq.family == Family::NORMAL)
        Y = simulateYNormal(simReq, X, mafs);
    info.setY(Y.col(0));


    if(STOP_RUNNING_THREAD){ Data empty; return empty; }

    printInfo("Starting tests...");

    std::vector<Test> tests;
    Test trueGT(Genotype::TRUE, simReq.testStatistic);
    Test expectedGT(Genotype::EXPECTED, simReq.testStatistic);
    Test calledGT(Genotype::CALL, simReq.testStatistic);

    tests.push_back(trueGT);
    tests.push_back(expectedGT);
    if(simReq.underNull() || fam == Family::NORMAL)
        tests.push_back(calledGT);

    int nboot = 0;
    if(simReq.useBootstrap)
        nboot = simReq.nboot;
    if(!simReq.useBootstrap)
        nboot = 1;

    size_t nthreads = simReq.nthreads;
    std::vector<ParallelTest> threads;
    Test null(Genotype::NONE, Statistic::NONE);
    for(size_t i = 0; i < nthreads; i++)
        threads.emplace_back(info, &Y);

    result.sampleInfo = info;

    // Record processing time
    auto finishTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finishTime - startTime;
    result.processingTime = elapsed.count();
    startTime = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < simReq.steps; i++){

        printInfo("Running step " + std::to_string(i+1) + " of " + std::to_string(simReq.steps) + ".");

        for(size_t j = 0; j < tests.size(); j++){
            tests[j].setSampleSize(simReq.nsamp(i));
            result.tests.push_back(tests[j]);

            //single thread
            if(nthreads <= 1){
                for(size_t k = 0; k < result.variants.size(); k++){
                    if(Y.cols() > 1)
                        result.sampleInfo.setY(Y.col(k));
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
                int startIdx = 0;
                while(counter < result.variants.size()){

                    if(STOP_RUNNING_THREAD){ Data empty; return empty; }

                    vs.push_back(&result.variants[counter]);
                    counter++;

                    if(vs.size() >= batchSize && threadsRunning < (nthreads-1)){
                        threads[threadsRunning].beginAssociationTest(startIdx, vs, tests[j], nboot);
                        threadsRunning++;
                        vs.clear();
                        startIdx = counter;
                    }
                }

                if(vs.size() > 0)
                    threads[threadsRunning].beginAssociationTest(startIdx, vs, tests[j], nboot);

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

    // Record evaluation time
    finishTime = std::chrono::high_resolution_clock::now();
    elapsed = finishTime - startTime;
    result.evaluationTime = elapsed.count();

    return result;
}
