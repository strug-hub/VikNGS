#pragma once
#include "Enums.h"

struct IntervalSet;

#include <string>
#include <vector>

struct Request {
private:
    ///input files
    std::string vcfDir;
    std::string sampleDir;
    std::string bedDir;
    std::string outputDir;

    bool simulation;
    bool keepFiltered;
    bool makePlot;
    bool retainGt;

    std::vector<Test> tests;

    IntervalSet* intervals;
    CollapseType collapse;
    int collapseSize;

    int nboot;
    bool stopEarly;
    int nthreads;
    int batchSize;

    //filtering paramaters
    int highLowCutOff;
    double mafCutoff;
    double missingThreshold;
    bool onlySNPsFilter;
    bool mustPASSFilter;
    int minPos;
    int maxPos;
    //filterChrName = "" means no filter
    std::string filterChrName;

public:

    inline void setAsSimulation(bool value) { simulation = value; }

    inline void setInputFiles(std::string vcfDir, std::string sampleDir){
        this->vcfDir = vcfDir;
        this->sampleDir = sampleDir;
    }
    inline void setCollapseFile(std::string bedDir){ this->bedDir = bedDir; }
    inline void setOutputDir(std::string outputDir){ this->outputDir = outputDir; }

    inline void addTest(Test t) { this->tests.push_back(t); }

    inline void setCollapseGene(){ collapse = CollapseType::COLLAPSE_GENE; }
    inline void setCollapseExon(){ collapse = CollapseType::COLLAPSE_EXON; }
    inline void setCollapse(int k) { collapse = CollapseType::COLLAPSE_K; collapseSize = k; }

    inline void setBootstrap(int value) { nboot = value; }
    inline void setStopEarly(bool value) { stopEarly = value; }
    inline void setNumberThreads(int nthreads) { this->nthreads = nthreads; }
    inline void setBatchSize(int size) { this->batchSize = size; }
    inline void setKeepFiltered(bool value) { this->keepFiltered = value; }
    inline void setMakePlot(bool value) { this->makePlot = value; if(!value) setRetainGenotypes(false); }
    inline void setRetainGenotypes(bool value) { this->retainGt = value; }

    inline void setMustPASS(bool value) { mustPASSFilter = value; }
    inline void setOnlySNPs(bool value) { onlySNPsFilter = value; }
    inline void setHighLowCutOff(int highLowCutOff) { this->highLowCutOff = highLowCutOff; }
    inline void setMafCutOff(double mafCutoff) { this->mafCutoff = mafCutoff; }
    inline void setMissingThreshold(double missingThreshold) { this->missingThreshold = missingThreshold; }

    inline void setMinPos(int min) { this->minPos = min; }
    inline void setMaxPos(int max) { this->maxPos = max; }
    inline void setChromosomeFilter(std::string chrom) { this->filterChrName = chrom; }
    inline void setIntervals(IntervalSet *is) { this->intervals = is; }

    inline bool requireExpectedGenotypes(){
        bool result = false;
        for(Test t : tests)
            result = result || t.needExpectedGenotypes();
        return result;
    }
    inline bool requireGenotypeCalls(){
        bool result = false;
        for(Test t : tests)
            result = result || t.needGenotypeCalls();
        return result;
    }
    inline bool requireVCFCalls(){
        bool result = false;
        for(Test t : tests)
            result = result || t.needVCFCalls();
        return result;
    }
    inline bool useCommon(){
        bool result = false;
        for(Test t : tests)
            result = result || t.isCommonTest();
        return result;
    }
    inline std::vector<Test> getTests() { return tests; }

    inline bool shouldCollapseBed() {return this->collapse != CollapseType::COLLAPSE_K && bedDir.size() > 0; }
    inline bool shouldCollapseK() { return this->collapse == CollapseType::COLLAPSE_K; }
    inline bool shouldCollapseGene() { return this->collapse == CollapseType::COLLAPSE_GENE; }
    inline bool shouldCollapseExon() { return this->collapse == CollapseType::COLLAPSE_EXON; }

    inline std::string getVCFDir() { return vcfDir; }
    inline std::string getBEDDir() { return bedDir; }
    inline std::string getSampleDir() { return sampleDir; }
    inline std::string getOutputDir() { return outputDir; }
    inline CollapseType getCollapseType() { return collapse; }
    inline int getCollapseSize() { return collapseSize; }
    inline int getNumberThreads() { return nthreads; }
    inline int getBatchSize() { return this->batchSize; }
    inline bool shouldPlot() { return this->makePlot; }
    inline bool shouldRetainGenotypes() { return this->retainGt; }

    inline int getHighLowCutOff() { return highLowCutOff; }
    inline bool mustPASS() { return mustPASSFilter; }
    inline bool onlySNPs() { return onlySNPsFilter; }

    inline double getMissingThreshold() { return missingThreshold; }
    inline double getMAFCutOff() { return mafCutoff; }

    inline bool filterByChromosome() { return filterChrName.size() > 0; }
    inline bool filterByMinPosition() { return minPos >= 0; }
    inline bool filterByMaxPosition() { return maxPos >= 0; }
    inline std::string getFilterChromosome() { return filterChrName; }
    inline int getMinPosition() { return minPos; }
    inline int getMaxPosition() { return maxPos; }
    inline int bootstrapSize() { return nboot; }
    inline bool useBootstrap() { return nboot>0; }
    inline bool useStopEarly() { return stopEarly; }

    inline bool shouldKeepFiltered() { return keepFiltered; }
    inline IntervalSet* getIntervals() { return intervals; }
    bool validate();
};

//-------------------------------------------------------------------------------------


Request getDefaultRequest();


