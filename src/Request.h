#pragma once
#include "vikNGS.h"
#include <string>

struct Request {
private:
	///input files
	std::string vcfDir;
	std::string sampleDir;
    std::string bedDir;
    std::string outputDir;

    bool simulation;

    std::vector<Test> tests;
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
    bool onlySNPs;
    bool mustPASS;
    int minPos;
    int maxPos;
    //filterChr = "" means no filter
    std::string filterChr;

public:

    inline void setAsSimulation(bool value) { simulation = value; }

    inline void setInputFiles(std::string vcfDir, std::string sampleDir){
        this->vcfDir = vcfDir;
        this->sampleDir = sampleDir;
    }
    inline void setCollapseFile(std::string bedDir){ this->bedDir = bedDir; }
    inline void setOutputDir(std::string outputDir){ this->outputDir = outputDir; }

    inline void doCommonRVS() { tests.push_back(Test::COMMON_LIKELIHOOD_RVS); }
    inline void doRareSkatRVS(int nboot) { tests.push_back(Test::RARE_LIKELIHOOD_RVS_SKAT); this->nboot = nboot; }
    inline void doRareCalphaRVS(int nboot) { tests.push_back(Test::RARE_LIKELIHOOD_RVS_CALPHA); this->nboot = nboot; }
    inline void doCommonCalls() { tests.push_back(Test::COMMON_REGULAR_GTCALL); }
    inline void doRareSkatCalls(int nboot) { tests.push_back(Test::RARE_REGULAR_GTCALL_SKAT); this->nboot = nboot; }
    inline void doRareCalphaCalls(int nboot) { tests.push_back(Test::RARE_REGULAR_GTCALL_CALPHA); this->nboot = nboot; }
    inline void doCommonTrue() { tests.push_back(Test::COMMON_REGULAR_TRUE); }
    inline void doRareSkatTrue(int nboot) { tests.push_back(Test::RARE_REGULAR_TRUE_SKAT); this->nboot = nboot; }
    inline void doRareCalphaTrue(int nboot) { tests.push_back(Test::RARE_REGULAR_TRUE_CALPHA); this->nboot = nboot; }

    inline void setCollapseGene(){ collapse = CollapseType::COLLAPSE_GENE; }
    inline void setCollapseExon(){ collapse = CollapseType::COLLAPSE_EXON; }
    inline void setCollapse(int k) { collapse = CollapseType::COLLAPSE_K; collapseSize = k; }

    inline void setStopEarly(bool value) { stopEarly = value; }
    inline void setNumberThreads(int nthreads) { this->nthreads = nthreads; }
    inline void setBatchSize(int size) { this->batchSize = size; }

    inline void setMustPASS(bool value) { mustPASS = value; }
    inline void setOnlySNPs(bool value) { onlySNPs = value; }
    inline void setHighLowCutOff(double mafCutoff) { this->highLowCutOff = highLowCutOff; }
    inline void setMafCutOff(int highLowCutOff) { this->highLowCutOff = highLowCutOff; }
    inline void setMissingThreshold(double missingThreshold) { this->missingThreshold = missingThreshold; }

    inline void setMinPos(int min) { this->minPos = min; }
    inline void setMaxPos(int max) { this->maxPos = max; }
    inline void setChromosomeFilter(std::string chrom) { this->filterChr = chrom; }

    inline bool useBootstrap(int i) { return useRare(i); }

    inline bool useCommon(int i) { return test[i] == Test::COMMON_LIKELIHOOD_NORVS ||
                                 test[i] == Test::COMMON_LIKELIHOOD_RVS ||
                                 test[i] == Test::COMMON_REGULAR_GTCALL ||
                                 test[i] == Test::COMMON_REGULAR_TRUE; }

    inline bool useRare(int i) { return test[i] == Test::RARE_LIKELIHOOD_NORVS_CALPHA ||
                                 test[i] == Test::RARE_LIKELIHOOD_NORVS_SKAT ||
                                 test[i] == Test::RARE_LIKELIHOOD_RVS_CALPHA ||
                                 test[i] == Test::RARE_LIKELIHOOD_RVS_SKAT ||
                                 test[i] == Test::RARE_REGULAR_GTCALL_CALPHA ||
                                 test[i] == Test::RARE_REGULAR_GTCALL_SKAT ||
                                 test[i] == Test::RARE_REGULAR_TRUE_CALPHA ||
                                 test[i] == Test::RARE_REGULAR_TRUE_SKAT; }

    inline bool useRVS(int i) { return test[i] == Test::COMMON_LIKELIHOOD_RVS ||
                                 test[i] == Test::RARE_LIKELIHOOD_RVS_SKAT ||
                                 test[i] == Test::RARE_LIKELIHOOD_RVS_CALPHA; }

    inline bool useRegular(int i) { return test[i] == Test::COMMON_REGULAR_GTCALL ||
                                 test[i] == Test::COMMON_REGULAR_TRUE ||
                                 test[i] == Test::RARE_REGULAR_GTCALL_CALPHA ||
                                 test[i] == Test::RARE_REGULAR_GTCALL_SKAT ||
                                 test[i] == Test::RARE_REGULAR_TRUE_CALPHA ||
                                 test[i] == Test::RARE_REGULAR_TRUE_SKAT; }

    inline bool shouldCollapseBed() { return collapseType != CollapseType::COLLAPSE_K && bedDir != ""; }
    inline bool shouldCollapseK() { return collapseType == CollapseType::COLLAPSE_K; }
    inline bool shouldCollapseGene() { return collapseType == CollapseType::COLLAPSE_GENE; }
    inline bool shouldCollapseExon() { return collapseType == CollapseType::COLLAPSE_EXON; }

    bool validate();
};

//-------------------------------------------------------------------------------------


Request getDefaultRequest();


