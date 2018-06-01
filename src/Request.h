#pragma once
#include "Log.h"
#include <string>
#include <fstream>

static const int COMMON_TEST = 1;
static const int RARE_TEST = 2;
static const int COLLAPSE_K = 0;
static const int COLLAPSE_GENE = 1;
static const int COLLAPSE_EXON = 2;
//static const int COLLAPSE_CODING = 3;

struct Request {

	///input files
	std::string vcfDir;
	std::string sampleDir;
	std::string bedDir = "";

	int highLowCutOff;
	int collapseType;

	//filtering paramaters
	double mafCutoff;
	double missingThreshold;
	bool onlySNPs;
	bool mustPASS;
	int minPos;
	int maxPos;

	std::string outputDir;

	int test;
	bool useBootstrap;
	int nboot;
	std::string rareTest;
	int collapse;	
	bool stopEarly;
    int nthreads;

	bool rvs;
    bool regularTest;

    int batchSize;
    bool retainVariants;

	inline bool useCommon() { return test == COMMON_TEST; }
    inline bool shouldCollapseBed() { return collapseType != COLLAPSE_K && bedDir != ""; }
    inline bool shouldCollapseK() { return collapseType == COLLAPSE_K; }
	inline bool shouldCollapseGene() { return collapseType == COLLAPSE_GENE; }
    inline bool shouldCollapseExon() { return collapseType == COLLAPSE_EXON; }
    //inline bool shouldCollapseCoding() { return collapseType == COLLAPSE_CODING; }

};

/**
Creates a Request object whether or not a file exists.

@param vcfDir directory pointing to a VCF file
@param sampleDir directory pointing to a sample info file
@throws Exception if vcfDir or sampleDir are not valid directories
@effect Creates a new static Request object with default parameters
*/
void initializeRequest(std::string vcfDir, std::string sampleDir);
void initializeRequest();

/**
The following functions modify the parameters of the static Request object.

@throws Exception parameter value is invalid
@effect Modifies parameters in Request object.
-------------------------------------------------------------------------------------
*/
void setCollapseFile(std::string bedDir);

void setCollapseGene();
void setCollapseExon();
//void setCollapseCoding();
void setCollapse(int k);
void setOnlySNPs(bool value);
void setRegularTest(bool value);
void setMustPASS(bool value); 
void setStopEarly(bool value);
void useRareTest(std::string type);
void useCommonTest();
void useBootstrap(int nboot);
void setRVS(bool value);
void setRetainVariants(bool value);
void setBatchSize(int size);
void setOutputDir(std::string outputDir);
void setMinPos(int min);
void setMaxPos(int max);

void setNumberThreads(int nthreads);
void setHighLowCutOff(int highLowCutOff);
void setMAFCutoff(double mafCutoff);
void setMissingThreshold(double missingThreshold);

//-------------------------------------------------------------------------------------

/**
Returns the Request object with all modified parameters applied

@return Request object
*/
Request getRequest();
