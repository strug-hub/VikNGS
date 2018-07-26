#include "Request.h"

static const std::string REQUEST_BUILDER = "request builder";
static Request request;

Request setDefaultParameters() {

	Request r;
	r.highLowCutOff = 30;
    r.collapseType = COLLAPSE_K;
	r.mafCutoff = 0.05;
	r.missingThreshold = 0.1;
	r.onlySNPs = true;
	r.mustPASS = true;

	r.outputDir = ".";

	r.minPos = -1;
	r.maxPos = -1;
    request.filterChr = "";

	r.test = COMMON_TEST;
	r.useBootstrap = false;
	r.collapse = 5;
	r.nboot = 0;
	r.nthreads = 1;

	r.stopEarly = false;

	r.rvs = true;
    r.regularTest = false;
    r.useTrueGenotypes = false;

    r.batchSize = 1000;
    r.retainVariants = false;

	return r;
}

/**
Checks whether or not a file exists.

@param dir Directory to check.
@return True if file exists.
*/
bool checkFileExists(const std::string& dir) {

	std::ifstream f(dir.c_str());
	if (f.good())
		return true;

	return false;
}

void initializeRequest(std::string vcfDir, std::string sampleDir) {

	request = setDefaultParameters();

	if (checkFileExists(vcfDir))
		request.vcfDir = vcfDir;
	else
		throwError(REQUEST_BUILDER, "Cannot find file at VCF directory.", vcfDir);
	
	if (checkFileExists(sampleDir))
		request.sampleDir = sampleDir;
	else
		throwError(REQUEST_BUILDER, "Cannot find file at sample info directory.", sampleDir);
}

void initializeRequest() {
    request = setDefaultParameters();
}
void setCollapse(int k){
    if(k < 2)
        throwError(REQUEST_BUILDER, "Collapse k value should be greater than 1.",
            std::to_string(k));

	request.collapse = k;
}

void setCollapseFile(std::string bedDir) {
	if (checkFileExists(bedDir))
		request.bedDir = bedDir;
	else
		throwError(REQUEST_BUILDER, "Cannot find file at BED directory.", bedDir);
}

void setOutputDir(std::string outputDir) {

	//todo : make sure this function works for directories and not just files
	if (checkFileExists(outputDir))
		request.outputDir = outputDir;
	else
		throwError(REQUEST_BUILDER, "Output directory is invalid.", outputDir);
}

void setCollapseGene() {
	request.collapseType = COLLAPSE_GENE;
}
void setCollapseExon() {
	request.collapseType = COLLAPSE_EXON;
}
/*
void setCollapseCoding() {
	request.collapseType = COLLAPSE_CODING;
}*/

void setOnlySNPs(bool value) {
	request.onlySNPs = value;
}

void setRVS(bool value) {
    request.rvs = value;
    if(value)
        request.regularTest = false;
}

void setRegularTest(bool value){
    request.regularTest = value;
    if(value)
        request.rvs=false;
}

void setMustPASS(bool value) {
	request.mustPASS = value;
}
void setStopEarly(bool value){
	request.stopEarly = value;
}

void useBootstrap(int nboot) {
	request.useBootstrap = true;
	if(nboot < 1)
		throwError(REQUEST_BUILDER, "Number of bootstrapping permutions should be greater than 0.",
			std::to_string(nboot));

	request.nboot = nboot;
}

void useRareTest(std::string type) {
	request.test = RARE_TEST;
	request.rareTest = type;
}
void useCommonTest() {
	request.test = COMMON_TEST;
}

void setMinPos(int min){
    if (min < 0)
        throwError(REQUEST_BUILDER, "Filter from POS should be greater than or equal to 0.",
            std::to_string(min));

	request.minPos = min;
}

void setMaxPos(int max){
    if (max < 1)
        throwError(REQUEST_BUILDER, "Filter from POS should be greater than 0.",
            std::to_string(max));

	request.maxPos = max;
}

void setFilterChr(std::string chr){
    request.filterChr = chr;
}

void useTrueGenotypes(){
    request.regularTest = true;
    request.useTrueGenotypes = true;
}


void setNumberThreads(int nthreads) {
	if (nthreads < 1)
		throwError(REQUEST_BUILDER, "Number of threads should be greater than 0.",
			std::to_string(nthreads));

	request.nthreads = nthreads;
}

void setRetainVariants(bool value) {
    request.retainVariants = value;
}

void setBatchSize(int batch) {
    if(batch < 1)
        throwError(REQUEST_BUILDER, "Batch size should be greater than 0.",
            std::to_string(batch));
    request.batchSize = batch;
}

void setHighLowCutOff(int highLowCutOff) {
	if (highLowCutOff < 1)
		throwError(REQUEST_BUILDER, "High-low read depth threshold should be greater than 0.",
			std::to_string(highLowCutOff));

	request.highLowCutOff = highLowCutOff;
}

void setMAFCutoff(double mafCutoff) {
	if (mafCutoff < 0 || mafCutoff > 0.5)
		throwError(REQUEST_BUILDER, "Minor allele frequency threshold should be a value between 0 and 0.5.",
			std::to_string(mafCutoff));

	request.mafCutoff = mafCutoff;
}

void setMissingThreshold(double missingThreshold) {
    if (missingThreshold < 0 || missingThreshold > 1)
        throwError(REQUEST_BUILDER, "Missing threshold should be a value between 0 and 1.",
			std::to_string(missingThreshold));

	request.missingThreshold = missingThreshold;
}

Request getRequest() { return request; }
