#include "Request.h"

#include <fstream>

static const std::string ERROR_SOURCE = "REQUEST_BUILDER";

Request getDefaultRequest() {

	Request r;
    r.setAsSimulation(false);
    r.setInputFiles("", "");
    r.setOutputDir(".");

    r.setCollapse(5);
    r.setStopEarly(false);
    r.setNumberThreads(1);
    r.setBatchSize(1000);
    r.setKeepFiltered(false);

    r.setMustPASS(true);
    r.setOnlySNPs(true);
    r.setHighLowCutOff(30);
    r.setMafCutOff(0.05);
    r.setMissingThreshold(0.1);

    r.setMaxPos(-1);
    r.setMinPos(-1);
    r.setChromosomeFilter("");

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

bool Request::validate() {

    if(!this->simulation && this->vcfDir.size() <= 0)
        throwError(ERROR_SOURCE, "No VCF file provided.");
    if(!this->simulation && this->sampleDir.size() <= 0)
        throwError(ERROR_SOURCE, "No sample info file provided.", sampleDir);
    if(this->tests.size() <= 0)
        throwError(ERROR_SOURCE, "No association tests specified.");

    if (!checkFileExists(vcfDir))
        throwError(ERROR_SOURCE, "Cannot find file at VCF directory.", vcfDir);
    if (!checkFileExists(sampleDir))
        throwError(ERROR_SOURCE, "Cannot find file at sample info directory.", sampleDir);
    if (bedDir.size() > 0 && !checkFileExists(bedDir))
        throwError(ERROR_SOURCE, "Cannot find file at BED directory.", bedDir);

    //todo : make sure this function works for directories and not just files
    if (!checkFileExists(outputDir))
        throwError(ERROR_SOURCE, "Output directory is invalid.", outputDir);


    for(size_t i = 0; i < tests.size(); i++){
        if (nboot < 1 && tests[i].isRareTest() )
            throwError(ERROR_SOURCE, "Requested a bootstrap test but number of permutions is less than 1.", std::to_string(nboot));

        if (shouldCollapseK() && tests[i].isRareTest() && collapseSize < 2)
            throwError(ERROR_SOURCE, "Please set variant collapse value to a number greater than 1.", std::to_string(collapseSize));
    }

    if (nthreads < 1)
        throwError(ERROR_SOURCE, "Number of threads should be greater than 0.", std::to_string(nthreads));
    if(batchSize < 1)
        throwError(ERROR_SOURCE, "Batch size should be greater than 0.", std::to_string(batchSize));
    if (highLowCutOff < 1)
        throwError(ERROR_SOURCE, "High-low read depth threshold should be greater than 0.", std::to_string(highLowCutOff));
    if (mafCutoff < 0 || mafCutoff > 0.5)
        throwError(ERROR_SOURCE, "Minor allele frequency threshold should be a value between 0 and 0.5.", std::to_string(mafCutoff));
    if (missingThreshold < 0 || missingThreshold > 1)
        throwError(ERROR_SOURCE, "Missing threshold should be a value between 0 and 1.", std::to_string(missingThreshold));

    return true;
}
