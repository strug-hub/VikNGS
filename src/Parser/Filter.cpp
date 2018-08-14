#include "Parser.h"
#include "Filter.h"

/*
Filters variants read from VCF file based variant specific columns of VCF file.

@param req Request object containing filtering parameters
@param chrom CHROM column value.
@param pos POS column value.
@param ref REF column value.
@param alt ALT column value.
@param filter FILTER column value.

@return Filter enum specifiying whether or not variant passes.
*/
Filter filterByVariantInfo(Request *req, std::string &chrom, std::string &pos, std::string &ref, std::string &alt, std::string &filter){

    int position = std::stoi(pos);
    if (req->filterByMinPosition() && position < req->getMinPosition())
        return Filter::IGNORE;
    if (req->filterByMaxPosition() && position > req->getMaxPosition())
        return Filter::IGNORE;
    if (req->filterByChromosome() && chrom != req->getFilterChromosome())
        return Filter::IGNORE;

    if (req->onlySNPs() && (!validBase(ref) || !validBase(alt)))
        return Filter::NOT_SNP;
    if (req->mustPASS() && filter != "PASS")
        return Filter::NO_PASS;

    return Filter::VALID;
}

/*
Filters variants read from VCF file based on genotype data. Filters each
genotype separately, fails if a single genotype fails.
Checks for missing data and filters data for common or rare test.

@param req Request object containing filtering parameters.
@param variant Single variant read from VCF file.
@param Y Vector of response variable.
@param family Specifies case-control or quantatitive for missing filter.

@return Filter enum based on whether or not variant passes the filters.
*/
Filter filterByGenotypes(Request *req, Variant &variant, VectorXd &Y, Family family) {

    std::vector<Genotype> genotypes;
    //must pass filter for all genotype sets

    for (Genotype gt : genotypes){
        Vector3d P = variant.getP(gt);
        if(!mafTest(P, req->getMAFCutOff(), req->useCommon()))
            return Filter::MAF;

        VectorXd X = variant.getGenotype(gt);
        if(!missingTest(X, Y, req->getMissingThreshold(), family))
            return Filter::MISSING_DATA;
    }

    return Filter::VALID;
}


/**
Counts the number of missing sample data for a variant and returns false if either
cases or controls is missing more than missingThreshold.

@param X Vector of genotypes.
@param Y Vector of response variable.
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@return True if variant is valid.
*/
bool missingTestCaseControl(VectorXd &X, VectorXd &Y, double missingThreshold){

    int nsamp = X.rows();

    double ncase = Y.sum();
    double ncontrol = nsamp - ncase;

    double missingCase = 0;
    double missingControl = 0;

    for (int i = 0; i < nsamp; i++) {

        if (std::isnan(X[i])){
            if(Y[i] < 1e4)
                missingControl++;
            else if(Y[i] > 0.9999)
                missingCase++;
        }
    }

    if(missingCase / ncase > missingThreshold)
        return false;
    if(missingControl / ncontrol > missingThreshold)
        return false;

    return true;
}

/**
Counts the number of missing sample data for a variant and returns false if genotypes
is missing more than missingThreshold.

@param X Vector of genotypes.
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@return True if variant is valid.
*/
bool missingTestQuantitative(VectorXd &X, double missingThreshold){

    int nsamp = X.rows();
    double nmissing= 0;

    for (int i = 0; i < nsamp; i++)
        if (std::isnan(X[i]))
            nmissing++;

    if (nmissing / (1.0*nsamp) > missingThreshold)
        return false;
    else
        return true;
}

/**

Filters variants based on minor allele frequency. Removes rare or common variants depending on the value of keepCommon.

@param P Genotype frequencies.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param keepCommon Indicates whether to keep common or rare variants, (common = true, rare = false).
@return True if variant is valid.
*/
inline bool mafTest(Vector3d &P, double mafCutoff, bool keepCommon) {

    double maf = 0.5 * P[1] + P[2];

    if (maf > 0.5)
        maf = 1 - maf;

    if ((keepCommon && maf > mafCutoff) || (!keepCommon && (maf < mafCutoff)))
        return true;
    else
        return false;

}






















std::string percent(int num, int denom){

    std::string pcnt = std::to_string(std::round(1000.0*num/denom)/10.0);
    std::string cut = "";
    for(int i =0; i<pcnt.size(); i++){
        cut += pcnt[i];
        if (pcnt[i] == '.' && pcnt.size() > (i+1)){
            cut += pcnt[i+1];
            break;
        }
    }

    return " (" + cut + "%) ";
}

void printFilterResults(Request &req, std::vector<std::string> variantInfo, std::vector<int> failCode, int total){

    int ncodes = 6;
    int nfiltered = variantInfo.size();

    std::vector<std::string> codeMap(ncodes);
    codeMap[PASS] = "PASS";
    codeMap[SNP_FAIL] = "Filtered for REF/ALT";
    codeMap[FILTER_FAIL] = "No PASS in FILTER";
    codeMap[MISSING_FAIL] = "Failed missing threshold";
    codeMap[HOMOZYGOUS_FAIL] = "No sample variability";
    codeMap[MAF_FAIL] = "Filtered for MAF";

    std::vector<int> codeCount(ncodes);

    printInfo(std::to_string(total - nfiltered) + "/" + std::to_string(total) +
              percent(total - nfiltered, total) + "variants remain after filtering.");

    for(int i = 0; i< failCode.size(); i++)
        codeCount[failCode[i]]++;

    for(int i = 0; i< codeMap.size(); i++){

        std::string s = codeCount[i] > 1 ? "s" : "";
        if(codeCount[i] < 1)
            continue;

        if(i == SNP_FAIL)
            printInfo(std::to_string(codeCount[i]) + " variant" + s + percent(codeCount[i],nfiltered) + "filtered by ALT and REF column (SNPs were retained).");
        if(i == FILTER_FAIL)
            printInfo(std::to_string(codeCount[i]) + " variant" + s + percent(codeCount[i],nfiltered) + "filtered by FILTER column (do not PASS).");
        if(i == HOMOZYGOUS_FAIL)
            printInfo(std::to_string(codeCount[i]) + " variant" + s + percent(codeCount[i],nfiltered) + "filtered due to homozygous call in samples (no variability).");
        if(i == MISSING_FAIL)
            printInfo(std::to_string(codeCount[i]) + " variant" + s + percent(codeCount[i],nfiltered) + "filtered by due to missing information (threshold = " + std::to_string(req.missingThreshold) + ").");
        if(i == MAF_FAIL){

            if (req.useCommon())
                 printInfo(std::to_string(codeCount[i]) + percent(codeCount[i],nfiltered) + "rare variant" + s +
                     " filtered (minor allele frequency less than " + std::to_string(req.mafCutoff) + ").");

             else
                 printInfo(std::to_string(codeCount[i]) + percent(codeCount[i],nfiltered) + "common variant" + s +
                     " filtered (minor allele frequency greater than " + std::to_string(req.mafCutoff) + ").");
        }

    }

    outputFiltered(variantInfo, failCode, codeMap, req.outputDir);
}

bool isIn(std::string &vcfLine, int minPos, int maxPos, std::string &chr){

    std::string vcfChr = "";
    std::string vcfPos = "";
    bool tab = false;

    for(int i = 0; i < vcfLine.size(); i++){

        if(vcfLine[i] == VCF_SEPARATOR){
            if(tab)
                break;

            tab = true;
            if (chr.size() > 0 && vcfChr != chr)
                return false;

            continue;
        }

        if(!tab)
            vcfChr += vcfLine[i];
        else
            vcfPos += vcfLine[i];
    }

    //possibility of throwing error
    int vcfPosValue = std::stoi(vcfPos);

    if((minPos > -1 && vcfPosValue < minPos) || (maxPos > -1 && vcfPosValue > maxPos))
        return false;

    return true;
}

int isIn(std::string vcfLine, int minPos, int maxPos, std::string &chr, std::vector<Interval> &intervals){

    std::string vcfChr = "";
    std::string vcfPos = "";
    bool tab = false;

    for(int i = 0; i < vcfLine.size(); i++){

        if(vcfLine[i] == VCF_SEPARATOR){
            if(tab)
                break;

            tab = true;
            if (chr.size() > 0 && vcfChr != chr)
                return false;

            continue;
        }

        if(!tab)
            vcfChr += vcfLine[i];
        else
            vcfPos += vcfLine[i];
    }

    //possibility of throwing error
    int vcfPosValue = std::stoi(vcfPos);

    if(vcfPosValue < minPos || vcfPosValue > maxPos)
        return false;

    return true;
}




int findInterval(std::vector<Interval> &intervals, Variant &variant){

    //todo: optimize!@!!!!!
    for(int i = 0; i < intervals.size(); i++){

        if(intervals[i].start <= variant.pos && intervals[i].end >= variant.pos
                && intervals[i].chr == variant.chr)
            return i;
    }
    return -1;
}




