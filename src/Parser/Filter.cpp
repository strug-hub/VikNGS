#include "Filter.h"
#include "../Request.h"
#include "../Variant.h"

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
Filter filterByGenotypes(Request* req, Variant& variant, VectorXd& Y, Family family) {

    std::vector<GenotypeSource> genotypes = variant.getAllGenotypes();
    //must pass filter for all genotype sets

    for (GenotypeSource gt : genotypes){
        Vector3d* P = variant.getP(gt);
        if(!mafTest(P, req->getMAFCutOff(), req->useCommon()))
            return Filter::MAF;

        VectorXd* X = variant.getGenotype(gt);
        if(!missingTest(X, Y, req->getMissingThreshold(), family))
            return Filter::MISSING_DATA;

        if(!checkVariability(X))
            return Filter::NO_VARIATION;
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
bool missingTestCaseControl(VectorXd* X, VectorXd &Y, double missingThreshold){

    int nsamp = X->rows();

    double ncase = Y.sum();
    double ncontrol = nsamp - ncase;

    double missingCase = 0;
    double missingControl = 0;

    for (int i = 0; i < nsamp; i++) {

        if (std::isnan(X->coeff(i))){
            if(Y[i] < 1e-4)
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
bool missingTestQuantitative(VectorXd* X, double missingThreshold){

    int nsamp = X->rows();
    double nmissing = 0;

    for (int i = 0; i < nsamp; i++)
        if (std::isnan(X->coeff(i)))
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
bool mafTest(Vector3d* P, double mafCutoff, bool keepCommon) {

    double maf = 0.5 * P->coeff(1) + P->coeff(2);

    if (maf > 0.5)
        maf = 1 - maf;

    if ((keepCommon && maf > mafCutoff) || (!keepCommon && (maf < mafCutoff)))
        return true;
    else
        return false;

}

bool checkVariability(VectorXd* X){

    double mean = 0;
    double n = 0;
    for(int i = 0; i < X->rows(); i++)
        if(!std::isnan((*X)[i])){
            mean += (*X)[i];
            n++;
        }

    if(n < 2)
        return false;

    mean = mean/n;
    double var = 0;

    for(int i = 0; i < X->rows(); i++)
        if(!std::isnan((*X)[i]))
            var += ((*X)[i] - mean) * ((*X)[i] - mean);

    var = var/(n-1);
    return var > 1e-6;
}
