#pragma once
#include "../InputParser.h"

/*
Verifies whether or not a string is a single nucleotide base (A, T, C or G).

@param base String to verify.
@return True if base is valid.
*/
inline bool validBase(std::string base) {
	return base == "T" || base == "A" || base == "C" || base == "G";
}

/*
Filters out non-SNP variants.

@param variant Variant read from VCF file.
@return True if variant is valid.
*/
inline bool snpTest(Variant &variant) {

    if (validBase(variant.alt) && validBase(variant.ref))
            return true;
    else
        return false;
}

/*
Removes variants that do not contain "PASS" in the FILTER column.

@param variant Variant read from VCF file.
@return True if variant is valid.
*/
inline bool filterTest(Variant &variant) {
    return variant.filter == "PASS";
}

/*
Filters homozygous variants (no variability across samples).

@param variant Variant read from VCF file.
@return True if variant is valid.
*/

inline bool homozygousTest(Variant &variant) {

    int j;

    VectorXd EG = variant.expectedGenotype;

    double mean = 0;
    double n = 0;
    for (j = 0; j < EG.size(); j++) {
        if (!std::isnan(EG[j])) {
            mean += EG[j];
            n++;
        }
    }
    mean = mean / n;

    double sd = 0;
    for (j = 0; j < EG.size(); j++)
        if (!std::isnan(EG[j]))
            sd += pow((EG[j] - mean), 2);

    if (1e-8*(n - 1) < sd)
        return true;
    else
        return false;
}


/*
Filters variants based on minor allele frequency. Removes rare or common variants depending on the value of keepCommon.

@param variant Lines read from VCF file.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param keepCommon Indicates whether to keep common or rare variants, (common = true, rare = false).
@return True if variant is valid.
*/

inline bool mafTest(Variant &variant, double mafCutoff, bool keepCommon) {

    double maf = 0.5 * variant.P[1] + variant.P[2];

    if (maf > 0.5)
        maf = 1 - maf;

    if ((keepCommon && maf > mafCutoff) || (!keepCommon && (maf < mafCutoff)))
        return true;
    else
        return false;
}

//case+control
//find percent missing for case and control separately
inline bool missingTestCaseControl(Variant &variant, VectorXd &Y, double missingThreshold){

    int nsamp = Y.rows();

    double ncase = Y.sum();
    double ncontrol = nsamp - ncase;

    double missingCase = 0;
    double missingControl = 0;

    for (int i = 0; i < nsamp; i++) {

        if (variant.likelihood[i].missing){
            if(Y[i] == 0)
                missingControl++;
            else if(Y[i] == 1)
                missingCase++;
        }
    }

    if(missingCase / ncase > missingThreshold)
        return false;
    if(missingControl / ncontrol > missingThreshold)
        return false;

    return true;
}

//quantitative
//find percent missing for all samples
inline bool missingTestQuantitative(Variant &variant, double missingThreshold){

    int nsamp = variant.likelihood.size();
    double nmissing= 0;

    for (int i = 0; i < nsamp; i++)
        if (variant.likelihood[i].missing)
            nmissing++;

    if (nmissing / (1.0*nsamp) > missingThreshold)
        return false;
    else
        return true;
}


/*
Counts the number of missing sample data for a variant and returns true if any group is missing more than missingThreshold.

@param variant VAriante read from VCF file.
@param Y Vector of response variable.
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@return True if variant is valid.
*/
inline bool missingTest(Variant &variant, VectorXd &Y, double missingThreshold, std::string family){

    if(family == "binomial")
        return missingTestCaseControl(variant, Y, missingThreshold);
    else
       return missingTestQuantitative(variant, missingThreshold);
}





