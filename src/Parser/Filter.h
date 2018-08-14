#pragma once
#include "../vikNGS.h"

#include <vector>
#include <string>

inline bool validBase(std::string &base) {
    return base == "T" || base == "A" || base == "C" || base == "G";
}

Filter filterByVariantInfo(Request * req, std::string &chrom, std::string &pos, std::string &ref, std::string &alt, std::string &filter);
Filter filterByGenotypes(Request *req, Variant &variant, VectorXd &Y, Family family);

inline bool mafTest(Vector3d &P, double mafCutoff, bool keepCommon);
bool missingTestCaseControl(VectorXd &X, VectorXd &Y, double missingThreshold);
bool missingTestQuantitative(VectorXd &X, double missingThreshold);

inline bool missingTest(VectorXd &X, VectorXd &Y, double missingThreshold, Family family){
    if(family == Family::BINOMIAL)
        return missingTestCaseControl(X, Y, missingThreshold);
    else
       return missingTestQuantitative(X, missingThreshold);
}














/*
Filters homozygous variants (no variability across samples).

@param variant Variant read from VCF file.
@return True if variant is valid.
*/

inline bool homozygousTest(Variant &variant) {

    int j;

    double mean = 0;
    double n = 0;
    for (j = 0; j < variant.expectedGenotype.size(); j++) {
        if (!std::isnan(variant.expectedGenotype[j])) {
            mean += variant.expectedGenotype[j];
            n++;
        }
    }
    mean = mean / n;

    double sd = 0;
    for (j = 0; j < variant.expectedGenotype.size(); j++)
        if (!std::isnan(variant.expectedGenotype[j]))
            sd += pow((variant.expectedGenotype[j] - mean), 2);

    if (1e-8*(n - 1) < sd)
        return true;
    else
        return false;
}




