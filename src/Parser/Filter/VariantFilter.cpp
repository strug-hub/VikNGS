#include "../InputParser.h"
#include "VariantFilterUtils.h"

/*
Filters variants read from VCF file based on missing sample data, failing to have 'PASS' in FILTER column
and whether the variant is an insertion or deletion (indel).

@param req Request object containing filtering parameters
@param variant Single variant read from VCF file.
@param G Vector of group ID.
@param Y Vector of response variable.
@return Code specifiying whether or not variant passes. 0 is pass, anything else is fail.
*/
int filterVariant(Request &req, Variant &variant, VectorXd &Y, std::string family) {
	
    if (req.mustPASS && !filterTest(variant))
        return FILTER_FAIL;

    if (req.onlySNPs && !snpTest(variant))
        return SNP_FAIL;

    if(!mafTest(variant, req.mafCutoff, req.useCommon()))
        return MAF_FAIL;

    if(!missingTest(variant, Y, req.missingThreshold, family))
        return MISSING_FAIL;

    if(!homozygousTest(variant))
        return HOMOZYGOUS_FAIL;

    return PASS;
}

void printFilterResults(Request &req, std::vector<std::string> variantInfo, std::vector<int> failCode) {

    std::vector<std::string> codeMap(6);
    codeMap[PASS] = "PASS";
    codeMap[SNP_FAIL] = "Filtered for REF/ALT";
    codeMap[FILTER_FAIL] = "No PASS in FILTER";
    codeMap[MISSING_FAIL] = "Failed missing threshold";
    codeMap[HOMOZYGOUS_FAIL] = "No sample variability";
    codeMap[MAF_FAIL] = "Filtered for MAF";

    outputFiltered(variantInfo, failCode, codeMap, req.outputDir);
}


//    std::string s = nremoved > 1 ? "s" : "";

    //printInfo(std::to_string(nremoved) + " variant" + s + " were filtered by ALT and REF column (SNPs were retained).");
    //printInfo(std::to_string(nremoved) + " variant" + s + " were filtered by FILTER column (do not PASS).");
    //printInfo(std::to_string(nremoved) + " variant" + s + " filtered due to homozygous call in all samples.");

   /* if (keepCommon) {
        printInfo(std::to_string(nremoved) + " rare variant" + s +
            " filtered (minor allele frequency less than " + std::to_string(mafCutoff) + ").");
    }
    else {
        printInfo(std::to_string(nremoved) + " common variant" + s +
            " filtered (minor allele frequency greater than " + std::to_string(mafCutoff) + ").");




        printInfo(std::to_string(nremoved) + " variant" + s +
            " filtered by due to missing information (threshold = " + std::to_string(missingThreshold) + ").");
    }
    outputFiltered(removedVariants, "Filtered by missing data", outputDir);
    */
