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

    printInfo(std::to_string(nfiltered) + "/" + std::to_string(total) +
              percent(nfiltered, total) + "variants were filtered.");

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


