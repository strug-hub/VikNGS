#include "../InputParser.h"
#include "VariantFilterUtils.h"

/*
Filters variants read from VCF file based on missing sample data, failing to have 'PASS' in FILTER column
and whether the variant is an insertion or deletion (indel).

@param req Request object containing filtering parameters
@param variants Lines read from VCF file.
@param G Vector of group ID.
@param Y Vector of response variable.
@return Filtered vector of Variants.
*/
std::vector<Variant> filterVariants(Request req, std::vector<Variant> &variants, VectorXd &Y, std::string family) {

	std::string outputDir = req.outputDir;
	
	if (req.mustPASS)
		variants = filterTest(variants, outputDir);

	if (req.onlySNPs)
		variants = snpTest(variants, outputDir);

    variants = filterMinorAlleleFrequency(variants, req.mafCutoff, req.useCommon(), outputDir);

    variants = missingTest(variants, Y, req.missingThreshold, family, outputDir);

    std::sort(variants.begin(), variants.end(), variantCompare);
	variants = removeDuplicates(variants, outputDir);

	variants = filterHomozygousVariants(variants, outputDir);
	variants = filterMinorAlleleFrequency(variants, req.mafCutoff, req.useCommon(), outputDir);
	
	int nfiltered = variants.size();
	if (nfiltered > 0) {
		std::string s = nfiltered > 1 ? "s" : "";
		printInfo(std::to_string(nfiltered) + " variant" + s + " left after filtering.");
	}

	return variants;
}
