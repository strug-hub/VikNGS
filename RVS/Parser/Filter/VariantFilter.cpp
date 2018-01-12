#include "../InputParser.h"
#include "VariantFilterUtils.h"

/*
Filters variants read from VCF file based on missing sample data, failing to have 'PASS' in FILTER column
and whether the variant is an insertion or deletion (indel).

@param req Request object containing filtering parameters
@param variants Lines read from VCF file.
@param G Vector of group ID.
@return Filtered vector of VCFLines.
*/
std::vector<VCFLine> filterVariants(Request req, std::vector<VCFLine> &variants, VectorXd &G) {
	
	int before = variants.size();

	if (req.mustPASS)
		variants = filterTest(variants);

	if (req.onlySNPs)
		variants = snpTest(variants);

	if (req.useCommon())
		variants = filterMinorAlleleFrequency(variants, req.mafCutoff, true);
	else
		variants = filterMinorAlleleFrequency(variants, req.mafCutoff, false);

	//assumes number of groups = highest group ID
	int ngroup = G.maxCoeff() + 1;

	variants = missingTest(variants, G, ngroup, req.missingThreshold);

	std::sort(variants.begin(), variants.end(), lineCompare);
	variants = removeDuplicates(variants);

	variants = filterHomozygousVariants(variants);
	variants = filterMinorAlleleFrequency(variants, req.mafCutoff, req.useCommon());
	
	int nfiltered = variants.size();
	if (nfiltered > 0) {
		std::string s = nfiltered > 1 ? "s" : "";
		printInfo(std::to_string(nfiltered) + " variant" + s + " left after filtering.");
	}

	return variants;
}