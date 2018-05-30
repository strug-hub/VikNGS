#include "InputParser.h"

static const std::string INPUT_PARSER = "input parser";

/*
Parses all input files and filters variants based on given criteria

Parameters contained in Request
------------------------
Parsing params:
@param vcfDir Directory of VCF file.
@param sampleDir Directory of sample directory file.
@param bedDir Directory of BED file.
@param highLowCutOff Cut-off for read depth. Anything >= value is considered high read depth.

Filtering params:
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@param onlySNPs If true, indels will be filtered based on REF and ALT columns.
@param mustPASS If true, variants where FILTER is not 'PASS' will be filtered.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param common Indicates whether to keep common or rare variants, (common = true, rare = false).
------------------------

Output params:
@param X Matrix of explanatory variable.
@param Y Vector of response variable.
@param Z Matrix of covariates.
@param G Vector of group ID.
@param P Vector with probability of 0, 1 or 2 minor alleles.
@param interval The indexes of variants within each interval.

@return TestInput object holding all the output parameters.
*/
TestInput parseAndFilter(Request req) {
	
	VectorXd Y, G; MatrixXd Z;
	std::map<int, int> readGroup;
	std::vector<std::vector<int>> interval;

    std::map<std::string, int> IDmap = getSampleIDMap(req.vcfDir);
	parseSampleLines(req, IDmap, Y, Z, G, readGroup);

    std::string family = determineFamily(Y);

    std::vector<Variant> variants = parseVCFLines(req, Y, family);

	if (variants.size() <= 0) 
		throwError(INPUT_PARSER, "No variants left after filtering step. No results to display.");
	
	if (req.shouldCollapseBed())
		interval = parseBEDLines(req.bedDir, variants, req.shouldCollapseCoding(), req.shouldCollapseExon());
	else if (!req.useCommon())
		interval = collapseEveryK(req.collapse, variants.size());

	MatrixXd X(variants[0].likelihood.size(), variants.size());
	for (int i = 0; i < variants.size(); i++)
		X.col(i) = variants[i].expectedGenotype;

	MatrixXd P(variants.size(), 3);
	for (int i = 0; i < variants.size(); i++)
		P.row(i) = variants[i].P;

    return buildTestInput(X, Y, Z, G, P, readGroup, interval, variants, family);
}
