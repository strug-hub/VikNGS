#pragma once
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

@return True if no issues detected during parsing or filtering.
@effect Modifies output params using information from input files.
*/
std::vector<VCFLine> parseAndFilter(Request req, MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd &P,
	std::vector<std::vector<int>> & interval) {
	
	std::map<std::string, int> IDmap = getSampleIDMap(req.vcfDir);
	parseSampleLines(req, IDmap, Y, Z, G, readGroup);

	std::vector<VCFLine> variants = parseVCFLines(req.vcfDir);
	variants = calculateExpectedGenotypes(variants);

	if (req.shouldCollapse())
		interval = parseBEDLines(req.bedDir, variants, req.shouldCollapseCoding(), req.shouldCollapseExon());

	variants = filterVariants(req, variants, G);

	if (variants.size() <= 0) 
		throwError(INPUT_PARSER, "No variants left after filtering step. No results to display.");
	

	MatrixXd x(variants[0].likelihood.size(), variants.size());
	for (int i = 0; i < variants.size(); i++)
		x.col(i) = variants[i].expectedGenotype;

	MatrixXd p(variants.size(), 3);
	for (int i = 0; i < variants.size(); i++)
		p.row(i) = variants[i].P;

	X = x;
	P = p;
	
	return variants;
}