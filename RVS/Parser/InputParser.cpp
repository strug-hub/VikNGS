#pragma once
#include "InputParser.h"
#include "../Request.h"

static const std::string PARSER_ERROR = "input parser";


/*
Parses a VCF file.

@param vcfDir Directory of VCF file.
@return Vector with information about each line of VCF file.
*/
std::vector<VCFLine> parseVCF(std::string vcfDir) {
	std::vector<VCFLine> variants;
	
	if (!fileExists(vcfDir)) {
		printError("Could not find VCF file at given directory: " + vcfDir);
		throw std::runtime_error("File not found");
	}

	try {
		variants = parseVCFLines(vcfDir);
	}
	catch (...) {
		throw;
	}

	return variants;
}

/*
Parses a sample data file.

@param sampleDir Directory of sample directory file.
@param vcfDir Directory of VCF file.
@param Y Vector of response variable.
@param Z Matrix of covariates.
@param G Vector of group ID.
@param readGroup Mapping of group ID to high or low read group.
@param highLowCutOff Cut-off for read depth. Anything >= value is considered high read depth.

@return None.
@effect Fills Y, Z, G and readGroup with information from sample file.
*/
void parseSample(Request req, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup) {

	if (!fileExists(req.sampleDir)) 
		throwError(PARSER_ERROR, "Cannot find file at sample info directory.", req.sampleDir);

	std::map<std::string, int> IDmap;

	try {
		IDmap = getSampleIDMap(req.vcfDir);
		if (IDmap.size() <= 0)
			throwError(PARSER_ERROR, "No sample IDs could be found in the VCF file header.");

		parseSampleLines(req.sampleDir, IDmap, Y, Z, G, readGroup, req.highLowCutOff);
	}
	catch (...) {
		throw;
	}
}

/*
Filters out variants based on given criteria.

@param variants Lines read from VCF file.
@param G Vector of group ID.
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@param onlySNPs If true, indels will be filtered based on REF and ALT columns.
@param mustPASS If true, variants where FILTER is not 'PASS' will be filtered.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param common Indicates whether to keep common or rare variants, (common = true, rare = false).

@return Filtered vector of VCFLines.
*/
std::vector<VCFLine> filter(Request req, std::vector<VCFLine> & variants, VectorXd & G) {
	
	variants = filterVariants(req, variants, G);

	std::sort(variants.begin(), variants.end(), lineCompare);
	variants = removeDuplicates(variants);

	variants = calculateExpectedGenotypes(variants);
	variants = filterHomozygousVariants(variants);
	variants = filterMinorAlleleFrequency(variants, req.mafCutoff, req.useCommon());

	int variantsLeft = variants.size();
	if (variantsLeft > 0) {
		std::string s = variantsLeft > 1 ? "s" : "";
		printInfo(std::to_string(variantsLeft) + " variant" + s + " left after filtering.");
	}

	return variants;
}

/*
Parses a BED file to collapse variants.

@param bedDir Directory of BED file.
@param variants Lines read from VCF file.
@param interval The indexes of variants within each interval.

@return None.
@effect Fills interval using information from BED file and variants.
*/
void parseBED(Request req, std::vector<VCFLine> & variants, std::vector<std::vector<int>> & interval) {
	try {
		if (!fileExists(req.bedDir))
			throwError(PARSER_ERROR, "Cannot find file at BED directory.", req.bedDir);
		else
			interval = parseBEDLines(req.bedDir, variants, req.shouldCollapseCoding(), req.shouldCollapseExon() );
	}
	catch (...) {
		printWarning("Failed to parse BED file, will not be used to collapse variants.");
	}
}

/*
Parses all input files and filters variants based on given criteria

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

Output params:
@param X Matrix of explanatory variable.
@param Y Vector of response variable.
@param Z Matrix of covariates.
@param G Vector of group ID.
@param P Vector with probability of 0, 1 or 2 minor alleles.
@param interval The indexes of variants within each interval.

@return True if no errors detected during parsing or filtering.
@effect Modifies output params using information from input files.
*/
bool parseAndFilter(Request req, MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd &P,
	std::vector<std::vector<int>> & interval) {

	try {
		std::vector<VCFLine> variants = parseVCF(req.vcfDir);
		parseSample(req, Y, Z, G, readGroup);

		if (req.shouldCollapse())
			parseBED(req, variants, interval);

		variants = filter(req, variants, G);

		if (variants.size() <= 0) {
			printWarning("No variants left after filtering!");
			return false;
		}

		MatrixXd x(variants[0].likelihood.size(), variants.size());
		for (int i = 0; i < variants.size(); i++)
			x.col(i) = variants[i].expectedGenotype;

		MatrixXd p(variants.size(), 3);
		for (int i = 0; i < variants.size(); i++)
			p.row(i) = variants[i].P;

		X = x;
		P = p;

	}
	catch (...) {
		return false;
	}
	
	return true;
}