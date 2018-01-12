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
Counts the number of missing sample data for a variant and returns true if any group is missing more than missingThreshold.

@param variants Lines read from VCF file.
@param G Vector of group ID.
@param ngroup Number of unique groups in G.
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@return Filtered vector of VCFLines, variants missing above threshold removed.
*/
std::vector<VCFLine> missingTest(std::vector<VCFLine> &variants, VectorXd &G, int ngroup, double missingThreshold);


/*
Removes variants that do not contain "PASS" in the FILTER column.

@param variants Lines read from VCF file.
@return Filtered vector of VCFLines.
*/
std::vector<VCFLine> filterTest(std::vector<VCFLine> &variants);


/*
Filters out non-SNP variants.

@param variants Lines read from VCF file.
@return Filtered vector of VCFLines, SNPs retained.
*/
std::vector<VCFLine> snpTest(std::vector<VCFLine> &variants);


/*
Filters variants based on minor allele frequency. Removes rare or common variants depending on the value of keepCommon.

@param variants Lines read from VCF file.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param keepCommon Indicates whether to keep common or rare variants, (common = true, rare = false).
@return Filtered vector of VCFLines (either common or rare variants removed).
*/
std::vector<VCFLine> filterMinorAlleleFrequency(std::vector<VCFLine> &variants, double mafCutoff, bool keepCommon);

/*
Filters duplicated variants (same position and chromosome).

@param variants Lines read from VCF file.
@requires variants to be sorted in order of chr, loc.
@return Filtered vector of VCFLines (duplicates removed).
*/
std::vector<VCFLine> removeDuplicates(std::vector<VCFLine> variants);


/*
Filters homozygous variants (no variability across samples).

@param variants Lines read from VCF file.
@return Filtered vector of VCFLines (homozygous calls removed).
*/
std::vector<VCFLine> filterHomozygousVariants(std::vector<VCFLine> &variants);