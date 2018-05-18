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
@param Y Vector of response variable.
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@return Filtered vector of Variants, variants missing above threshold removed.
*/
std::vector<Variant> missingTest(std::vector<Variant> &variants, VectorXd &Y, double missingThreshold, std::string family, std::string outputDir);


/*
Removes variants that do not contain "PASS" in the FILTER column.

@param variants Lines read from VCF file.
@return Filtered vector of Variants.
*/
std::vector<Variant> filterTest(std::vector<Variant> &variants, std::string outputDir);


/*
Filters out non-SNP variants.

@param variants Lines read from VCF file.
@return Filtered vector of Variants, SNPs retained.
*/
std::vector<Variant> snpTest(std::vector<Variant> &variants, std::string outputDir);


/*
Filters variants based on minor allele frequency. Removes rare or common variants depending on the value of keepCommon.

@param variants Lines read from VCF file.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param keepCommon Indicates whether to keep common or rare variants, (common = true, rare = false).
@return Filtered vector of Variants (either common or rare variants removed).
*/
std::vector<Variant> filterMinorAlleleFrequency(std::vector<Variant> &variants, double mafCutoff, bool keepCommon, std::string outputDir);

/*
Filters duplicated variants (same position and chromosome).

@param variants Lines read from VCF file.
@requires variants to be sorted in order of chr, loc.
@return Filtered vector of Variants (duplicates removed).
*/
std::vector<Variant> removeDuplicates(std::vector<Variant> variants, std::string outputDir);


/*
Filters homozygous variants (no variability across samples).

@param variants Lines read from VCF file.
@return Filtered vector of Variants (homozygous calls removed).
*/
std::vector<Variant> filterHomozygousVariants(std::vector<Variant> &variants, std::string outputDir);
