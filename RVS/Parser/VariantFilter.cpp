#include "InputParser.h"

/*
Verifies whether or not a string is a single nucleotide base (A, T, C or G).

@param base String to verify.
@return True if base is valid.
*/
inline bool validBase(std::string base) {
	return base == "T" ||  base == "A" || base == "C" || base == "G";
}

/*
Counts the number of missing sample data for a variant and returns true if any group is missing more than missingThreshold.

@param variant Line read from VCF file.
@param G Vector of group ID.
@param ngroup Number of unique groups in G.
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@return True if base is missing data above the threshold.
*/
inline double missingTest(VCFLine variant, VectorXd &G, int ngroup, double missingThreshold) {
	
	//TODO: calculate percent missing by group?

	std::vector<double> counter(ngroup, 0);
	std::vector<double> n(ngroup, 0);

	for (int i = 0; i < variant.likelihood.size(); i++) {

		n[G[i]]++;
		if (variant.likelihood[i].missing)
			counter[G[i]]++;
	}

	for (int i = 0; i < ngroup; i++)
		if (counter[i] / n[i] > missingThreshold)
			return false;

	return true;
}

/*
Filters variants read from VCF file based on missing sample data, failing to have 'PASS' in FILTER column
and whether the variant is an insertion or deletion (indel).

@param req Request object containing filtering parameters
@param variants Lines read from VCF file.
@param G Vector of group ID.
@return Filtered vector of VCFLines.
*/
std::vector<VCFLine> filterVariants(Request req, std::vector<VCFLine> variants, VectorXd &G) {
	int failCounter = 0;
	int missingCounter = 0;
	int indelCounter = 0;

	//assumes number of groups = highest group ID
	int ngroup = G.maxCoeff() + 1;

	std::vector<VCFLine> filteredVariants;
	VCFLine variant;

	for (int i = 0; i < variants.size(); i++) {
		variant = variants[i];

		if (req.mustPASS && variant.filter != "PASS") {
			failCounter++;
			continue;
		}

		if (req.onlySNPs && !validBase(variant.alt) || !validBase(variant.ref)) {
			indelCounter++;
			continue;
		}

		if (!missingTest(variant, G, ngroup, req.missingThreshold)) {
			missingCounter++;
			continue;
		}

		filteredVariants.push_back(variant);
	}

	if (req.mustPASS && failCounter > 0) {
		std::string s = failCounter > 1 ? "s" : "";
		printInfo(std::to_string(failCounter) + " variant" + s + 
			" were filtered by FILTER column (do not PASS).");
	}
	if (req.onlySNPs && indelCounter > 0) {
		std::string s = indelCounter > 1 ? "s" : "";
		printInfo(std::to_string(indelCounter) + " variant" + s + 
			" were filtered by ALT and REF column (indel variant" + s + " removed).");
	}
	if (missingCounter > 0) {
		std::string s = missingCounter > 1 ? "s" : "";
		printInfo(std::to_string(missingCounter) + " variant" + s + 
			" filtered by due to missing information (threshold = " + std::to_string(req.missingThreshold) + ").");
	}

	return filteredVariants;
}

/*
Filters duplicated variants (same position and chromosome).

@param variants Lines read from VCF file.
@requires variants to be sorted in order of chr, loc.
@return Filtered vector of VCFLines (duplicates removed).
*/
std::vector<VCFLine> removeDuplicates(std::vector<VCFLine> variants) {

	std::vector<VCFLine> filteredVariants;
	VCFLine lastVariant;

	if (variants.size() > 0) {
		lastVariant = variants[0];
		filteredVariants.push_back(lastVariant);
	}

	for (int i = 0; i < variants.size(); i++) {
		if (variants[i].loc == lastVariant.loc && variants[i].chr == lastVariant.chr)
			continue;

		lastVariant = variants[i];
		filteredVariants.push_back(lastVariant);
	}

	int nremoved =  variants.size() - filteredVariants.size();

	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s +
			" filtered by due to duplication (multiple variants at same genomic position).");
	}

	return filteredVariants;
}

/*
Filters homozygous variants (no variability across samples).

@param variants Lines read from VCF file.
@return Filtered vector of VCFLines (homozygous calls removed).
*/std::vector<VCFLine> filterHomozygousVariants(std::vector<VCFLine> &variants) {

	int i, j;

	double mean, n, sd;
	std::vector<VCFLine> filteredVariants;

	for (i = 0; i < variants.size(); i++) {
		VectorXd EG = variants[i].expectedGenotype;

		//check and filter if variant is homozygous
		mean = 0;
		n = 0;
		for (j = 0; j < EG.size(); j++) {
			if (!isnan(EG[j])) {
				mean += EG[j];
				n++;
			}
		}
		mean = mean / n;

		sd = 0;
		for (j = 0; j < EG.size(); j++)
			if (!isnan(EG[j]))
				sd += pow((EG[j] - mean), 2);

		if (1e-8*(n - 1) < sd)
			filteredVariants.push_back(variants[i]);
	}
	
	int nremoved = variants.size() - filteredVariants.size();
	
	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s +
			" filtered due to homozygous call in all samples.");
	}

	return filteredVariants;
}

/*
Filters variants based on minor allele frequency. Removes rare or common variants depending on the value of common.

@param variants Lines read from VCF file.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param common Indicates whether to keep common or rare variants, (common = true, rare = false).
@return Filtered vector of VCFLines (either common or rare variants removed).
*/
std::vector<VCFLine> filterMinorAlleleFrequency(std::vector<VCFLine> &variants, double mafCutoff, bool common) {
	double maf;
	std::vector<VCFLine> filteredVariants;

	for (size_t i = 0; i < variants.size(); i++) {
		maf = 0.5 * variants[i].P[1] + variants[i].P[2];
		if (maf > 0.5)
			maf = 1 - maf;

		if (common && maf > mafCutoff || !common && (maf < mafCutoff))
			filteredVariants.push_back(variants[i]);
	}
		

	int nremoved = variants.size() - filteredVariants.size();

	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";

		if (common) {
			printInfo(std::to_string(nremoved) + " rare variant" + s +
				" filtered (minor allele frequency < " + std::to_string(mafCutoff) + ").");
		}
		else {
			printInfo(std::to_string(nremoved) + " common variant" + s +
				" filtered (minor allele frequency > " + std::to_string(mafCutoff) + ").");
		}
	}

	return filteredVariants;
}
