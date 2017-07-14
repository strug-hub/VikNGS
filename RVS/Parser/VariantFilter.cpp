#include "stdafx.h"
#include "InputParser.h"

inline bool validBase(std::string base) {
	return base == "T" ||  base == "A" || base == "C" || base == "G";
}

//TODO: calculate percent missing by group?
inline double missingTest(VCFLine variant, VectorXd &G, int ngroup, double missingThreshold) {
	
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

std::vector<VCFLine> filterVariants(std::vector<VCFLine> variants, VectorXd &G, double missingThreshold) {
	int filterCounter = 0;
	int missingCounter = 0;
	int indelCounter = 0;

	int ngroup = G.maxCoeff() + 1;

	std::vector<VCFLine> filteredVariants;
	VCFLine variant;

	for (int i = 0; i < variants.size(); i++) {
		variant = variants[i];

		if (variant.filter != "PASS") {
			filterCounter++;
			continue;
		}

		if (!validBase(variant.alt) || !validBase(variant.ref)) {
			indelCounter++;
			continue;
		}

		if (!missingTest(variant, G, ngroup, missingThreshold)) {
			missingCounter++;
			continue;
		}

		filteredVariants.push_back(variant);
	}

	std::cout << filterCounter;
	std::cout << " Variant(s) were filtered by FILTER column.\n";
	std::cout << indelCounter;
	std::cout << " Variant(s) were filtered by ALT and REF column.\n";
	std::cout << missingCounter;
	std::cout << " Variant(s) were filtered by missing threshold.\n";

	return filteredVariants;
}

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

	std::cout << variants.size() - filteredVariants.size();
	std::cout << " Variant(s) were filtered by duplication.\n";

	return filteredVariants;
}

std::vector<VCFLine> filterHomozygousVariants(std::vector<VCFLine> &variants) {

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
	
	std::cout << variants.size() - filteredVariants.size();
	std::cout << " Variant(s) were removed because of homozygous call in all samples.\n";

	return filteredVariants;
}

/*
Calculates the minor allele frequency (MAF) from conditional expected genotype probability.

@param snps A vector of SNPs.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param common Indicates common or rare variants, (common = true, rare = false).
@return None.
@effect Sets maf for each SNP in snps.
*/
std::vector<VCFLine> filterMinorAlleleFrequency(std::vector<VCFLine> &variants, double mafCutoff, bool common) {
	double maf;
	std::vector<VCFLine> filteredVariants;

	for (size_t i = 0; i < variants.size(); i++) {
		maf = 0.5 * variants[i].P[1] + variants[i].P[2];
		if (maf > 0.5)
			maf = 1 - maf;

		if (common && maf >= mafCutoff || !common && (maf < mafCutoff))
			filteredVariants.push_back(variants[i]);
	}
		
	std::cout << variants.size() - filteredVariants.size();
	std::cout << " Variant(s) filtered by minor allele frequency.\n";

	return filteredVariants;
}
