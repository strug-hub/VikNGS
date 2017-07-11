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
	std::cout << " SNPs were filtered by FILTER column.\n";
	std::cout << indelCounter;
	std::cout << " SNPs were filtered by ALT and REF column.\n";
	std::cout << missingCounter;
	std::cout << " SNPs were filtered by missing threshold.\n";

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
	std::cout << " SNPs were filtered by duplication\n";
	std::cout << filteredVariants.size();
	std::cout << " SNPs remain after filtering\n";

	return filteredVariants;
}