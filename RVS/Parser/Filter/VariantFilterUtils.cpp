#include "VariantFilterUtils.h"

std::vector<VCFLine> missingTest(std::vector<VCFLine> &variants, VectorXd &G, int ngroup, double missingThreshold){

	int before = variants.size();
	std::vector<VCFLine> filteredVariants;

	for (int i = 0; i < variants.size(); i++) {

		//TODO: calculate percent missing by group?

		std::vector<double> counter(ngroup, 0);
		std::vector<double> n(ngroup, 0);

		for (int i = 0; i < variants[i].likelihood.size(); i++) {

			n[G[i]]++;
			if (variants[i].likelihood[i].missing)
				counter[G[i]]++;
		}


		bool valid = true;

		for (int i = 0; i < ngroup; i++) {
			if (counter[i] / n[i] > missingThreshold) {
				valid = false;
				break;
			}
		}

		if (valid)
			filteredVariants.push_back(variants[i]);
	}

	int nremoved = before - filteredVariants.size();
	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s +
			" filtered by due to missing information (threshold = " + std::to_string(missingThreshold) + ").");
	}

	return filteredVariants;
}


std::vector<VCFLine> filterTest(std::vector<VCFLine> &variants) {
	int before = variants.size();
	std::vector<VCFLine> filteredVariants;

	for (int i = 0; i < variants.size(); i++) {
		if (variants[i].filter == "PASS")
			filteredVariants.push_back(variants[i]);
	}

	int nremoved = before - filteredVariants.size();
	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s + " were filtered by FILTER column (do not PASS).");
	}

	return filteredVariants;
}

std::vector<VCFLine> snpTest(std::vector<VCFLine> &variants) {

	int before = variants.size();
	std::vector<VCFLine> filteredVariants;

	for (int i = 0; i < variants.size(); i++) {
		if (validBase(variants[i].alt) && validBase(variants[i].ref))
			filteredVariants.push_back(variants[i]);
	}

	int nremoved = before - filteredVariants.size();
	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s + " were filtered by ALT and REF column (SNPs were retained).");

	}

	return filteredVariants;

}

std::vector<VCFLine> filterMinorAlleleFrequency(std::vector<VCFLine> &variants, double mafCutoff, bool keepCommon) {
	double maf;
	int before = variants.size();
	std::vector<VCFLine> filteredVariants;

	for (size_t i = 0; i < variants.size(); i++) {
		maf = 0.5 * variants[i].P[1] + variants[i].P[2];
		if (maf > 0.5)
			maf = 1 - maf;

		if ((keepCommon && maf > mafCutoff) || (!keepCommon && (maf < mafCutoff)))
			filteredVariants.push_back(variants[i]);
	}

	int nremoved = before - filteredVariants.size();

	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";

		if (keepCommon) {
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

std::vector<VCFLine> removeDuplicates(std::vector<VCFLine> variants) {

	std::vector<VCFLine> filteredVariants;
	int before = variants.size();

	VCFLine lastVariant;

	if (before > 0) {
		lastVariant = variants[0];
		filteredVariants.push_back(lastVariant);
	}

	for (int i = 0; i < variants.size(); i++) {
		if (variants[i].loc == lastVariant.loc && variants[i].chr == lastVariant.chr)
			continue;

		lastVariant = variants[i];
		filteredVariants.push_back(lastVariant);
	}

	int nremoved = before - filteredVariants.size();

	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s +
			" filtered by due to duplication (multiple variants at same genomic position).");
	}

	return filteredVariants;
}

std::vector<VCFLine> filterHomozygousVariants(std::vector<VCFLine> &variants) {

	int i, j;

	double mean, n, sd;
	std::vector<VCFLine> filteredVariants;
		int before = variants.size();


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

