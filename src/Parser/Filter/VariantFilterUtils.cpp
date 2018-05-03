#include "VariantFilterUtils.h"
#include "../../Output/OutputHandler.h"

//case+control
//find percent missing for case and control separately
std::vector<Variant> missingTestCaseControl(std::vector<Variant> &variants, VectorXd &Y, double missingThreshold, std::string outputDir){

    std::vector<Variant> filteredVariants;
    std::vector<Variant> removedVariants;

    int nvariants = variants.size();
    int nsamp = Y.rows();

    double ncase = Y.sum();
    double ncontrol = nsamp - ncase;

    for (int i = 0; i < nvariants; i++) {

        double missingCase = 0;
        double missingControl = 0;

        for (int j = 0; j < nsamp; j++) {

            if (variants[i].likelihood[j].missing){
                if(Y[j] == 0)
                    missingControl++;
                else if(Y[j] == 1)
                    missingCase++;
            }
        }

        bool valid = true;

        if(missingCase / ncase > missingThreshold)
            valid = false;
        if(missingControl / ncontrol > missingThreshold)
            valid = false;

        if (valid)
            filteredVariants.push_back(variants[i]);
        else
            removedVariants.push_back(variants[i]);
    }

    int nremoved = removedVariants.size();
    if (nremoved > 0) {
        std::string s = nremoved > 1 ? "s" : "";
        printInfo(std::to_string(nremoved) + " variant" + s +
            " filtered by due to missing information (threshold = " + std::to_string(missingThreshold) + ").");
    }
    outputFiltered(removedVariants, "Filtered by missing data", outputDir);
    return filteredVariants;
}

//quantitative
//find percent missing for all samples
std::vector<Variant> missingTestQuantitative(std::vector<Variant> &variants, VectorXd &Y, double missingThreshold, std::string outputDir){

    std::vector<Variant> filteredVariants;
    std::vector<Variant> removedVariants;

    int nvariants = variants.size();
    int nsamp = Y.rows();

    for (int i = 0; i < nvariants; i++) {

        double nmissing= 0;

        for (int j = 0; j < nsamp; j++)
            if (variants[i].likelihood[j].missing)
                nmissing++;

        if (nmissing / (1.0*nsamp) > missingThreshold)
            removedVariants.push_back(variants[i]);
        else
            filteredVariants.push_back(variants[i]);
    }

    int nremoved = removedVariants.size();
    if (nremoved > 0) {
        std::string s = nremoved > 1 ? "s" : "";
        printInfo(std::to_string(nremoved) + " variant" + s +
            " filtered by due to missing information (threshold = " + std::to_string(missingThreshold) + ").");
    }
    outputFiltered(removedVariants, "Filtered by missing data", outputDir);
    return filteredVariants;
}

std::vector<Variant> missingTest(std::vector<Variant> &variants, VectorXd &Y, double missingThreshold, std::string family, std::string outputDir){

    if(family == "binomial")
        return missingTestCaseControl(variants, Y, missingThreshold, outputDir);
    else
       return missingTestQuantitative(variants, Y, missingThreshold, outputDir);
}

std::vector<Variant> filterTest(std::vector<Variant> &variants, std::string outputDir) {

    std::vector<Variant> filteredVariants;
    std::vector<Variant> removedVariants;

	for (int i = 0; i < variants.size(); i++) {
		if (variants[i].filter == "PASS")
			filteredVariants.push_back(variants[i]);
		else
            removedVariants.push_back(variants[i]);
	}

	int nremoved = removedVariants.size();
	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s + " were filtered by FILTER column (do not PASS).");
	}

	outputFiltered(removedVariants, "Filtered by PASS", outputDir);
	return filteredVariants;
}

std::vector<Variant> snpTest(std::vector<Variant> &variants, std::string outputDir) {

    std::vector<Variant> filteredVariants;
    std::vector<Variant> removedVariants;

	for (int i = 0; i < variants.size(); i++) {
		if (validBase(variants[i].alt) && validBase(variants[i].ref))
			filteredVariants.push_back(variants[i]);
		else
            removedVariants.push_back(variants[i]);
	}

	int nremoved = removedVariants.size();
	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s + " were filtered by ALT and REF column (SNPs were retained).");

	}
	outputFiltered(removedVariants, "Filtered for REF/ALT",outputDir);

	return filteredVariants;
}

std::vector<Variant> filterMinorAlleleFrequency(std::vector<Variant> &variants, double mafCutoff, bool keepCommon, std::string outputDir) {
	double maf;
    std::vector<Variant> filteredVariants;
    std::vector<Variant> removedVariants;

	for (size_t i = 0; i < variants.size(); i++) {
		maf = 0.5 * variants[i].P[1] + variants[i].P[2];

		if (maf > 0.5)
			maf = 1 - maf;

		if ((keepCommon && maf > mafCutoff) || (!keepCommon && (maf < mafCutoff)))
			filteredVariants.push_back(variants[i]);
		else
            removedVariants.push_back(variants[i]);
	}

	int nremoved = removedVariants.size();

	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";

		if (keepCommon) {
			printInfo(std::to_string(nremoved) + " rare variant" + s +
                " filtered (minor allele frequency less than " + std::to_string(mafCutoff) + ").");
		}
		else {
			printInfo(std::to_string(nremoved) + " common variant" + s +
                " filtered (minor allele frequency greater than " + std::to_string(mafCutoff) + ").");
		}
	}
	outputFiltered(removedVariants, "Filtered for minor allele frequency",outputDir);
	return filteredVariants;
}

std::vector<Variant> removeDuplicates(std::vector<Variant> variants, std::string outputDir) {

    std::vector<Variant> filteredVariants;
    std::vector<Variant> removedVariants;

    Variant lastVariant;

	if (filteredVariants.size() > 0) {
		lastVariant = variants[0];
		filteredVariants.push_back(lastVariant);
	}

	for (int i = 0; i < variants.size(); i++) {
		if (variants[i].pos == lastVariant.pos && variants[i].chr == lastVariant.chr){
            removedVariants.push_back(variants[i]);
			continue;
		}

		lastVariant = variants[i];
		filteredVariants.push_back(lastVariant);
	}

	int nremoved = removedVariants.size();

	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s +
			" filtered by due to duplication (multiple variants at same genomic position).");
	}

	outputFiltered(removedVariants, "Filtered for duplication",outputDir);
	return filteredVariants;
}

std::vector<Variant> filterHomozygousVariants(std::vector<Variant> &variants, std::string outputDir) {

	int i, j;

	double mean, n, sd;
    std::vector<Variant> filteredVariants;
    std::vector<Variant> removedVariants;

	for (i = 0; i < variants.size(); i++) {
		VectorXd EG = variants[i].expectedGenotype;

		//check and filter if variant is homozygous
		mean = 0;
		n = 0;
		for (j = 0; j < EG.size(); j++) {
			if (!std::isnan(EG[j])) {
				mean += EG[j];
				n++;
			}
		}
		mean = mean / n;

		sd = 0;
		for (j = 0; j < EG.size(); j++)
			if (!std::isnan(EG[j]))
				sd += pow((EG[j] - mean), 2);

		if (1e-8*(n - 1) < sd)
			filteredVariants.push_back(variants[i]);
		else
            removedVariants.push_back(variants[i]);
	}

	int nremoved = removedVariants.size();

	if (nremoved > 0) {
		std::string s = nremoved > 1 ? "s" : "";
		printInfo(std::to_string(nremoved) + " variant" + s +
			" filtered due to homozygous call in all samples.");
	}
	outputFiltered(removedVariants, "Filtered for lack of variability", outputDir);
	return filteredVariants;
}

