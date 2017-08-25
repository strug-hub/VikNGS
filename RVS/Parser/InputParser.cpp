#pragma once
#include "InputParser.h"
#include <fstream>

/*
Uses EM algorithm to estimate the genotype frequencies in the sample

@param likelihood A vector with genotype likelihoods.
@return A vector with three doubles to stand for probability of 0, 1 or 2 minor alleles.
*/
VectorXd calcEM(std::vector<GenotypeLikelihood> &likelihood) {
	double p = 0.15;
	double q = 0.15;
	double qn = 1;
	double pn = 0;
	double dn = 0;
	double d = 0;

	int glCounter = 0;

	double Ep;
	double Eq;
	double pD;

	int k = 0;

	while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

		d = 1 - p - q;
		glCounter = 0;
		Ep = 0;
		Eq = 0;

		for (int i = 0; i < likelihood.size(); i++) {
			if (likelihood[i].missing)
				continue;

			glCounter++;

			pD = 1 / (p * likelihood[i].L00 + q * likelihood[i].L01 + d * likelihood[i].L11);
			Ep += p * likelihood[i].L00 * pD;
			Eq += q * likelihood[i].L01 * pD;
		}

		pn = p;
		qn = q;
		dn = 1 - q - p;
		p = Ep / glCounter;
		q = Eq / glCounter;

		k++;
		if (k == 1000)
			break;
	}

	VectorXd freq(3);
	freq[0] = p;
	freq[1] = q;
	freq[2] = 1 - p - q;

	return freq;
}

/*
Calculates the conditional expected genotype probability E( P(G_ij | D_ij) )
given the genotype likelihoods P(D_ij | G_ij = g) and frequencies.

E( P(G_ij | D_ij) ) = sum from g=0 to 2 P(G_ij = g | D_ij),
where P(G_ij = g | D_ij) = P(D_ij | G_ij = g) * P(G_ij = g)/P(D_ij).

@param snp SNP with genotype likelihoods P(D_ij | G=AA, Aa or aa)} for one locus.
@return A vector containing conditional expectation probability E( P(G_ij | D_ij) ).
*/
VectorXd calcEG(std::vector<GenotypeLikelihood> &likelihood, VectorXd &p) {

	double m0;
	double m1;
	double m2;
	double m;

	VectorXd EG(likelihood.size());

	for (int i = 0; i < likelihood.size(); i++) {

		if (!likelihood[i].missing) {
			m0 = likelihood[i].L00 * p[0];
			m1 = likelihood[i].L01 * p[1];
			m2 = likelihood[i].L11 * p[2];
			m = 1 / (m0 + m1 + m2);

			EG[i] = m1*m + 2 * m2*m;
		}
		else
			EG[i] = NAN;
	}

	return EG;
}

/*
Generates the expected probabilities of the genotypes E(G_ij | D_ij).
Using the genotype likelihood for case and controls from VCF file to generates the population frequency by calling function calcEM
and then use it to calculate the expected genotype probabilities E(G_ij | D_ij) by calling function calcEG.
Variants with homozygous call in the whole sample (standard deviation of their E(G_ij | D_ij) < 10^4) are removed.

@param snps A vector of SNPs.
@return None.
@effect Sets p for each SNP in snps and removes SNPs from snps with homozygous calls.
*/
std::vector<VCFLine> calculatedExpectedGenotypes(std::vector<VCFLine> &variants) {

	for (size_t i = 0; i < variants.size(); i++) {
		variants[i].P = calcEM(variants[i].likelihood);
		variants[i].expectedGenotype = calcEG(variants[i].likelihood, variants[i].P);
	}

	return variants;
}

inline bool fileExists(const std::string& name) {

	std::ifstream f(name.c_str());
	if (f.good())
		return true;
	
	std::string message = "File does not exist or cannot be opened: ";
	message.append(name);
	printError(message);
	return false;
}

bool parseInput(std::string vcfDir, std::string infoDir, std::string bedDir, double mafCutoff, bool common,
	MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd &P,
	std::vector<std::vector<int>> & interval) {

	if (!fileExists(vcfDir))
		return false;
	if (!fileExists(infoDir))
		return false;

	std::map<std::string, int> IDmap = getSampleIDMap(vcfDir);

	if (IDmap.size() <= 0) {
		printError("No samples could be found in the VCF file header! Exiting.");
		return false;
	}

	bool valid = parseInfo(infoDir, IDmap, Y, Z, G, readGroup);
	if (!valid) {
		printError("Failed to parse sample data file. Exiting.");
		return false;
	}
	std::vector<VCFLine> variants;
	try {
		variants = parseVCFLines(vcfDir);
	}
	catch (...) {
		printError("Failed to parse VCF. Exiting.");
		return false;
	}

	variants = filterVariants(variants, G, 0.2);
	std::sort(variants.begin(), variants.end(), lineCompare);
	variants = removeDuplicates(variants);
	variants = calculatedExpectedGenotypes(variants);
	variants = filterHomozygousVariants(variants);
	variants = filterMinorAlleleFrequency(variants, mafCutoff, common);
	std::cout << variants.size();
	std::cout << " SNPs remain after filtering\n";

	if (variants.size() <= 0)
		return false;

	MatrixXd x(variants[0].likelihood.size(), variants.size());
	for (int i = 0; i < variants.size(); i++) 
		x.col(i) = variants[i].expectedGenotype;
	
	MatrixXd p(variants.size(), 3);
	for (int i = 0; i < variants.size(); i++)
		p.row(i) = variants[i].P;
	
	X = x;
	P = p;
	
	interval = parseIntervals(bedDir, variants);

	for (int i = 0; i < interval.size(); i++) {
		for (int j = 0; j < interval[i].size(); j++) {

			std::cout << interval[i][j];
			std::cout << "   ";

		}

		std::cout << "\n";

	}

	return true;
}