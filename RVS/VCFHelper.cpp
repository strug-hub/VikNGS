#include "stdafx.h"
#include "RVS.h"
#include "MemoryMapped/MemoryMapped.h"
#include "InputParser.h"

#include <iostream>  
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <map>






/*
Separates the case and control IDs

@param vcfDir Full directory path of the VCF file.
@param sampleInfoDir Full directory path of the sampleInfo file.
@param ncolID The number of columns before the sample IDs start in the last line of the headers in the VCF file.
@return A vector of Samples with information parsed from sampleInfo file

std::vector<Sample> getSampleInfo(std::string vcfDir, std::string sampleInfoDir, int ncolID) {

	//open VCF file and find the line with column names
	MemoryMapped vcf(vcfDir);
	int pos = 0;
	std::vector<std::string> IDs = parseHeader(vcf, pos);
	vcf.close();


	//open sampleInfo file
	MemoryMapped sampleInfo(sampleInfoDir);
	pos = 0;

	int nIDs = IDs.size() - ncolID;
	std::vector<Sample> sample(nIDs);

	//extract every sample from sampleInfo file, assumes one sample per line
	std::string ID;
	std::vector<std::string> lineSplit;
	int groupIndex = 0;

	while (true) {

		lineSplit = readLine(sampleInfo, pos);

		if (lineSplit.size() < 2)
			break;

		Sample s;
		s.ID = lineSplit[0];
		s.y = std::stod(lineSplit[1]);
		s.groupID = lineSplit[2];
		s.groupIndex = groupIndex;
		groupIndex++;

		std::vector<std::string> cov;
		for (size_t i = 4; i < lineSplit.size(); i++)
			cov.push_back(lineSplit[i]);

		s.covariates = handleCovariates(cov);

		std::string depth = lineSplit[3];

		if (depth.find("H") != std::string::npos ||
			depth.find("h") != std::string::npos)
			s.hrg = 1;

		int index = findIndex(s.ID, IDs) - ncolID;

		sample[index] = s;
	}

	std::vector<std::string> uniqueGroups;
	std::map<std::string, int> groupIDMap;
	int index = 0;

	for (size_t i = 0; i < sample.size(); i++) {
	
		std::string g = sample[i].groupID;
		if (std::find(uniqueGroups.begin(), uniqueGroups.end(), g) == uniqueGroups.end()){
			uniqueGroups.push_back(sample[i].groupID);
			groupIDMap[g] = index;
			index++;
		}

		sample[i].groupIndex = groupIDMap[g];

	}

	//TODO: dummyproofing parsing here before finalized product


	
	//if (lineCounter != countCase)
	//std::cout << "Warning: " + std::to_string(lineCounter) + " lines were counted in case ID file but only " +
	//std::to_string(countCase) + " were found to correspond to columns in .vcf file.\n";


	sampleInfo.close();
	return sample;
}
*/
	

/*
Generates the expected probabilities of the genotypes E(G_ij | D_ij).
Using the genotype likelihood for case and controls from VCF file to generates the population frequency by calling function calcEM
and then use it to calculate the expected genotype probabilities E(G_ij | D_ij) by calling function calcEG.
Variants with homozygous call in the whole sample (standard deviation of their E(G_ij | D_ij) < 10^4) are removed.

@param snps A vector of SNPs.
@return None.
@effect Sets p for each SNP in snps and removes SNPs from snps with homozygous calls.
*/
void getExpGeno(std::vector<SNP> &snps) {

	std::vector<double> EG;
	std::vector<size_t> filterIndex;

	double mean;
	double sdSum;

	double n;
	int filterCount = 0;

	for (size_t i = 0; i < snps.size(); i++) {
		std::vector<double> p;

		//calculates E(G_ij | D_ij)
		p = calcEM(snps[i]);
		snps[i].p = p;
		EG = calcEG(snps[i]);
		snps[i].EG = EG;

		//check and filter if variant is homozygous
		mean = 0;
		n = 0;
		for (size_t j = 0; j < EG.size(); j++) {
			if (!isnan(EG[j])) {
				mean += EG[j];
				n++;
			}
		}
		mean = mean / n;

		sdSum = 0;
		for (size_t j = 0; j < EG.size(); j++)
			if (EG[j] != NULL)
				sdSum += pow((EG[j] - mean), 2);

		if (1e-8*(n - 1) > sdSum)
			filterIndex.push_back(i);
	}

	for (size_t i = 0; i < filterIndex.size(); i++) {
		snps.erase(snps.begin() + filterIndex[i]);
		filterCount++;
	}

	if (filterCount > 0) {
		std::cout << filterCount;
		std::cout << " SNPs were removed because of homozygous call in all samples\n";
	}
	return;
}

/*
Calculates the minor allele frequency (MAF) from conditional expected genotype probability.

@param snps A vector of SNPs.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param common Indicates common or rare variants, (common = true, rare = false).
@return None.
@effect Sets maf for each SNP in snps.
*/
void getExpMAF(std::vector<SNP> &snps, double mafCut, bool common) {
	double maf;
	int mafCounter = 0;


	for (size_t i = 0; i < snps.size(); i++) {
		maf = 0.5 * snps[i].p[1] + snps[i].p[2];
		if (maf > 0.5)
			maf = 1 - maf;

		if (common) {
			if (maf >= mafCut) {
				mafCounter++;
				snps[i].maf = maf;
			}
		}
		else {
			if (maf < mafCut) {
				mafCounter++;
				snps[i].maf = maf;
			}
		}
	}

	std::cout << mafCounter;
	std::cout << " out of ";
	std::cout << snps.size();
	std::cout << " variants satisfy the MAF condition provided\n";

	//remove SNPs that failed MAF condition
	int i = 0;
	while (i < snps.size()) {
		if (isnan(snps[i].maf))
			snps.erase(snps.begin() + i);
		else
			i++;
	}
}







/*
Uses EM algorithm to estimate the genotype frequencies in the sample

@param snp A SNP with genotype likelihoods.
@return A vector with three doubles to stands for probability of 0, 1 or 2 minor alleles.
*/
std::vector<double> calcEM(SNP &snp) {
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

		for (int i = 0; i < snp.gl.size(); i++) {
			if (isnan(snp.L00(i)))
				continue;

			glCounter++;

			pD = 1 / (p * snp.L00(i) + q * snp.L01(i) + d * snp.L11(i));
			Ep += p * snp.L00(i) * pD;
			Eq += q * snp.L01(i) * pD;
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

	std::vector<double> freq;
	freq.push_back(p);
	freq.push_back(q);
	freq.push_back(1 - p - q);

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
std::vector<double> calcEG(SNP &snp) {

	double m0;
	double m1;
	double m2;
	double m;

	std::vector<double> EG;

	for (int i = 0; i < snp.gl.size(); i++) {

		if (!isnan(snp.L00(i))) {
			m0 = snp.L00(i) * snp.p[0];
			m1 = snp.L01(i) * snp.p[1];
			m2 = snp.L11(i) * snp.p[2];
			m = 1 / (m0 + m1 + m2);

			EG.push_back(m1*m + 2 * m2*m);
		}
		else
			EG.push_back(NAN);

	}

	return EG;
}