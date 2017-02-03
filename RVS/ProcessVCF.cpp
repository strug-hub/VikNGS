#include "stdafx.h"
#include "RVS.h"
#include "MemoryMapped/MemoryMapped.h"

#include <iostream>  
#include <string>
#include <vector>

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
			if (EG[j] != NULL) {
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
		if (snps[i].maf == NULL) 
			snps.erase(snps.begin() + i);
		else 
			i++;
	}
}

/*
Gets the expected probabilities of the genotypes E(G_ij | D_ij) from a VCF file
*/
std::vector<SNP> vcf_process(std::string vcfDir, std::string caseIDDir,
	double mafCut, bool common, std::vector<bool> IDmap) {
	
	std::vector<SNP> snps = parseAndFilter(vcfDir, 9, 0.2, IDmap);

	if (snps.size() == 0)
		std::cout << "No SNPs left after filtering .vcf file\n";

	getExpGeno(snps);
	
	if (snps.size() == 0)
		std::cout << "No SNPs left removing homozygous variants\n";

	getExpMAF(snps, mafCut, true);

	if (snps.size() == 0)
		std::cout << "No SNPs left after applying MAF condition\n";

	std::cout << "==========\n";

	return snps;
}

int main() {
	//TODO: take as input from command line
	//---------------------------------------
	bool common = true;
	double mafCut = 0.05;
	std::string vcfDir = "C:/Users/Scott/Desktop/RVS-master/example/example_1000snps.vcf";
	std::string caseIDDir = "C:/Users/Scott/Desktop/RVS-master/example/caseID.txt";
	//---------------------------------------

	std::vector<bool> IDmap = getIDs(vcfDir, caseIDDir, 9);

	std::vector<SNP> snps = vcf_process(vcfDir, caseIDDir, mafCut, common, IDmap);

	calcMeanVar(IDmap, snps);
	std::vector<double> pvals = RVSasy(snps, IDmap, true);
	//RVSbtrap(snps, IDmap, true, 1000);

	std::cout << "Chr\tLoc\tMAF\tp-value\n";
	for (size_t i = 0; i < snps.size(); i++) {
		if (snps[i].maf != NULL) {
			std::cout << snps[i].chr;
			std::cout << '\t';
			std::cout << snps[i].loc;
			std::cout << '\t';
			std::cout << snps[i].maf;
			std::cout << '\t';
			std::cout << pvals[i];
			std::cout << '\n';
		}
	}

	RVSrare(snps, IDmap, 10000);

	//keep console open while debugging
	//TODO: be sure to remove eventually!
	while (true) {}
	return 0;
}