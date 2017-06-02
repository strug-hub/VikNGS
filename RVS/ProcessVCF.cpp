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
Parses a VCF file to get the expected probabilities of the genotypes E(G_ij | D_ij) 
and expected minor allele frequency.
@param vcfDir Path of VCF file.
@param mafCut Cut-off value for minor allele frequency
@param sample Vector with sample information.
@return Vector of SNPs parsed from VCF file.
*/
std::vector<SNP> processVCF(std::string vcfDir, std::string sampleInfoDir, double mafCut, std::vector<Sample> sample) {
	
	std::vector<SNP> snps = parseAndFilter(vcfDir, 9, 0.2, sample);

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

/*
Creates a vector of Groups based on the sample information

@param sample Vector with sample information.
@param snps Vector of SNPs. TODO: is this actually needed?
@return Vector of Groups.
*/
std::vector<Group> calcGroups(std::vector<Sample> &sample, std::vector<SNP> &snps) {
	std::vector<Group> group;
	std::vector<bool> done;
	for (size_t i = 0; i < sample.size(); i++)
		done.push_back(false);

	for (size_t i = 0; i < sample.size(); i++){
		if (!done[i]) {
			Group g;
			g.ID = sample[i].groupID;
			g.hrg = sample[i].hrg;
			
			for (size_t j = i; j < sample.size(); j++) 
				if (sample[j].groupID == g.ID && sample[j].hrg == g.hrg) {
					g.index.push_back(j);
					done[j] = true;
				}
		
			group.push_back(g);
		}
	}

	return group;
}


#include <iostream>
#include <fstream>

void generateForR(std::vector<Sample> sample, std::vector<SNP> snps) {
	std::ofstream X("C:/Users/Scott/Desktop/RVS-master/example/X.txt");
	std::ofstream Y("C:/Users/Scott/Desktop/RVS-master/example/Y.txt");
	std::ofstream P("C:/Users/Scott/Desktop/RVS-master/example/P.txt");
	std::ofstream M("C:/Users/Scott/Desktop/RVS-master/example/M.txt");
	std::ofstream Z("C:/Users/Scott/Desktop/RVS-master/example/Z.txt");

	if (M.is_open())
	{
		M << "ID\thrg\n";

		for (size_t i = 0; i < sample.size(); i++) {

			M << sample[i].groupID;
			M << '\t';
			M << sample[i].hrg;
			M << '\n';
		}
		M.close();
	}


	if (Y.is_open())
	{
		for (size_t i = 0; i < sample.size(); i++) {
			Y << sample[i].y;
			Y << '\n';
		}
		Y.close();
	}

	if (X.is_open())
	{
		for (size_t i = 0; i < snps.size(); i++) {
			for (size_t j = 0; j < snps[i].EG.size(); j++) {

				if (isnan(snps[i].EG[j]))
					X << "NA";
				else
					X << snps[i].EG[j];

				if(j <snps[i].EG.size()-1)
					X << '\t';
			}
			X << '\n';
		}

		X.close();
	}
	if (P.is_open())
	{
		for (size_t i = 0; i < snps.size(); i++) {
			P << snps[i].p[0];
			P << '\t';
			P << snps[i].p[1];
			P << '\t';
			P << snps[i].p[2];
			P << '\n';
		}
		
		P.close();
	}

	if (Z.is_open())
	{
		for (size_t i = 0; i < sample.size(); i++) {

			for (size_t j = 0; j < sample[i].numeric_cov.size(); j++) {

				Z << sample[i].numeric_cov[j];
				if (sample[i].factor_cov.size() == 0 &&
					sample[i].numeric_cov.size() - 1 != j)
					Z << '\t';
			}

			for (size_t j = 0; j < sample[i].factor_cov.size(); j++) {

				Z << sample[i].factor_cov[j];

				if (sample[i].factor_cov.size() - 1 != j)
					Z << '\t';
			}

			Z << '\n';
		}

		Z.close();
	}
}


int main() {
	//TODO: take as input from command line
	//---------------------------------------
	double mafCut = 0.05;
	std::string vcfDir = "C:/Users/Scott/Desktop/RVS-master/example/example_1000snps.vcf";
	std::string sampleInfoDir = "C:/Users/Scott/Desktop/RVS-master/example/sampleInfo.txt";
	std::string bedDir = "C:/Users/Scott/Desktop/RVS-master/example/chr11.bed";
	//---------------------------------------


	//TODO: check to see if file can be opened when another application is using it (excel)
	//TODO: test windows vs unix EOF characters, doesn't seem to work well with windows

	std::vector<Interval> collapse = getIntervals(bedDir);
	std::vector<Sample> sample = getSampleInfo(vcfDir, sampleInfoDir, 9);
	std::vector<SNP> snps = processVCF(vcfDir, sampleInfoDir, mafCut, sample);
	std::vector<Group> group = calcGroups(sample, snps);

	collapseVariants(snps, collapse);

	CovariateRegression(snps[0], sample);

	generateForR( sample, snps);

	std::vector<bool> IDmap;
	for (size_t i = 0; i < snps.size(); i++) {
		IDmap.push_back(sample[i].y);
	}


	std::vector<double> pvals = RVSasy(snps, sample, group);
	//std::vector<double> pvals = RVSbtrap(snps, sample, group, 10000, true, true);

	/*
	for (size_t i = 2; i <= 6; i++) {
	int nboot = pow(10, i);
	auto t = startTime();
	RVSbtrap(snps, sample, nboot, true);
	endTime(t, std::to_string(nboot));
	}
	*/

	std::cout << "Case\tChr\tLoc\tMAF\tp-value\n";
	for (size_t i = 0; i < snps.size(); i++) {
		if (!isnan(snps[i].maf)) {
			std::cout << sample[i].y;
			std::cout << '\t';
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

//	auto t = startTime();
//	pvals = RVSbtrap(snps, sample, 1000000, true, true);
//	endTime(t, "btrp=1000000");


	
	for (size_t i = 0; i < snps.size(); i += 5) {
		std::vector<SNP> snps2eval;

		for (size_t j = 0; j <= 4; j++) {
			snps2eval.push_back(snps[i+j]);

		}
		
		std::vector<double> p = RVSrare(snps2eval, sample, group, 100000);

		std::cout << p[0];
		std::cout << "\t";
		std::cout << p[1];
		std::cout << "\n";


	}

	

	std::cout << "done...>";

	//keep console open while debugging
	//TODO: be sure to remove eventually!
	while (true) {}
	return 0;
}