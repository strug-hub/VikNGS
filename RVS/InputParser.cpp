#include "stdafx.h"
#include "RVS.h"
#include "InputParser.h"


#include <set>


std::vector<SNP> parseInput(std::string vcfDir, std::string infoDir, double mafCut, 
	VectorXd &Y, VectorXd &G, VectorXd &H, MatrixXd &Z) {

	std::map<std::string, int> IDmap = getSampleIDMap(vcfDir);
	parseInfo(infoDir, IDmap, Y, G, H, Z);
	
	std::vector<VCFLine> variants = parseVCFLines(vcfDir);

	variants = filterVariants(variants, G, 0.2);
	std::sort(variants.begin(), variants.end(), lineCompare);
	variants = removeDuplicates(variants);

		
	for (int i = 0; i < variants.size(); i++) {

		//variants[i].print();

	}


	std::vector<SNP> s;
	return s;
}


/*
Parses a VCF file to get the expected probabilities of the genotypes E(G_ij | D_ij)
and expected minor allele frequency.
@param vcfDir Path of VCF file.
@param mafCut Cut-off value for minor allele frequency
@param sample Vector with sample information.
@return Vector of SNPs parsed from VCF file.

std::vector<SNP> processVCF(std::string vcfDir, std::string sampleInfoDir, double mafCut, std::vector<Sample> sample) {

	std::vector<SNP> snps = parseVCF(vcfDir, 9, 0.2, sample);

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
*/






