#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <string>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>

/*
void generateForR(std::vector<Sample> sample, std::vector<SNP> snps) {
	std::ofstream X("C:/Users/Scott/Desktop/RVS-master/example/X.txt");
	std::ofstream Y("C:/Users/Scott/Desktop/RVS-master/example/Y.txt");
	std::ofstream P("C:/Users/Scott/Desktop/RVS-master/example/P.txt");
	std::ofstream M("C:/Users/Scott/Desktop/RVS-master/example/M.txt");
	std::ofstream Z("C:/Users/Scott/Desktop/RVS-master/example/Z.txt");

	int precise = 12;

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
			Y << std::setprecision(precise) << sample[i].y;
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
					X << std::setprecision(precise) << snps[i].EG[j];

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
			P << std::setprecision(precise) << snps[i].p[0];
			P << '\t';
			P << std::setprecision(precise) << snps[i].p[1];
			P << '\t';
			P << std::setprecision(precise) << snps[i].p[2];
			P << '\n';
		}
		
		P.close();
	}

	if (Z.is_open())
	{
		for (size_t i = 0; i < sample.size(); i++) {

			for (size_t j = 0; j < sample[i].covariates.size(); j++) {
				Z << std::setprecision(precise) << sample[i].covariates[j];
				
				if(j+1 != sample[i].covariates.size())
					Z << '\t';
			}

			Z << '\n';
		}

		Z.close();
	}
}
*/

int main() {

	//simulate();

	//TODO: take as input from command line
	//---------------------------------------
	double mafCutoff = 0.05;
	std::string vcfDir = "C:/Users/Scott/Desktop/RVS-master/example/example_1000snps.vcf";
	std::string infoDir = "C:/Users/Scott/Desktop/RVS-master/example/sampleInfo.txt";
	std::string bedDir = "C:/Users/Scott/Desktop/RVS-master/example/chr11.bed";
	//---------------------------------------


	//TODO: check to see if file can be opened when another application is using it (excel)
	//TODO: test windows vs unix EOF characters, doesn't seem to work well with windows
	
	//TODO:
	//std::vector<Interval> collapse = getIntervals(bedDir);

	VectorXd Y, G; MatrixXd X, Z, P;
	std::map<int, int> readGroup;

	bool valid = parseInput(vcfDir, infoDir, mafCutoff, true, X, Y, Z, G, readGroup, P);
	if (!valid)
		return 0;

	std::vector<double> pvals = runCommonTest(X, Y, G, readGroup, P);
	//std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P);
	//std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P, 1000, true);

	std::cout << "Common Test p-values\n";
	for (size_t i = 0; i < pvals.size(); i++) {
			std::cout << pvals[i];
			std::cout << '\n';
	}

	std::vector<std::vector<double>> pval = runRareTest(X, Y, Z, G, readGroup, P, 50000, true);

	std::cout << "Rare Test p-values\n";
	std::cout << pval[0][0];
	std::cout << '\t';
	std::cout << pval[0][1];


	/*

	collapseVariants(snps, collapse);

	generateForR(sample, snps);
	auto t = startTime();
	pvals = RVSbtrap(snps, sample, 1000000, true, true);
	endTime(t, "btrp=1000000");


	*/

	//keep console open while debugging
	//TODO: be sure to remove eventually!

	std::cout << "done...>";

	while (true) {}
	return 0;
}