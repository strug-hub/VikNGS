#include "stdafx.h"
#include "RVS.h"
#include "MemoryMapped/MemoryMapped.h"

#include <iostream>  
#include <string>
#include <vector>

#include "TestSet.h"

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
			g.groupIndex = sample[i].groupIndex;

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
#include <iomanip>

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
	generateForR(sample, snps);


	std::vector<Group> group = calcGroups(sample, snps);

	collapseVariants(snps, collapse);

	generateForR(sample, snps);

	std::vector<bool> IDmap;
	for (size_t i = 0; i < snps.size(); i++) {
		IDmap.push_back(sample[i].y);
	}

	std::vector<double> pvals = runCommonTest(snps, sample, group);
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

	std::vector<SNP> toTest;
	toTest.push_back(snps[0]);
	toTest.push_back(snps[1]);
	toTest.push_back(snps[2]);
	toTest.push_back(snps[3]);
	toTest.push_back(snps[4]);


	std::cout << runRareTest(toTest, sample, group, 100, true);
	
//	auto t = startTime();
//	pvals = RVSbtrap(snps, sample, 1000000, true, true);
//	endTime(t, "btrp=1000000");

	

	std::cout << "done...>";

	//keep console open while debugging
	//TODO: be sure to remove eventually!
	while (true) {}
	return 0;
}