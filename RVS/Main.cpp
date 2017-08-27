#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <string>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>


void generateForR(MatrixXd X, VectorXd Y, MatrixXd Z, VectorXd G, MatrixXd P, std::map<int, int> readGroup) {
	std::ofstream Xfile("C:/Users/Scott/Desktop/RVS-master/example/X.txt");
	std::ofstream Yfile("C:/Users/Scott/Desktop/RVS-master/example/Y.txt");
	std::ofstream Pfile("C:/Users/Scott/Desktop/RVS-master/example/P.txt");
	std::ofstream Mfile("C:/Users/Scott/Desktop/RVS-master/example/M.txt");
	std::ofstream Zfile("C:/Users/Scott/Desktop/RVS-master/example/Z.txt");

	int precise = 12;

	if (Mfile.is_open())
	{
		Mfile << "ID\thrg\n";

		for (size_t i = 0; i < G.rows(); i++) {

			Mfile << G[i];
			Mfile << '\t';
			Mfile << readGroup[G[i]];
			Mfile << '\n';
		}
		Mfile.close();
	}


	if (Yfile.is_open())
	{
		for (size_t i = 0; i < Y.rows(); i++) {
			Yfile << std::setprecision(precise) << Y[i];
			Yfile << '\n';
		}
		Yfile.close();
	}

	if (Xfile.is_open())
	{
		for (size_t i = 0; i < X.cols(); i++) {
			for (size_t j = 0; j < X.rows(); j++) {

				if (isnan(X(j, i)))
					Xfile << "NA";
				else
					Xfile << std::setprecision(precise) << X(j, i);

				if (j < X.rows() - 1)
					Xfile << '\t';
			}
			Xfile << '\n';
		}

		Xfile.close();
	}

	if (Pfile.is_open())
	{
		for (size_t i = 0; i < P.rows(); i++) {
			Pfile << std::setprecision(precise) << P(i, 0);
			Pfile << '\t';
			Pfile << std::setprecision(precise) << P(i, 1);
			Pfile << '\t';
			Pfile << std::setprecision(precise) << P(i, 2);
			Pfile << '\n';
		}

		Pfile.close();
	}

	if (Zfile.is_open())
	{
		for (size_t i = 0; i < Z.rows(); i++) {
			for (size_t j = 0; j < Z.cols(); j++) {
				Zfile << std::setprecision(precise) << Z(i, j);

				if (j < Z.cols() - 1)
					Zfile << '\t';
			}

			Zfile << '\n';
		}

		Zfile.close();
	}
}







int main() {



	//TODO: take as input from command line
	//---------------------------------------
	bool simulation = false;
	bool common = true;
	
	int highLowCutOff = 30;
	bool collapseCoding = false;
	bool collapseExon = true;

	//filtering paramaters
	double mafCutoff = 0.05;
	double missingThreshold = 0.2;
	bool onlySNPs = true;
	bool mustPASS = true;


	bool rvs = true;



	///input files
	std::string vcfDir = "C:/Users/Scott/Desktop/RVS-master/example/example_1000snps.vcf";
	std::string infoDir = "C:/Users/Scott/Desktop/RVS-master/example/sampleInfo.txt";
	//std::string bedDir = "";
	std::string bedDir = "C:/Users/Scott/Desktop/RVS-master/example/chr11.bed";
	//---------------------------------------

	//TODO: check to see if file can be opened when another application is using it (excel)
	//TODO: test windows vs unix EOF characters, doesn't seem to work well with windows

	VectorXd Y, G; MatrixXd X, Z, P;
	std::map<int, int> readGroup;
	std::vector<std::vector<int>> interval;

	if (simulation) {
		simulate(X, Y, G, readGroup, P);
		//std::vector<double> pvals = runCommonTest(X, Y, G, readGroup, P, 1000);
		std::vector<double> pvals = runCommonTest(X, Y, G, readGroup, P, 1000, true);

		std::cout << "Common Test p-values\n";
		for (size_t i = 0; i < pvals.size(); i++) {
			std::cout << pvals[i];
			std::cout << '\n';
		}
	}
	else {
		bool valid = parseAndFilter(vcfDir, infoDir, bedDir, 
			highLowCutOff, collapseCoding, collapseExon,
			missingThreshold, onlySNPs, mustPASS, mafCutoff, common,
			X, Y, Z, G, readGroup, P, interval);
		if (!valid) {
			return 0;
		}


		//generateForR(X, Y, Z, G, P, readGroup);

		if (common) {
			//std::vector<double> pvals = runCommonTest(X, Y, G, readGroup, P, 1000);
			std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P);
			//std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P, 1000, true);


			std::cout << "Common Test p-values\n";
			for (size_t i = 0; i < pvals.size(); i++) {
				std::cout << pvals[i];
				std::cout << '\n';
			}
		}
		else {

			std::vector<std::vector<double>> pval = runRareTest(X, Y, G, readGroup, P, 20000, true);
			//std::vector<std::vector<double>> pval = runRareTest(X, Y, G, readGroup, P, Z, 20000, true);

			std::cout << "Rare Test p-values\n";
			std::cout << pval[0][0];
			std::cout << '\t';
			std::cout << pval[0][1];

		}
	}


	//keep console open while debugging
	//TODO: be sure to remove eventually!
	std::cout << "\ndone...>";

	while (true) {}
	return 0;
}