#pragma once
#include "stdafx.h"
#include "Parser/InputParser.h"

#include <string>
#include <vector>
#include <map>

#include <algorithm>
#include <iostream>  
#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::DiagonalMatrix;

//========================================================
// Interval struct
//========================================================
/*
struct Interval {
	int start = 0;
	int end = 0;
	std::vector<int> indexes;
	std::string chr;

	std::vector<SNP> getInInterval(std::vector<SNP> &snps) {
		std::vector<SNP> in;
		for (size_t i = 0; i < indexes.size(); i++)
			in.push_back(snps[indexes[i]]);

		return in;
	}

	bool addIfIn(SNP snp, int index) {

		if (snp.loc >= start && snp.loc <= end && chr.compare(snp.chr) == 0) {
			indexes.push_back(index);
			return true;
		}
		
		return false;
	}

};
*/

//========================================================
// functions
//========================================================

//VCFParser.cpp
bool parseInput(std::string vcfDir, std::string infoDir, double mafCut, bool common,
	MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd &P);
std::vector<std::string> parseHeader(MemoryMapped &, int &);

//BEDParser.cpp
//std::vector<Interval> getIntervals(std::string);
//void collapseVariants(std::vector<SNP> &, std::vector<Interval> &);

//VectorHelper.cpp
VectorXd extractRows(VectorXd &v, VectorXd &where, double equals);
MatrixXd extractRows(MatrixXd &m, VectorXd &where, double equals);
VectorXd whereNAN(VectorXd &X, VectorXd &Y, MatrixXd &Z); 
VectorXd whereNAN(VectorXd &Y, MatrixXd &Z);
VectorXd whereNAN(VectorXd &X);
VectorXd whereNAN(VectorXd &X, VectorXd &Y);
double average(std::vector<VectorXd> v);
double variance(VectorXd &v);

//StatisticsHelper.cpp
VectorXd getBeta(VectorXd &X, VectorXd &Y, MatrixXd &Z);
std::vector<VectorXd> fitModel(VectorXd &beta, std::vector<VectorXd> &y, std::vector<MatrixXd> &z, std::string distribution);
double calcRobustVar(VectorXd p);
double chiSquareOneDOF(double);
std::vector<double> randomSample(std::vector<double> &, int);
VectorXd CovariateRegression(VectorXd &Y, MatrixXd &Z);
int generateRandomInteger(int from, int to);
double generateRandomDouble(double from, double to);
double randomNormal(double mean, double sd);
MatrixXd covariance(MatrixXd &M);
MatrixXd correlation(MatrixXd &M);
double pnorm(double x);

//CommonTest.cpp
std::vector<double> runCommonTest(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot=0, bool rvs=true);
std::vector<double> runCommonTest(MatrixXd &X, VectorXd &Y, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot=0, bool rvs=true);

//RareTest.cpp
std::vector<std::vector<double>> runRareTest(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, bool rvs = true);
std::vector<std::vector<double>> runRareTest(MatrixXd &X, VectorXd &Y, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, bool rvs = true);

//CompQuadForm.cpp
double qfc(std::vector<double>, double, int);

//Simulation.cpp
void simulate();

//========================================================
// timing functions
//========================================================

#include <time.h>

inline clock_t startTime() {
	clock_t t;
	t = clock();
	return t;
}

inline void endTime(clock_t t, std::string label) {
	t = clock() - t;
	std::cout << "Timer " + label + " took ";
	std::cout << (float)t / CLOCKS_PER_SEC;
	std::cout << " seconds\n";
}

inline double endTime(clock_t t) {
	t = clock() - t;
	return (double)t / CLOCKS_PER_SEC;
}