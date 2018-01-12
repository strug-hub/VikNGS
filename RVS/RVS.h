#pragma once
#include "stdafx.h"
#include "Parser/InputParser.h"
#include "Request.h"

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>  

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Vector3d;
using Eigen::DiagonalMatrix;

struct SimulationRequestGroup {
	int n; //The number of samples in this group
	bool isHrg;
	bool isCase;
	double meanDepth;
	double sdDepth;

	//for debugging
	void print() {
		std::cout << "\tn = " + std::to_string(n) + "\n";
		std::cout << "\tread depth = ";
		std::cout << (isHrg ? "high\n" : "low\n");
		std::cout << "\tgroup = ";
		std::cout << (isCase ? "case\n" : "control\n");

		std::cout << "\tmean depth = " + std::to_string(meanDepth) + "\n";
		std::cout << "\tdepth sd = " + std::to_string(sdDepth) + "\n";
	}
};

struct SimulationRequest {
	int npop; //The number of population
	double prevalence; //A decimal between[0, 1], prevalence rate of the disease.

	int nsnp;  //Integer.The number of variants or bases.
	double me; //The mean error rate of sequencing.
	double sde;  //The standard deviation for the error rate.

	double oddsRatio;  //Under H0
	double maf;

	std::vector<SimulationRequestGroup> groups;

	std::string test;
	bool useBootstrap;
	int nboot;

	//for debugging
	void print() {
		std::cout << "npop = " + std::to_string(npop) + "\n";
		std::cout << "prevalence = " + std::to_string(prevalence) + "\n";
		std::cout << "nsnp = " + std::to_string(nsnp) + "\n";
		std::cout << "me = " + std::to_string(me) + "\n";
		std::cout << "sde = " + std::to_string(sde) + "\n";
		std::cout << "oddsRatio = " + std::to_string(oddsRatio) + "\n";
		std::cout << "MAF = " + std::to_string(maf) + "\n";

		for (int i = 0; i < groups.size(); i++) {
			std::cout << "group " + std::to_string(i) + ":\n";
			groups[i].print();
		}
	}
};




//========================================================
// functions
//========================================================

SimulationRequest newSimulationRequest(std::string npop, std::string prevalence,
	std::string nsnp, std::string me, std::string sde, std::string oddsRatio,
	std::string maf, std::vector<SimulationRequestGroup> groups,
	std::string test, bool useBootstrap, std::string nboot);

SimulationRequestGroup newSimulationRequestGroup(int groupID, std::string n, std::string isCase,
	std::string isHrg, std::string meanDepth, std::string sdDepth);

std::vector<double> startSimulation (SimulationRequest req);
std::vector<double> startVikNGS(Request req);

//InputParser
std::vector<VCFLine> parseAndFilter(Request req, MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, 
	std::map<int, int> &readGroup, MatrixXd &P, std::vector<std::vector<int>> & interval);

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
int randomInt(int from, int to);
double randomDouble(double from, double to);
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
std::vector<double> runRareTest(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, std::string test = "calpha", int collapseNumber = 5, bool rvs = true);
std::vector<double> runRareTest(MatrixXd &X, VectorXd &Y, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, std::string test = "calpha", int collapseNumber = 5, bool rvs = true);

//CompQuadForm.cpp
double qfc(std::vector<double>, double, int);

//Simulation.cpp
void simulate(SimulationRequest req, MatrixXd &X, VectorXd &Y, VectorXd &G, std::map<int, int> &readGroup, MatrixXd &P);

//OutputHandler.cpp
void outputPvals(std::vector<double> pvalues, std::string outputDir);

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
