#pragma once

#include "Parser/InputParser.h"
#include "Request.h"
#include "Variant.h"
#include "TestInput.h"
#include "Output/OutputHandler.h"

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


//========================================================
// functions
//========================================================

std::vector<Variant> startVikNGS(Request req);
std::vector<Variant> runTest(TestInput &input, Request &req);

//InputParser
TestInput parseSampleLines(Request req);
std::vector<Variant> processVCF(TestInput input, Request req);
TestInput parseInfo(Request req);

//CommonTest.cpp
std::vector<Variant> runCommonTest(Request &req, TestInput &input);
//RareTest.cpp
std::vector<Variant> runRareTest(Request req, TestInput input);


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
