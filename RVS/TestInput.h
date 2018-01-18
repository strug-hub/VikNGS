#pragma once
#include "Log.h"
#include <vector>
#include <map>

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct TestInput {

	VectorXd Y, G; MatrixXd X, Z, P;
	std::map<int, int> readGroup;
	std::vector<std::vector<int>> interval;

	std::vector<std::string> variantInfo;

	inline bool hasCovariates() { return Z.rows() > 0 && Z.cols() > 0; }
	inline bool isCollapsed() { return interval.size() > 0; }
	inline int getNumberOfGroups() { return 1 + (int)G.maxCoeff(); }
};

/**
Creates a TestInput object.
*/
inline TestInput buildTestInput(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, MatrixXd &P,
	std::map<int, int> &readGroup, std::vector<std::vector<int>> &interval, std::vector<std::string> &info) {

	TestInput input;
	input.X = X;
	input.Y = Y;
	input.Z = Z;
	input.G = G;
	input.P = P;

	input.readGroup = readGroup;
	input.interval = interval;
	input.variantInfo = info;

	return input;
}