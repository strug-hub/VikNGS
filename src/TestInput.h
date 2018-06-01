#pragma once
#include "Log.h"
#include "Variant.h"
#include "Parser/BED/Interval.h"

#include <vector>
#include <map>

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct TestInput {

	VectorXd Y, G; MatrixXd X, Z, P;
	std::map<int, int> readGroup;
        std::vector<Interval> intervals;
        std::vector<std::vector<int>> collapse;

        std::vector<Variant> variants;
        std::string family;

	inline bool hasCovariates() { return Z.rows() > 0 && Z.cols() > 0; }
        inline int countCovariates() { return Z.cols() -1; }
        inline bool hasIntervals() { return intervals.size() > 0; }
	inline int getNumberOfGroups() { return 1 + (int)G.maxCoeff(); }
};

/**
Creates a TestInput object.
*/
inline TestInput buildTestInput(VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup,
                                std::vector<Interval> intervals, std::string family) {

	TestInput input;
	input.Y = Y;
	input.Z = Z;
	input.G = G;
	input.readGroup = readGroup;
        input.family = family;
        input.intervals = intervals;
	return input;
}

/**
Creates a TestInput object.
*/
inline TestInput buildTestInput(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, MatrixXd &P,
        std::map<int, int> &readGroup, std::vector<Variant> &variants, std::string family) {

        TestInput input;
        input.X = X;
        input.Y = Y;
        input.Z = Z;
        input.G = G;
        input.P = P;

        input.readGroup = readGroup;
        input.variants = variants;
        input.family = family;
        return input;
}

inline TestInput addVariants(TestInput t, std::vector<Variant> &variants, std::vector<std::vector<int>> &collapse) {

        MatrixXd X(variants[0].likelihood.size(), variants.size());
        for (int i = 0; i < variants.size(); i++)
                X.col(i) = variants[i].expectedGenotype;

        MatrixXd P(variants.size(), 3);
        for (int i = 0; i < variants.size(); i++)
                P.row(i) = variants[i].P;

        t.variants = variants;
        t.X = X;
        t.P = P;
        t.collapse = collapse;
        return t;
}
