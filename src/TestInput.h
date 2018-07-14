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

        VectorXd Y, G; MatrixXd Z;
	std::map<int, int> readGroup;
        std::vector<Interval> intervals;
        std::vector<std::vector<int>> collapse;

        std::vector<Variant> variants;
        std::string family;

	inline bool hasCovariates() { return Z.rows() > 0 && Z.cols() > 0; }
        inline int countCovariates() { return Z.cols() -1; }
        inline bool hasIntervals() { return intervals.size() > 0; }
	inline int getNumberOfGroups() { return 1 + (int)G.maxCoeff(); }

        inline MatrixXd getP(){
            MatrixXd P(variants.size(), 3);
            for (int i = 0; i < variants.size(); i++)
                    P.row(i) = variants[i].P;
            return P;
        }

        inline MatrixXd getX(bool useRegular, bool useTrueGenotypes){

            MatrixXd X(variants[0].likelihood.size(), variants.size());

            if(!useRegular){
            for (int i = 0; i < variants.size(); i++)
                    X.col(i) = variants[i].expectedGenotype;
            }
            else if(useRegular && useTrueGenotypes){
                for (int i = 0; i < variants.size(); i++)
                    X.col(i) = variants[i].trueGenotype;
            }
            else{
                for (int i = 0; i < variants.size(); i++)
                    X.col(i) = variants[i].genotypeCalls;
            }

            return X;
        }

        inline void addCollapse(std::vector<std::vector<int>> c) {
                collapse = c;
        }
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
inline TestInput buildTestInput(VectorXd &Y, MatrixXd &Z, VectorXd &G,
        std::map<int, int> &readGroup, std::vector<Variant> &variants, std::string family) {

        TestInput input;
        input.Y = Y;
        input.Z = Z;
        input.G = G;
        input.readGroup = readGroup;
        input.variants = variants;
        input.family = family;
        return input;
}

inline TestInput addVariants(TestInput t, std::vector<Variant> &variants, std::vector<std::vector<int>> &collapse) {

        t.variants = variants;
        t.collapse = collapse;
        return t;
}


