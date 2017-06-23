#pragma once
#include "TestSet.h"

class MultiTestSet {
	public:
		std::vector<TestSet> testset;

		MultiTestSet(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group);

		VectorXd getRobustVariance();

		VectorXd getYmHRG();
		VectorXd getYmLRG();

		VectorXd getScoreVector();
		VectorXd getBootstrapScoreVector();

		inline int length() { return testset.size(); }

		inline int ngroups() { return testset[0].length(); }
		inline bool isHRG(int index) { return testset[0].isHRG(index); }
		MatrixXd getXMatrix(int index);
		MatrixXd getBootstrapXMatrix(int index);

		TestSet & operator [](int i) { return testset[i]; }
};