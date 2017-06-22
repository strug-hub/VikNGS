#pragma once
#include "TestSet.h"

class MultiTestSet {
public:
	std::vector<TestSet> testset;

	MultiTestSet(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group);

	VectorXd getRobustVariance();
	std::vector <double> calculateYbar();

	VectorXd getYm(bool hrg);

	inline int length() { return testset.size(); }
	inline int ngroups() { return testset[0].length(); }
	inline bool isHRG(int index) { return testset[0].isHRG(index); }

	MatrixXd getCovariateMatrix(int index);
	MatrixXd getCorrelationMatrix(int index);

	TestSet & operator [](int i) { return testset[i]; }
};