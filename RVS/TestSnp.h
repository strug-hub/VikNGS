#pragma once
#include "TestGroup.h"

class TestSnp {
private:

	std::vector<VectorXd> X;
	std::vector<VectorXd> Y;
	std::vector<MatrixXd> Z;
	double robustVar;
	std::vector<bool> hrg;


	VectorXd getBeta(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group);

	//for rare tests:
	//after filtering for NAN in X, Y, and Z
	int nhrg;
	int nlrg;
	//after filtering for NAN in Y and Z
	int nhrg_filterz;
	int nlrg_filterz;

public:
	TestSnp() {};
	TestSnp(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group, bool rare);
	TestSnp(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group)
		: TestSnp(snp, sample, group, false) {}

	std::vector<TestHRG> groupHRG;
	std::vector<TestLRG> groupLRG;

	inline double getRobustVariance() { return robustVar; }
	inline int length() { return groupHRG.size() + groupLRG.size(); }
	TestGroup & operator [](int i) { return get(i); }

	inline TestGroup & get(int i) {
		int index = i - groupHRG.size();
		if (index < 0)
			return groupHRG[i];
		else
			return groupLRG[index];
	}

	//rare tests:
	double getScore();
	double getBootstrapScore();

	double getYmLRG();
	double getYmHRG();
	inline VectorXd getX(int index) { return get(index).getX_filterz(); }
	inline bool isHRG(int index) { return index < groupHRG.size(); }

};

