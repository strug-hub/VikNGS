#pragma once
#include "stdafx.h"
#include "Eigen/Dense"

#include <vector>
#include <random>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;



class TestGroup {
private:
	bool hrg = false;
	std::string id;
	int index;
	VectorXd Xboot;
	VectorXd Xcenter;
	bool isCentered = false;

public:
	VectorXd X;
	VectorXd Y;
	VectorXd Ybar;
	MatrixXd Z;
	std::vector<double> P;

	TestGroup();
	TestGroup(std::string id, int index, bool hrg);
	void bootstrapX();
	void centerX();

	double score();
	double variance(bool rvs);
	double bootScore();
	double bootVariance();

	inline bool isHRG() { return(hrg); }
	inline int length() { return X.rows(); }
};

class TestSet {
	private:
		VectorXd Y;
		VectorXd X;
		MatrixXd Z;
		void calculateYbar();

	public:
		TestSet(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group);

		std::vector<TestGroup> groups;

		inline int length() { return groups.size(); }
		TestGroup & operator [](int i) { return groups[i]; }
};