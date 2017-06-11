#include "stdafx.h"
#include "RVS.h"
#include "TestSet.h"

void TestSet::calculateYbar() {
	VectorXd beta = CovariateRegression(Y, Z);

	for (size_t i = 0; i < groups.size(); i++)
		groups[i].Ybar = fitModel(beta, groups[i].Z);
}

TestSet::TestSet(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group) {

	size_t i, j, k;
	int index, length, c;
	bool flag;
	std::vector<TestGroup> grps(group.size());

	VectorXd bigx(sample.size());
	VectorXd bigy(sample.size());
	MatrixXd bigz(sample.size(), sample[0].covariates.size() + 1);
	int d = 0;

	for (i = 0; i < group.size(); i++) {
		index = group[i].groupIndex;
		length = group[i].index.size();
		grps[index] = TestGroup(group[i].ID, index, group[i].hrg);

		VectorXd x(length);
		VectorXd y(length);
		MatrixXd z(length, sample[0].covariates.size() + 1);

		c = 0;
		for (j = 0; j < length; j++) {
			index = group[i].index[j];

			if (!isnan(snp.EG[index]) && !isnan(sample[index].y)) {

				x[c] = snp.EG[index];
				bigx[d] = snp.EG[index];
				y[c] = sample[index].y;
				bigy[d] = sample[index].y;

				flag = false;

				for (size_t k = 0; k < sample[index].covariates.size(); k++) {
					if (isnan(sample[index].covariates[k])) {
						flag = true;
						break;
					}
					z(c, k + 1) = sample[index].covariates[k];
					bigz(d, k + 1) = sample[index].covariates[k];
				}

				if (!flag) {
					z(c, 0) = 1;
					bigz(d, 0) = 1;
					c++;
					d++;
				}
			}
		}
		index = group[i].groupIndex;
		grps[index].X = x.block(0, 0, c, 1);
		grps[index].Y = y.block(0, 0, c, 1);;
		grps[index].Z = z.block(0, 0, c, sample[0].covariates.size() + 1);
		grps[index].P = snp.p;
	}

	groups = grps;
	X = bigx.block(0, 0, d, 1);
	Y = bigy.block(0, 0, d, 1);;
	Z = bigz.block(0, 0, d, sample[0].covariates.size() + 1);

	calculateYbar();
}

void TestGroup::centerX() {
	VectorXd xmid(length());

	double mean = X.mean();

	for (size_t i = 0; i < length(); i++)
		xmid[i] = X[i] - mean;

	Xcenter = xmid;
	isCentered = true;
}

void TestGroup::bootstrapX() {
	if (!isCentered)
		centerX();

	VectorXd xrand(length());

	for (size_t i = 0; i < length(); i++) 
		xrand[i] = Xcenter[generateRandomNumber(0, length() - 1)];

	Xboot = xrand;
}

TestGroup::TestGroup(){}
TestGroup::TestGroup(std::string id, int index, bool hrg) {
	this->id = id;
	this->index = index;
	this->hrg = hrg;
}

double TestGroup::score() {
	double score = 0;
	for (size_t i = 0; i < length(); i++)
		score += (Y[i] - Ybar[i]) * X[i];

	return score;
}

double TestGroup::variance(bool rvs) {
	double variance = 0;
	if (rvs && isHRG()) {
		for (size_t i = 0; i < length(); i++)
			variance += pow(Y[i] - Ybar[i], 2);
		variance = variance * calcRobustVar(P);
	}
	else{
		for (size_t i = 0; i < length(); i++)
			variance += pow(Y[i] - Ybar[i], 2);
		variance = variance * var(X);
	}

	return variance;
}

double TestGroup::bootScore() {
	double score = 0;
	for (size_t i = 0; i < length(); i++)
		score += (Y[i] - Ybar[i]) * Xboot[i];

	return score;
}

double TestGroup::bootVariance() {
	double variance = 0;
	for (size_t i = 0; i < length(); i++)
		variance += pow(Y[i] - Ybar[i], 2);
	variance = variance * var(Xboot);

	return variance;
}