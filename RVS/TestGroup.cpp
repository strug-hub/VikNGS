#include "stdafx.h"
#include "RVS.h"
#include "TestGroup.h"


//---------------------------------
//	TestGroup
//---------------------------------

TestGroup::TestGroup(SNP &snp, std::vector<Sample> &sample, Group &group, bool rare) {
	int nobs = group.index.size();
	int ncov = sample[0].covariates.size();

	int index;
	size_t i,j;

	VectorXd x(nobs);
	VectorXd y(nobs);
	MatrixXd z(nobs, ncov + 1);
	
	for (i = 0; i < nobs; i++) {
		index = group.index[i];

		x[i] = snp.EG[index];
		y[i] = sample[index].y;
		for (j = 0; j < ncov; j++){
			z(i, j + 1) = sample[index].covariates[j];
			z(i, 0) = 1;
		}
	}

	this->X_nofilter = x;
	this->Y_nofilter = y;
	this->Z_nofilter = z;

	filterNAN();
	if (rare)
		filterNAN_z();

	centerX();
	this->robustVar = calcRobustVar(snp.p);

}

void TestGroup::filterNAN() {
	size_t i, j;
	int nobs = X_nofilter.rows();
	int ncov = Z_nofilter.cols();
	std::vector<bool> toRemove(nobs, false);

	for (i = 0; i < nobs; i++) {
		if (isnan(X_nofilter[i]) || isnan(Y_nofilter[i]))
			toRemove[i] = true;

		else {
			for (j = 0; j < ncov; j++)
				if (isnan(Z_nofilter(i,j)))
					toRemove[i] = true;
		}
	}

	VectorXd x(nobs);
	VectorXd y(nobs);
	MatrixXd z(nobs, ncov + 1);

	int c = 0;
	for (i = 0; i < nobs; i++) {
		if (!toRemove[i]) {
			x[c] = X_nofilter[i];
			y[c] = Y_nofilter[i];

			for (j = 0; j < ncov; j++)
				z(c, j) = Z_nofilter(i, j);

			c++;
		}
	}

	this->X = x.block(0, 0, c, 1);
	this->Y = y.block(0, 0, c, 1);
	this->Z = z.block(0, 0, c, ncov);
}

void TestGroup::filterNAN_z() {
	size_t i, j;
	int nobs = X_nofilter.rows();
	int ncov = Z_nofilter.cols();
	std::vector<bool> toRemove(nobs, false);

	for (i = 0; i < nobs; i++) {
		if (isnan(Y_nofilter[i]))
			toRemove[i] = true;

		else {
			for (j = 0; j < ncov; j++)
				if (isnan(Z_nofilter(i, j)))
					toRemove[i] = true;
		}
	}

	VectorXd x(nobs);
	VectorXd y(nobs);
	MatrixXd z(nobs, ncov + 1);

	int c = 0;
	for (i = 0; i < nobs; i++) {
		if (!toRemove[i]) {
			x[c] = X_nofilter[i];
			y[c] = Y_nofilter[i];

			for (j = 0; j < ncov; j++)
				z(c, j) = Z_nofilter(i, j);

			c++;
		}
	}

	this->X_filterz = x.block(0, 0, c, 1);
	this->Y_filterz = y.block(0, 0, c, 1);
	this->Z_filterz = z.block(0, 0, c, ncov);

	this->X_filterz = nanToZero(X_filterz);
	this->Y_filterz = nanToZero(Y_filterz);
	this->Z_filterz = nanToZero(Z_filterz);
}

void TestGroup::centerX() {
	Xcenter = (X.array() - X.mean());
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

double TestGroup::bootVariance() {
	double variance = Ycenter.array().pow(2).sum();
	return variance*var(Xboot);
}

void TestGroup::fitModel(VectorXd &beta, std::string distribution) {

	VectorXd meanValue = Z * beta;

	if (distribution == "binom")
		meanValue = 1 / (1 + exp(-meanValue.array()));

	Ycenter = Y - meanValue;
}


//---------------------------------
//	TestHRG
//---------------------------------

double TestHRG::variance(bool rvs) {
	double variance = Ycenter.array().pow(2).sum();

	if (rvs)
		return variance * robustVar;
	else
		return variance = variance * var(X);
}

//---------------------------------
//	TestLRG
//---------------------------------

double TestLRG::variance(bool rvs) {
	double variance = Ycenter.array().pow(2).sum();
		return variance = variance * var(X);
}
