#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <math.h> 

class CommonTestObject {
private:
	//NANs have been removed
	std::vector<VectorXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;
	std::vector<int> readDepth;
	std::vector<VectorXd> ycenter;
	double robustVar;

	void filterNAN_xyz() {
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(x[i], y[i], z[i]);
			x[i] = extractRows(x[i], toRemove, 0);
			y[i] = extractRows(y[i], toRemove, 0);
			z[i] = extractRows(z[i], toRemove, 0);
		}
	}

	void filterNAN_xy() {
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(x[i], y[i]);
			x[i] = extractRows(x[i], toRemove, 0);
			y[i] = extractRows(y[i], toRemove, 0);
		}
	}

public:
	CommonTestObject(VectorXd &X, VectorXd &Y, MatrixXd &Z,
		std::vector<VectorXd> &x, std::vector<VectorXd> &y, std::vector<MatrixXd> &z,
		std::vector<int> &readDepth, VectorXd &P) {

		this->x = x;
		this->y = y;
		this->z = z;
		this->readDepth = readDepth;

		filterNAN_xyz();

		VectorXd beta = getBeta(X, Y, Z);
		ycenter = fitModel(beta, this->y, this->z, "norm");
		for (int i = 0; i < size(); i++)
			xcenter.push_back(this->x[i].array() - this->x[i].mean());

		robustVar = calcRobustVar(P);
	}

	CommonTestObject(std::vector<VectorXd> &x, std::vector<VectorXd> &y, 
		std::vector<int> &readDepth, VectorXd &P) {

		this->x = x;
		this->y = y;
		this->readDepth = readDepth;

		filterNAN_xy();

		double meanY = average(this->y);
		
		for (int i = 0; i < size(); i++) {
			xcenter.push_back(this->x[i].array() - this->x[i].mean());
			ycenter.push_back(this->y[i].array() - meanY);
		}

		robustVar = calcRobustVar(P);
	}

	inline double getScore(int i) { return (ycenter[i].array() * x[i].array()).sum(); }

	inline double getVariance(int i, bool rvs) {
			double var = ycenter[i].array().pow(2).sum();
			if (rvs && readDepth[i] == 1)
				return var * robustVar;
			else
				return var * variance(x[i]);
	}

	inline int size() { return x.size(); }

	//bootstrapping ------------------------------
	std::vector<VectorXd> xcenter;
	std::vector<VectorXd> xboot;

	void bootstrap() {
		std::vector<VectorXd> newxboot;
		int i, j, length;

		for (i = 0; i < size(); i++) {
			length = x[i].rows();
			VectorXd xrand(length);
			for (j = 0; j < length; j++) 
				xrand[j] = xcenter[i][generateRandomInteger(0, length - 1)];
			
			newxboot.push_back(xrand);

		}
		xboot = newxboot;
	}
	inline double bootScore(int i) { return (ycenter[i].array() * xboot[i].array()).sum(); }
	inline double bootVariance(int i) { return ycenter[i].array().pow(2).sum() * variance(xboot[i]); }
	//bootstrapping ------------------------------
	
};

/*
RVS using asymptotic distribution for score test statistic. This functions includes
two association tests which differ in the estimation of the variance of score. It uses
var_case(E(G | D)) when rvs = false and var_case(G) when rvs = true.

Reference: http://www.ncbi.nlm.nih.gov/pubmed/22570057
Reference: http://www.ncbi.nlm.nih.gov/pubmed/24733292
@param CommonTestObject Test object.
@param rvs Indicates how to calculate variance
@return p-value
*/
double commonAsymptotic(CommonTestObject t, bool rvs) {

	double variance = 0;
	double score = 0;

	for (int i = 0; i < t.size(); i++) {
		score += t.getScore(i);
		variance += t.getVariance(i, rvs);
	}
	double xxx = pow(score, 2);
	return chiSquareOneDOF(pow(score, 2) / variance);
}


/*
Uses RVS to test associaton by bootstrap, given phenotype, expected values of genotypes,
estimated genotype frequency and number of bootstrap iterations.

@param snps Vector of SNPs.
@param sample Vector with sample information.
@param nboot Number of bootstrap iterations.
@param earlyStop Stop bootstrapping with p-value is large.
@param rvs Indicates how to calculate variance of score statistic.
@return Vector of p-values for the SNPs.
*/
double commonBootstrap(CommonTestObject t, int nboot, bool earlyStop, bool rvs) {
	int i;
	double variance = 0;
	double score = 0;

	for (int i = 0; i < t.size(); i++) {
		score += t.getScore(i);
		variance += t.getVariance(i, rvs);
	}

	double tobs = pow(score, 2) / variance;
	double tcount = 0;
	double average = 0;

	for (int n = 0; n < nboot; n++) {
		variance = 0;
		score = 0;

		t.bootstrap();

		for (i = 0; i < t.size(); i++) {
			score += t.bootScore(i);
			variance += t.bootVariance(i);
		}

		if (tobs <= pow(score, 2) / variance)
			tcount++;

		average += score;
	}

	return (tcount + 1) / (nboot + 1);

	/*
	if (earlyStop && bootcount > 1000) {
	pstar = a / ((bootcount + c)*(1 + c));
	if (bootScoreCount / bootcount > pstar) {
	std::cout << "early stop";
	break;
	}
	*/
}

std::vector<double> runCommonTest(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, bool rvs) {

	int i,j;
	std::vector<double> pvals;

	std::vector<MatrixXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;
	std::vector<int> rd;

	int ngroups = 1 + (int)G.maxCoeff();

	for (i = 0; i < ngroups; i++) {
		x.push_back(extractRows(X, G, i));
		y.push_back(extractRows(Y, G, i));
		z.push_back(extractRows(Z, G, i));
		rd.push_back(readGroup[i]);
	}

	for (i = 0; i < X.cols(); i++) {

		std::vector<VectorXd> x_i;
		for (j = 0; j < ngroups; j++) 
			x_i.push_back(x[j].col(i));

		VectorXd X_i = X.col(i);
		VectorXd P_i = P.row(i);
		CommonTestObject t(X_i, Y, Z, x_i, y, z, rd, P_i);

		if (nboot > 0)
			pvals.push_back(commonBootstrap(t, nboot, false, rvs));
		else
			pvals.push_back(commonAsymptotic(t, rvs));
	}

	return pvals;

}

std::vector<double> runCommonTest(MatrixXd &X, VectorXd &Y, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, bool rvs) {

	int i, j;
	std::vector<double> pvals;

	std::vector<MatrixXd> x;
	std::vector<VectorXd> y;
	std::vector<int> rd;

	int ngroups = 1 + (int)G.maxCoeff();

	for (i = 0; i < ngroups; i++) {
		x.push_back(extractRows(X, G, i));
		y.push_back(extractRows(Y, G, i));
		rd.push_back(readGroup[i]);
	}

	for (i = 0; i < X.cols(); i++) {

		std::vector<VectorXd> x_i;
		for (j = 0; j < ngroups; j++)
			x_i.push_back(x[j].col(i));

		VectorXd P_i = P.row(i);
		CommonTestObject t(x_i, y, rd, P_i);

		if (nboot > 0)
			pvals.push_back(commonBootstrap(t, nboot, false, rvs));
		else
			pvals.push_back(commonAsymptotic(t, rvs));
	}

	return pvals;

}
