#include "stdafx.h"
#include "RVS.h"
#include <iostream>

class RareTestObject {
private:
	std::vector<VectorXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;

	//filter all NAN
	std::vector<VectorXd> x__;
	std::vector<VectorXd> y__;
	std::vector<MatrixXd> z__;

	//filter NAN from Y and Z only
	std::vector<VectorXd> x_;
	std::vector<VectorXd> y_;
	std::vector<MatrixXd> z_;

	std::vector<int> readDepth;
	std::vector<VectorXd> ycenter;
	double robustVar;
	double nhrd_;
	double nlrd_;
	double nhrd__;
	double nlrd__;

	void filterNAN() {
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(x[i], y[i], z[i]);
			x__.push_back(extractRows(x[i], toRemove, 0));
			y__.push_back(extractRows(y[i], toRemove, 0));
			z__.push_back(extractRows(z[i], toRemove, 0));
		}
	}

	void filterNAN_ZYOnly() {
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(y[i], z[i]);
			x_.push_back(extractRows(x[i], toRemove, 0));
			y_.push_back(extractRows(y[i], toRemove, 0));
			z_.push_back(extractRows(z[i], toRemove, 0));
		}
	}

	void countRD() {
		for (int i = 0; i < size(); i++) {
			if (readDepth[i] == 1) {
				nhrd_ += x_[i].rows();
				nhrd__ += x__[i].rows();
			}
			else {
				nlrd_ += x_[i].rows();
				nlrd__ += x__[i].rows();
			}
		}
	}

public:
	RareTestObject(VectorXd &X, VectorXd &Y, MatrixXd &Z,
		std::vector<VectorXd> &x, std::vector<VectorXd> &y, std::vector<MatrixXd> &z,
		std::vector<int> &readDepth, VectorXd &P) {

		this->x = x;
		this->y = y;
		this->z = z;
		this->readDepth = readDepth;

		filterNAN();
		filterNAN_ZYOnly();
		countRD();

		VectorXd beta = getBeta(X, Y, Z);
		this->ycenter = fitModel(beta, y__, z__, "norm");
		
		for (int i = 0; i < size(); i++)
			this->xcenter.push_back(this->x[i].array() - this->x[i].mean());

		this->robustVar = calcRobustVar(P);
	}

	inline double getRobustVar() { return robustVar; }
	inline double getYm(int depth) {
		double ym = 0;
		for (int i = 0; i < size(); i++)
			if (readDepth[i] == depth)
				ym += ycenter[i].array().pow(2).sum();

		if(depth ==1)
			ym = ym / nhrd__ * nhrd_;
		else
			ym = ym / nlrd__ * nlrd_;

		return sqrt(ym);
	}

	inline double getScore() {
		double score = 0;
		for (int i = 0; i < size(); i++)
			score += (ycenter[i].array() * x__[i].array()).sum();
		return score;
	}

	inline bool isHRG(int i) { return readDepth[i] == 1; }
	inline VectorXd getX(int i) { return x_[i]; }
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
RVS analysis with rare variants for one group. Get p-values from RVS
with CAST and C-alpha (resampling: bootstrap; variance estimate : robust)

@param snps Vector of SNPs.
@param sample Vector with sample information.
@param group Vector with group information.
@param nboot Number of bootstrap iterations.
@param rvs Indicates how to calculate variance of score statistic.
@param njoint Number of SNPs grouping together for one test, default is 5.

checkHomMatrix parameters:
@param hom 1 or 2; 1 means making changes with the 1st non-NAN element, 2 means making changes with a random non-NAN element
@param multiplier Value 1 or 2; 1 is dividing by 2 and 2 is multiplying by 2.

@return Vector with two p-values. First element is linear p-value (CAST), second is quadratic p-value (C-alpha).

std::vector<double> RVSrare(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group, int nboot, bool rvs, int njoint, int hom, int multiplier) {
	
	SNP snp;
	//std::vector<double> tobs = testStatistic(snps, sample, group, rvs);
	std::vector<double> tsamp;
	std::vector<SNP> samples = cloneAll(snps);

	//SNPs > groups > EG values
	std::vector<std::vector<std::vector<double>>> X;
	std::vector<double> temp;

	double tcountLinear = 0;
	double tcountQuadratic = 0;
	double bootcount = 0;

	size_t i, j, k, l, h;

	//create vectors to sample from (x - xbar)
	for (i = 0; i < snps.size(); i++) {
		snp = snps[i];
		std::vector<std::vector<double>> xsnp;

		for (j = 0; j < group.size(); j++) {

			std::vector<double> xeg;
			double xbar = meanX(snp, group[j]);

			for (l = 0; l < group[j].index.size(); l++) {
				k = group[j].index[l];
				if (!isnan(snp.EG[k]))
					xeg.push_back(snp.EG[k] - xbar);
				//TODO: ask if NAs should be sampled?? (they aren't in common test)
				else
					xeg.push_back(snp.EG[k]);
			}
			xsnp.push_back(xeg);
		}
		X.push_back(xsnp);
	}
	
	//start bootstrapping
	for (h = 0; h < nboot; h++) {
		
		for (i = 0; i < snps.size(); i++) {
			for (j = 0; j < group.size(); j++) {
				temp = randomSample(X[i][j], X[i][j].size());
				
				for (l = 0; l < group[j].index.size(); l++) {
					k = group[j].index[l];
					samples[i].EG[k] = temp[l];
				}
			}
		}
		
	//	tsamp = testStatistic(samples, sample, group, false);

	//	if (tsamp[0] <= tobs[0])
	//		tcountLinear++;
	//	if (tsamp[1] <= tobs[1])
	//		tcountQuadratic++;

		bootcount++;

	}

	return{ (tcountLinear+1)/(bootcount+1), (tcountQuadratic+1)/(bootcount+1) };
}
*/
std::vector<double> getTestStatistics(MatrixXd &diagS, VectorXd &score) {

	VectorXd e = diagS.eigenvalues().real();
	std::vector<double> eigenvals(e.data(), e.data() + e.size());

	//CAST -> linear test
	double cast = pnorm(score.sum() / sqrt(diagS.sum()));
	//C-alpha -> quadratic test
	double calpha = qfc(eigenvals, score.array().pow(2).sum(), score.rows());

	return{ cast, calpha };
}


std::vector<double> rareTest(std::vector<RareTestObject> &t, int nboot, bool rvs) {

	int i, j;
	int nsnp = t.size();
	MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

	VectorXd var(nsnp);
	VectorXd ym_hrd(nsnp);
	VectorXd ym_lrd(nsnp);

	for (i = 0; i < nsnp; i++){
		var[i] = sqrt(t[i].getRobustVar());
		ym_hrd[i] = t[i].getYm(1);
		ym_lrd[i] = t[i].getYm(0);
	}

	MatrixXd diagVar = var.asDiagonal();
	MatrixXd diagYm_hrd = ym_hrd.asDiagonal();
	MatrixXd diagYm_lrd = ym_lrd.asDiagonal();

	for (i = 0; i < t[0].size(); i++) {
		VectorXd x_0 = t[0].getX(i);
		MatrixXd x(x_0.rows(), nsnp);
		x.col(0) = x_0;

		for (j = 1; j < nsnp; j++)
			x.col(j) = t[j].getX(i);

		if (t[0].isHRG(i)) {
			if (rvs) {
				MatrixXd var_hrd = diagVar.transpose() * correlation(x) * diagVar;
				diagS += diagYm_hrd * var_hrd * diagYm_hrd;
			}
			else
				diagS += diagYm_hrd * covariance(x) * diagYm_hrd;
		}
		else
			diagS += diagYm_lrd * covariance(x) * diagYm_lrd;
	}

	VectorXd score(nsnp);
	for (i = 0; i < nsnp; i++)
		score[i] = t[i].getScore();

	std::vector<double> tobs = getTestStatistics(diagS, score);

	return tobs;
	//TODO: bootstrapping!!
}


std::vector<std::vector<double>> runRareTest(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, bool rvs) {

	int i, j, l;
	std::vector<std::vector<double>> pvals;

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

	std::vector<RareTestObject> t;

	
	//for (i = 0; i < X.cols(); i++) {
	
	for (i = 0; i < 5; i++) {
		std::vector<VectorXd> x_i;
		for (j = 0; j < ngroups; j++)
			x_i.push_back(x[j].col(i));

		VectorXd X_i = X.col(i);
		VectorXd P_i = P.row(i);
		RareTestObject t_i(X_i, Y, Z, x_i, y, z, rd, P_i);
		t.push_back(t_i);
	}
	
	pvals.push_back(rareTest(t, nboot, rvs));
	return pvals;
}

