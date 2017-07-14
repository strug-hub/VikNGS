#include "stdafx.h"
#include "RVS.h"
#include <iostream>

class RareTestObject {
private:
	std::vector<VectorXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;
	std::vector<int> readDepth;
	std::vector<VectorXd> ycenter;
	double robustVar;

	void filterNAN() {
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(x[i], y[i], z[i]);
			x[i] = extractRows(x[i], toRemove, 0);
			y[i] = extractRows(y[i], toRemove, 0);
			z[i] = extractRows(z[i], toRemove, 0);
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

		VectorXd beta = getBeta(X, Y, Z);
		this->ycenter = fitModel(beta, this->y, this->z, "norm");
		for (int i = 0; i < size(); i++)
			this->xcenter.push_back(this->x[i].array() - this->x[i].mean());

		this->robustVar = calcRobustVar(P);
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
/*
std::vector<double> bootstrap(MultiTestSet &ts, int nboot, std::vector<double> &tobs) {
	size_t i;
	int n = ts.length();
	std::vector<double> tBoot;
	int tCountLinear = 0;
	int tCountQuadratic = 0;
	int bootCount = 0;

	MatrixXd diagYmHRG = ts.getYmHRG().asDiagonal();
	MatrixXd diagYmLRG = ts.getYmLRG().asDiagonal();
	MatrixXd X;
	MatrixXd diagS;

	for (size_t h = 0; h < nboot; h++) {
		diagS = MatrixXd::Constant(n, n, 0);

		for (i = 0; i < ts.ngroups(); i++) {
			X = ts.getBootstrapXMatrix(i);
			diagS += diagYmLRG * covariance(X) * diagYmLRG;
		}

	//	tBoot = getTestStatistics(diagS, sample, group, false);

			if (tBoot[0] <= tobs[0])
				tCountLinear++;
			if (tBoot[1] <= tobs[1])
				tCountQuadratic++;

			bootCount++;
		
	}

	return{ (tCountLinear + 1.0) / (bootCount + 1.0), 
		    (tCountQuadratic + 1.0) / (bootCount + 1.0) };
}

std::vector<double> RVSrare(MultiTestSet &ts, int nboot, bool rvs) {
	size_t i, j;
	int n = ts.length();
	MatrixXd diagS = MatrixXd::Constant(n, n, 0);

	VectorXd var = sqrt(ts.getRobustVariance().array());
	MatrixXd diagVar = var.asDiagonal();

	VectorXd YmHRG = ts.getYmHRG();
	MatrixXd diagYmHRG = YmHRG.asDiagonal();

	VectorXd YmLRG = ts.getYmLRG();
	MatrixXd diagYmLRG = YmLRG.asDiagonal();

	MatrixXd X;

	for (i = 0; i < ts.ngroups(); i++) {
		X = ts.getXMatrix(i);

		if (ts.isHRG(i)) {
			if (rvs) {
				MatrixXd varHRG = diagVar.transpose() * correlation(X) * diagVar;
				diagS += diagYmHRG * varHRG * diagYmHRG;
			}
			else
				diagS += diagYmHRG * covariance(X) * diagYmHRG;
		}
		else 
			diagS += diagYmLRG * covariance(X) * diagYmLRG;
	}

	std::vector<double> tobs = getTestStatistics(diagS, ts.getScoreVector());

	std::cout << "\n";
	std::cout << tobs[0];
	std::cout << "\n";
	std::cout << tobs[1];

	return tobs;
	//TODO: bootstrapping!!
}
*/

std::vector<double> rareTest(RareTestObject &t, int nboot, bool rvs) {
	
	std::vector<double> tobs;
	size_t i, j;
	int n = t.size();
	MatrixXd diagS = MatrixXd::Constant(n, n, 0);
	/*
	VectorXd var = sqrt(ts.getRobustVariance().array());
	MatrixXd diagVar = var.asDiagonal();

	VectorXd YmHRG = ts.getYmHRG();
	MatrixXd diagYmHRG = YmHRG.asDiagonal();

	VectorXd YmLRG = ts.getYmLRG();
	MatrixXd diagYmLRG = YmLRG.asDiagonal();

	MatrixXd X;

	for (i = 0; i < ts.ngroups(); i++) {
		X = ts.getXMatrix(i);

		if (ts.isHRG(i)) {
			if (rvs) {
				MatrixXd varHRG = diagVar.transpose() * correlation(X) * diagVar;
				diagS += diagYmHRG * varHRG * diagYmHRG;
			}
			else
				diagS += diagYmHRG * covariance(X) * diagYmHRG;
		}
		else
			diagS += diagYmLRG * covariance(X) * diagYmLRG;
	}

	std::vector<double> tobs = getTestStatistics(diagS, ts.getScoreVector());

	std::cout << "\n";
	std::cout << tobs[0];
	std::cout << "\n";
	std::cout << tobs[1];
	*/
	return tobs;
	//TODO: bootstrapping!!
}


std::vector<std::vector<double>> runRareTest(MatrixXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, MatrixXd P,
	int nboot, bool rvs) {

	int i, j;
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
		RareTestObject t(X_i, Y, Z, x_i, y, z, rd, P_i);

		pvals.push_back(rareTest());

	}

	return pvals;
}

