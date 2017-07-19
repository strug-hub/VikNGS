#include "stdafx.h"
#include "RVS.h"
#include <iostream>

class RareTestObject {
private:
	VectorXd Y;
	MatrixXd Z;

	//filter NAN from Y and Z only
	//this x stays the same after every bootstrap!
	std::vector<VectorXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;

	std::vector<VectorXd> x_;
	std::vector<VectorXd> y_;
	std::vector<MatrixXd> z_;

	std::vector<int> readDepth;
	std::vector<VectorXd> xcenter;
	std::vector<VectorXd> ycenter;
	double robustVar;
	double nhrd;
	double nlrd;
	double nhrd_;
	double nlrd_;

	void filterNAN_yz() {
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(y[i], z[i]);
			x[i] = extractRows(x[i], toRemove, 0);
			y[i] = extractRows(y[i], toRemove, 0);
			z[i] = extractRows(z[i], toRemove, 0);
		}
	}

	void filterNAN() {
		std::vector<VectorXd> fx;
		std::vector<VectorXd> fy;
		std::vector<MatrixXd> fz;
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(x[i], y[i], z[i]);
			fx.push_back(extractRows(x[i], toRemove, 0));
			fy.push_back(extractRows(y[i], toRemove, 0));
			fz.push_back(extractRows(z[i], toRemove, 0));
		}
		x_ = fx;
		y_ = fy;
		z_ = fz;
	}

	void filterNAN_YZ() {
		for (int i = 0; i < size(); i++) {
			VectorXd toRemove = whereNAN(Y, Z);
			Y = extractRows(Y, toRemove, 0);
			Z = extractRows(Z, toRemove, 0);
		}
	}

	void countRD() {
		nhrd = 0;
		nhrd_ = 0;
		nlrd = 0;
		nlrd_ = 0;

		for (int i = 0; i < size(); i++) {
			if (readDepth[i] == 1) {
				nhrd += x[i].rows();
				nhrd_ += x_[i].rows();
			}
			else {
				nlrd += x[i].rows();
				nlrd_ += x_[i].rows();
			}
		}
	}


public:
	RareTestObject(VectorXd &X, VectorXd &Y, MatrixXd &Z,
		std::vector<VectorXd> &x, std::vector<VectorXd> &y, std::vector<MatrixXd> &z,
		std::vector<int> &readDepth, VectorXd &P) {

		this->Y = Y;
		this->Z = Z;
		filterNAN_YZ();

		this->x = x;
		this->y = y;
		this->z = z;
		filterNAN_yz();

		this->readDepth = readDepth;
		this->robustVar = calcRobustVar(P);

		filterNAN();
		countRD();

		VectorXd beta = getBeta(X, Y, Z);
		ycenter = fitModel(beta, y_, z_, "norm");

		for (int i = 0; i < size(); i++)
			this->xcenter.push_back(this->x[i].array() - x_[i].mean());
	}

	inline double getRobustVar() { return robustVar; }
	inline double getYm(int depth) {
		double ym = 0;
		for (int i = 0; i < size(); i++)
			if (readDepth[i] == depth)
				ym += ycenter[i].array().pow(2).sum();

		if(depth == 1)
			ym = ym / nhrd_ * nhrd;
		else
			ym = ym / nlrd_ * nlrd;

		return sqrt(ym);
	}

	inline double getScore() {
		double score = 0;
		for (int i = 0; i < size(); i++)
			score += (ycenter[i].array() * x_[i].array()).sum();
		return score;
	}

	inline bool isHRG(int i) { return readDepth[i] == 1; }
	inline VectorXd getX(int i) { return x[i]; }
	inline int size() { return x.size(); }

	//bootstrapping ------------------------------
	
	void bootstrap() {
		std::vector<VectorXd> xboot;
		int i, j, length;
		int c = 0;
		VectorXd Xnew(Y.rows());

		for (i = 0; i < size(); i++) {
			length = xcenter[i].rows();
			VectorXd xrand(length);
			for (j = 0; j < length; j++) {
				xrand[j] = xcenter[i][generateRandomInteger(0, length - 1)];
				Xnew[c] = xrand[j];
				c++;
			}

			xboot.push_back(xrand);
		}
		
		this->x = xboot;
		filterNAN();
		countRD();

		VectorXd beta = getBeta(Xnew, Y, Z);
		ycenter = fitModel(beta, y_, z_, "norm");
	}

	//bootstrapping ------------------------------

};

std::vector<double> getTestStatistics(std::vector<RareTestObject> &t, bool rvs) {
	int nsnp = t.size();
	int i, j;

	MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

	VectorXd var(nsnp);
	VectorXd ym_hrd(nsnp);
	VectorXd ym_lrd(nsnp);

	for (i = 0; i < nsnp; i++) {
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

	VectorXd e = diagS.eigenvalues().real();
	std::vector<double> eigenvals(e.data(), e.data() + e.size());

	//CAST -> linear test
	double cast = pnorm(score.sum() / sqrt(diagS.sum()));
	//C-alpha -> quadratic test
	double calpha = qfc(eigenvals, score.array().pow(2).sum(), score.rows());

	return{ cast, calpha };
}

std::vector<double> rareTest(std::vector<RareTestObject> &t, int nboot, bool rvs) {
	int i, j, h;
	int nsnp = t.size();

	std::vector<double> tobs = getTestStatistics(t, rvs);
	double linearObs = tobs[0];
	double quadObs = tobs[1];

	//start bootstrapping!

	double bootCount = 0;
	double linearCount = 0;
	double quadCount = 0;
	std::vector<double> tsamp;

	for (h = 0; h < nboot; h++) {
		for (i = 0; i < nsnp; i++)
			t[i].bootstrap();

		tsamp = getTestStatistics(t, false);

		if (tsamp[0] <= linearObs)
			linearCount++;
		if (tsamp[1] <= quadObs)
			quadCount++;

		bootCount++;
	}

	return{ (linearCount + 1) / (bootCount + 1), (quadCount + 1) / (bootCount + 1) };
}

//0.00153997      0.0362593
//first 5
//0.667967        0.940141
//next 5
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

