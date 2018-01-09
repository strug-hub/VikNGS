#pragma once
#include "../RVS.h"

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

		double ybar = average(this->y);

		for (int i = 0; i < size(); i++) {
			xcenter.push_back(this->x[i].array() - this->x[i].mean());
			ycenter.push_back(this->y[i].array() - ybar);
		}

		robustVar = calcRobustVar(P);
	}

	inline double getScore(int i) {	return (ycenter[i].array() * x[i].array()).sum(); }

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
				xrand[j] = xcenter[i][randomInt(0, length - 1)];

			newxboot.push_back(xrand);

		}
		xboot = newxboot;
	}
	inline double bootScore(int i) { return (ycenter[i].array() * xboot[i].array()).sum(); }
	inline double bootVariance(int i) { return ycenter[i].array().pow(2).sum() * variance(xboot[i]); }
	//bootstrapping ------------------------------

};