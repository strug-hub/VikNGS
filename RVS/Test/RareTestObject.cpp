#include "RareTestObject.h"

void RareTestObject::filterNAN() {
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

void RareTestObject::filterNAN_yz() {
	for (int i = 0; i < size(); i++) {
		VectorXd toRemove = whereNAN(y[i], z[i]);
		x[i] = extractRows(x[i], toRemove, 0);
		y[i] = extractRows(y[i], toRemove, 0);
		z[i] = extractRows(z[i], toRemove, 0);
	}
}

void RareTestObject::filterNAN_xy() {
	std::vector<VectorXd> fx;
	std::vector<VectorXd> fy;
	for (int i = 0; i < size(); i++) {
		VectorXd toRemove = whereNAN(x[i], y[i]);
		fx.push_back(extractRows(x[i], toRemove, 0));
		fy.push_back(extractRows(y[i], toRemove, 0));
	}
	x_ = fx;
	y_ = fy;
}

void RareTestObject::filterNAN_YZ() {
	for (int i = 0; i < size(); i++) {
		VectorXd toRemove = whereNAN(Y, Z);
		Y = extractRows(Y, toRemove, 0);
		Z = extractRows(Z, toRemove, 0);
	}
}

void RareTestObject::countRD() {
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

void RareTestObject::covBootstrap() {
	std::vector<VectorXd> xboot;
	int i, j, length;
	int c = 0;
	VectorXd Xnew(Y.rows());

	for (i = 0; i < size(); i++) {
		length = xcenter[i].rows();
		VectorXd xrand(length);
		for (j = 0; j < length; j++) {
			xrand[j] = xcenter[i][randomInt(0, length - 1)];
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

void RareTestObject::noCovBootstrap() {
	std::vector<VectorXd> xboot;
	int i, j, length;

	for (i = 0; i < size(); i++) {
		length = xcenter[i].rows();
		VectorXd xrand(length);
		for (j = 0; j < length; j++)
			xrand[j] = xcenter[i][randomInt(0, length - 1)];

		xboot.push_back(xrand);
	}

	this->x = xboot;
	filterNAN_xy();
	countRD();

	double ybar = average(y_);
	
	for (int i = 0; i < size(); i++) 
		ycenter[i] = y_[i].array() - ybar;
	
}