#include "VectorHelper.h"

static const std::string VECTOR_HELPER = "vector helper";

/*
if (v.rows() != where.rows()) {
std::string message = "Internal error in function extractRows: number of rows in input vector (" +
std::to_string(v.rows()) + ") should match number of rows in condition vector (" +
std::to_string(where.rows()) + ").";

throwError(VECTOR_HELPER, message);
}
*/

VectorXd extractRows(VectorXd &v, VectorXd &where, double equals) {
	
	VectorXd subset(v.rows());
	int c = 0;

	for (int i = 0; i < where.rows(); i++) {
		if (where[i] == equals) {
			subset[c] = v[i];
			c++;
		}
	}

	return subset.block(0, 0, c, 1);
}

MatrixXd extractRows(MatrixXd &m, VectorXd &where, double equals) {

	MatrixXd subset(m.rows(), m.cols());
	int c = 0;
	int i, j;

	for (i = 0; i < where.rows(); i++) {
		if (where[i] == equals) {
			for (j = 0; j < m.cols(); j++)
				subset(c, j) = m(i, j);
			c++;
		}
	}

	return subset.block(0, 0, c, m.cols());
}

VectorXd whereNAN(VectorXd &X, VectorXd &Y, MatrixXd &Z) {

	int nobs = Y.rows();
	int ncov = Z.cols();	

	VectorXd isNAN(nobs);
	
	for (int i = 0; i < nobs; i++) {
		isNAN[i] = 0;

		if (isnan(X[i]) || isnan(Y[i]))
			isNAN[i] = 1;
		else {
			for (int j = 0; j < ncov; j++)
				if (isnan(Z(i, j)))
					isNAN[i] = 1;
		}
	}

	return isNAN;
}

VectorXd whereNAN(VectorXd &Y, MatrixXd &Z) {
	int nobs = Y.rows();
	int ncov = Z.cols();

	VectorXd isNAN(nobs);

	for (int i = 0; i < nobs; i++) {
		isNAN[i] = 0;

		if (isnan(Y[i]))
			isNAN[i] = 1;
		else {
			for (int j = 0; j < ncov; j++)
				if (isnan(Z(i, j)))
					isNAN[i] = 1;
		}
	}

	return isNAN;
}

VectorXd whereNAN(VectorXd &X, VectorXd &Y) {
	int nobs = X.rows();
	int ncov = Y.cols();
	VectorXd toRemove(nobs);

	for (int i = 0; i < nobs; i++) {
		toRemove[i] = 0;
		if (isnan(X[i]) || isnan(Y[i]))
			toRemove[i] = 1;
	}
	return toRemove;
}

VectorXd whereNAN(VectorXd &X) {
	VectorXd toRemove(X.rows());
	for (int i = 0; i < X.rows(); i++) {
		toRemove[i] = 0;

		if (isnan(X[i]))
			toRemove[i] = 1;
		else 
			toRemove[i] = 0;
	}

	return toRemove;
}

double average(std::vector<VectorXd> v) {
	double sum = 0;
	double n = 0;
	for (int i = 0; i < v.size(); i++) {
		n += v[i].size();
		sum += v[i].sum();
	}
	return sum / n;
}

double variance(VectorXd &v) {

	double mean = v.mean();
	double var = (v.array() - mean).pow(2).sum();

	return var/(v.rows()-1);
}

double center(VectorXd &v) {

	double mean = v.mean();
	double var = (v.array() - mean).pow(2).sum();

	return var / (v.rows() - 1);
}

MatrixXd nanToZero(MatrixXd &M) {

	for (size_t i = 0; i < M.cols(); i++) 
		for (size_t j = 0; j < M.rows(); j++) 
			if (isnan(M(j, i)))
				M(j, i) = 0;

	return  M;
}

VectorXd nanToZero(VectorXd &V) {

	for (size_t i = 0; i < V.rows(); i++)
			if (isnan(V[i]))
				V[i] = 0;

	return V;
}