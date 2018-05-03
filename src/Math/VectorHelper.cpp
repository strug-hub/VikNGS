#include "MathHelper.h"

static const std::string MATH_HELPER = "math helper";


VectorXd concat(std::vector<VectorXd> &v) {

    VectorXd cat = v[0];

    for (int i = 1; i < v.size(); i++) {
        VectorXd temp(cat.rows() + v[i].rows());
        temp << cat, v[i];
        cat = temp;
    }

    return cat;
}


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

		if (std::isnan(X[i]) || std::isnan(Y[i]))
			isNAN[i] = 1;
		else {
			for (int j = 0; j < ncov; j++)
				if (std::isnan(Z(i, j)))
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

		if (std::isnan(Y[i]))
			isNAN[i] = 1;
		else {
			for (int j = 0; j < ncov; j++)
				if (std::isnan(Z(i, j)))
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
		if (std::isnan(X[i]) || std::isnan(Y[i]))
			toRemove[i] = 1;
	}
	return toRemove;
}

VectorXd whereNAN(VectorXd &X) {
	VectorXd toRemove(X.rows());
	for (int i = 0; i < X.rows(); i++) {
		toRemove[i] = 0;

		if (std::isnan(X[i]))
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

double variance(std::vector<VectorXd> v) {

    double mean = average(v);
    double sum = 0;
    double n = 0;

    for (int i = 0; i < v.size(); i++) {

        for (int j = 0; j < v[i].size(); j++) {

            sum += (v[i][j] - mean)*(v[i][j] - mean);
            n++;
        }
    }

    //not this is n, not n-1
    return sum/(n);
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
			if (std::isnan(M(j, i)))
				M(j, i) = 0;

	return  M;
}

VectorXd nanToZero(VectorXd &V) {

	for (size_t i = 0; i < V.rows(); i++)
			if (std::isnan(V[i]))
				V[i] = 0;

	return V;
}
