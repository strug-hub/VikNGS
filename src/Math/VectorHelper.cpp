#include "MathHelper.h"

static const std::string MATH_HELPER = "math helper";


VectorXd concatenate(std::vector<VectorXd> &v) {

    VectorXd cat = v[0];

    for (int i = 1; i < v.size(); i++) {
        VectorXd temp(cat.rows() + v[i].rows());
        temp << cat, v[i];
        cat = temp;
    }

    return cat;
}

MatrixXd concatenate(std::vector<MatrixXd> &m) {

    MatrixXd cat = m[0];

    for (int i = 1; i < m.size(); i++) {
        VectorXd temp(cat.rows() + m[i].rows(), m[i].cols());
        temp << cat, m[i];
        cat = temp;
    }

    return cat;
}


MatrixXd replaceNAN(MatrixXd & M, double value) {
    int nrow = M.rows();
    int ncol = M.cols();

    for (int i = 0; i < nrow; i++)
        for (int j = 0; j < ncol; j++)
                if (std::isnan(M(i, j)))
                    M(i, j) = value;

    return M;
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

VectorXd whereNAN(VectorXd &V1, VectorXd &V2, MatrixXd &M) {

    int nrow = V1.rows();
    int ncol = M.cols();

    VectorXd isNAN(nrow);
	
    for (int i = 0; i < nrow; i++) {
		isNAN[i] = 0;

        if (std::isnan(V1[i]) || std::isnan(V2[i]))
			isNAN[i] = 1;
		else {
            for (int j = 0; j < ncol; j++)
                if (std::isnan(M(i, j)))
					isNAN[i] = 1;
		}
	}

	return isNAN;
}

VectorXd whereNAN(VectorXd &V, MatrixXd &M) {
    int nrow = V.rows();
    int ncol = M.cols();

    VectorXd isNAN(nrow);

    for (int i = 0; i < nrow; i++) {
		isNAN[i] = 0;

        if (std::isnan(V[i]))
			isNAN[i] = 1;
		else {
            for (int j = 0; j < ncol; j++)
                if (std::isnan(M(i, j)))
					isNAN[i] = 1;
		}
	}

	return isNAN;
}

//todo: print warning if many rows are removed due to NA?
VectorXd whereNAN(VectorXd &V1, VectorXd &V2) {
    int nrow = V1.rows();
    VectorXd toRemove(nrow);

    for (int i = 0; i < nrow; i++) {
		toRemove[i] = 0;
        if (std::isnan(V1[i]) || std::isnan(V2[i]))
			toRemove[i] = 1;
	}
	return toRemove;
}

VectorXd whereNAN(VectorXd &V) {
    int nrow = V.rows();
    VectorXd toRemove(nrow);
    for (int i = 0; i < nrow; i++) {
		toRemove[i] = 0;

        if (std::isnan(V[i]))
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
