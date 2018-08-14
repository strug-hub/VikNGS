#include "Math.h"

static const std::string ERROR_SOURCE = "VECTOR_HELPER";

/**
Finds the mean of the values in a vector of vectors.

@param v Vector of vectors which is used to calculate the mean.
@return The mean value of all the vectors.
*/
double average(std::vector<VectorXd> v) {
    double sum = 0;
    double n = 0;
    for (size_t i = 0; i < v.size(); i++) {
        n += v[i].size();
        sum += v[i].sum();
    }
    return sum / n;
}

/**
Finds the variance of the values in a vector of vectors.
Note: uses n, not n-1 for final division

@param v Vector of vectors which is used to calculate the variance.
@return Total variance of vectors.
*/
double variance(std::vector<VectorXd> v) {

    double mean = average(v);
    double sum = 0;
    double n = 0;

    for (size_t i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            sum += (v[i][j] - mean)*(v[i][j] - mean);
            n++;
        }
    }

    //not this is n, not n-1
    return sum/(n);
}

/**
Finds the variance of the values in a vector.

@param v Vector which is used to calculate the variance.
@return Variance of vector v.
*/
double variance(VectorXd &v) {

    double mean = v.mean();
    double var = (v.array() - mean).pow(2).sum();

    return var/(v.rows()-1);
}


VectorXd concatenate(std::vector<VectorXd> &v) {

    VectorXd cat = v[0];

    for (size_t i = 1; i < v.size(); i++) {
        VectorXd temp(cat.rows() + v[i].rows());
        temp << cat, v[i];
        cat = temp;
    }

    return cat;
}

MatrixXd concatenate(std::vector<MatrixXd> &m) {

    MatrixXd cat = m[0];

    for (size_t i = 1; i < m.size(); i++) {
        MatrixXd temp(cat.rows() + m[i].rows(), m[i].cols());
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

VectorXd replaceNAN(VectorXd & V, double value) {
    int nrow = V.rows();

    for (int i = 0; i < nrow; i++)
        if (std::isnan(V[i]))
            V[i] = value;

    return V;
}


double center(VectorXd &v) {

	double mean = v.mean();
	double var = (v.array() - mean).pow(2).sum();

	return var / (v.rows() - 1);
}



/**
Creates a subset of vector v by taking all indexes i such that where[i] == equals.
@param v The vector to subset.
@param where The vector to condition the subset on.
@param equals Specifies which indexes to retain.

@requires Vector v and vector where must have the same number of rows
@return A subset of v.
*/
VectorXd extractRows(VectorXd &v, VectorXi &where, int equals) {

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

/**
Creates a subset of matrix m (rows) by taking all indexes i such that where[i] == equals.
@param m The matrix to extract rows from.
@param where The vector to condition the subset on.
@param equals Specifies which indexes to retain.

@requires Matrix m and vector where must have the same number of rows.
@return A subset of the rows of m.
*/
MatrixXd extractRows(MatrixXd &m, VectorXi &where, int equals) {

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

/**
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V1, V2 or M
- 0 otherwise

@param V1 Vector to check for NANs.
@param V2 Vector to check for NANs.
@param M Matrix to check for NANs (rows).

@requires V1, V2, M have the same number of rows.
@return A vector specifing where NANs appear in V1, V2 or M.
*/
VectorXi whereNAN(VectorXd &V1, VectorXd &V2, MatrixXd &M) {

    int nrow = V1.rows();
    int ncol = M.cols();

    VectorXi isNAN(nrow);

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

/**
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V or M
- 0 otherwise

@param V Vector to check for NANs.
@param M Matrix to check for NANs (rows).

@requires V, M have the same number of rows.
@return A vector specifing where NANs appear in V or M.
*/
VectorXi whereNAN(VectorXd &V, MatrixXd &M) {
    int nrow = V.rows();
    int ncol = M.cols();

    VectorXi isNAN(nrow);

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

/**
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V1 or V2
- 0 otherwise

@param V1 Vector to check for NANs.
@param V2 Vector to check for NANs.

@requires V1, V2 have the same number of rows.
@return A vector specifing where NANs appear in V1 or V2.
*/
//todo: print warning if many rows are removed due to NA?
VectorXi whereNAN(VectorXd &V1, VectorXd &V2) {
    int nrow = V1.rows();
    VectorXi toRemove(nrow);

    for (int i = 0; i < nrow; i++) {
        toRemove[i] = 0;
        if (std::isnan(V1[i]) || std::isnan(V2[i]))
            toRemove[i] = 1;
    }
    return toRemove;
}

/**
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V
- 0 otherwise

@param V Vector to check for NANs.
@return A vector specifing where NANs appear in V.
*/
VectorXi whereNAN(VectorXd &V) {
    int nrow = V.rows();
    VectorXi toRemove(nrow);
    for (int i = 0; i < nrow; i++) {
        toRemove[i] = 0;

        if (std::isnan(V[i]))
            toRemove[i] = 1;
        else
            toRemove[i] = 0;
    }

    return toRemove;
}

