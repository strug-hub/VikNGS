#include "Math.h"

static const std::string ERROR_SOURCE = "VECTOR_HELPER";

double variance(VectorXd& X, VectorXi& G, int group){

    double mean = 0;
    double sum = 0;
    double n = 0;

    for (int i = 0; i < X.rows(); i++) {
        if(G[i] == group){
            sum += X[i];
            n++;
        }
    }

    mean = sum / n;
    sum = 0;

    for (int i = 0; i < X.rows(); i++)
        if(G[i] == group)
            sum += (X[i] - mean)*(X[i] - mean);

    //this is n, not n-1
    return sum/(n-1);
}

/**
Determines if there is variation in the values in a vector.

@param v Vector
@return true if variance is present.
*/
bool hasVariance(VectorXd &v) {

    int j;

    double mean = 0;
    double n = 0;
    for (j = 0; j < v.rows(); j++) {
        if (!std::isnan(v[j])) {
            mean += v[j];
            n++;
        }
    }
    mean = mean / n;

    double sd = 0;
    for (j = 0; j < v.size(); j++)
        if (!std::isnan(v[j]))
            sd += pow((v[j] - mean), 2);

    if (1e-8*n < sd)
        return true;

    return false;
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

/**
Finds the variance of the values in a column of a matrix.

@param M Matrix which is used to calculate the variance.
@param column Column index of M where the variance will be calculated.
@return Variance of column in M.
*/
double variance(MatrixXd& M, int column) {

    double mean = M.col(column).mean();
    double var = (M.col(column).array() - mean).pow(2).sum();

    return var/(M.rows()-1);
}

MatrixXd subtractGroupMean(MatrixXd& M, VectorXi& G){

    MatrixXd Mcenter(M.rows(), M.cols());
    int groupID;

    for(int j = 0; j < M.cols(); j++){
        std::map<int, double> mean;
        std::map<int, double> n;

        for (int i = 0; i < M.rows(); i++){
            groupID = G[i];
            if(mean.count(groupID) < 1){
                mean[groupID] = M(i,j);
                n[groupID] = 1;
            }
            else{
                mean[groupID] += M(i,j);
                n[groupID]++;
            }
        }

        for (std::pair<int, double> e : mean)
            mean[e.first] = e.second / n[e.first];

        for (int i = 0; i < M.rows(); i++)
            Mcenter(i,j) = M(i,j) - mean[G[i]];
    }

    return Mcenter;
}

std::vector<VectorXd> splitIntoGroups(VectorXd& v, VectorXi& g, int gSize){
    int ngroups;
    if(gSize > 0)
        ngroups = gSize;
    else
        ngroups = 1 + g.maxCoeff();

    std::vector<VectorXd> result(static_cast<size_t>(ngroups));

    for (int i = 0; i < ngroups; i++)
        result[i] = extractRows(v, g, i);

    return result;
}

std::vector<MatrixXd> splitIntoGroups(MatrixXd& m, VectorXi& g, int gSize){
    int ngroups;
    if(gSize > 0)
        ngroups = gSize;
    else
        ngroups = 1 + g.maxCoeff();

    std::vector<MatrixXd> result(static_cast<size_t>(ngroups));

    for (int i = 0; i < ngroups; i++)
        result[i] = extractRows(m, g, i);

    return result;
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
Creates a subset of vector v by taking all indexes i such that where[i] == equals.
@param v The vector to subset.
@param where The vector to condition the subset on.
@param equals Specifies which indexes to retain.

@requires Vector v and vector where must have the same number of rows
@return A subset of v.
*/
VectorXi extractRows(VectorXi &v, VectorXi &where, int equals) {

    VectorXi subset(v.rows());
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
- the vector contains a 1 at if a NAN exists in a given row of M
- 0 otherwise

@param M Matrix to check for NANs (rows).
@return A vector specifing where NANs appear in M
*/
VectorXi whereNAN(MatrixXd &M) {

    int nrow = M.rows();
    int ncol = M.cols();

    VectorXi isNAN(nrow);

    for (int i = 0; i < nrow; i++) {
        isNAN[i] = 0;

        for (int j = 0; j < ncol; j++)
            if (std::isnan(M(i, j)))
                isNAN[i] = 1;
    }

    return isNAN;
}

/**
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V
- 0 otherwise

@param V Vector to check for NANs.
@return A vector specifing where NANs appear in V
*/
VectorXi whereNAN(VectorXd &V) {

    int nrow = V.rows();

    VectorXi isNAN(nrow);

    for (int i = 0; i < nrow; i++)
        if (std::isnan(V[i]))
            isNAN[i] = 1;
        else
            isNAN[i] = 0;

    return isNAN;
}
