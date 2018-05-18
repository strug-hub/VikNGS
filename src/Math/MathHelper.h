#pragma once
#include <vector>
#include <random>
#include "../Log.h"
#include "../Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;


VectorXd concatenate(std::vector<VectorXd> &v);
MatrixXd concatenate(std::vector<MatrixXd> &m);

/*
Creates a subset of vector v by taking all indexes i such that where[i] == equals.
@param v The vector to subset.
@param where The vector to condition the subset on.
@param equals Specifies which indexes to retain.

@requires Vector v and vector where must have the same number of rows
@return A subset of v.
*/
VectorXd extractRows(VectorXd &v, VectorXd &where, double equals);

MatrixXd replaceNAN(MatrixXd & M, double value);

/*
Creates a subset of matrix m (rows) by taking all indexes i such that where[i] == equals.
@param m The matrix to extract rows from.
@param where The vector to condition the subset on.
@param equals Specifies which indexes to retain.

@requires Matrix m and vector where must have the same number of rows.
@return A subset of the rows of m.
*/
MatrixXd extractRows(MatrixXd &m, VectorXd &where, double equals);

/*
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V1, V2 or M
- 0 otherwise

@param V1 Vector to check for NANs.
@param V2 Vector to check for NANs.
@param M Matrix to check for NANs (rows).

@requires V1, V2, M have the same number of rows.
@return A vector specifing where NANs appear in V1, V2 or M.
*/
VectorXd whereNAN(VectorXd &V1, VectorXd &V2, MatrixXd &M);

/*
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V or M
- 0 otherwise

@param V Vector to check for NANs.
@param M Matrix to check for NANs (rows).

@requires V, M have the same number of rows.
@return A vector specifing where NANs appear in V or M.
*/
VectorXd whereNAN(VectorXd &V, MatrixXd &M);

/*
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V
- 0 otherwise

@param V Vector to check for NANs.
@return A vector specifing where NANs appear in V.
*/
VectorXd whereNAN(VectorXd &V);

/*
Creates a vector such that:
- the vector contains a 1 at if a NAN exists in a given row of V1 or V2
- 0 otherwise

@param V1 Vector to check for NANs.
@param V2 Vector to check for NANs.

@requires V1, V2 have the same number of rows.
@return A vector specifing where NANs appear in V1 or V2.
*/
VectorXd whereNAN(VectorXd &V1, VectorXd &V2);

/*
Finds the mean of the values in a vector of vectors.

@param v Vector of vectors which is used to calculate the mean.
@return The mean value of all the vectors.
*/
double average(std::vector<VectorXd> v);

/*
Finds the variance of the values in a vector.

@param v Vector which is used to calculate the variance.
@return Variance of vector v.
*/
double variance(VectorXd &v);

/*
Finds the variance of the values in a vector of vectors.
Note: uses n, not n-1 for final division

@param v Vector of vectors which is used to calculate the variance.
@return Total variance of vectors.
*/
double variance(std::vector<VectorXd> v);



//StatisticsHelper.cpp
VectorXd getBeta(VectorXd &X, VectorXd &Y, MatrixXd &Z, std::string family);
VectorXd getBeta(VectorXd &Y, MatrixXd &Z, std::string family);
std::vector<VectorXd> fitModel(VectorXd &beta, std::vector<VectorXd> &y, std::vector<MatrixXd> &z, std::string distribution);
VectorXd shuffleWithoutReplacement(std::vector<VectorXd> &v);
MatrixXd shuffleWithReplacement(std::vector<MatrixXd> &m);

/*
Calculates the robust variance of E(G | D). var(x) = E(x^2) - E(x)^2

@param p Genotype frequency from EM algorithm
@return Robust variance of E(G | D)
*/
double calcRobustVar(VectorXd &p);

/*
Finds p-value for test statistic using a chi-squared distribution with one degree of freedom
using chi-squared probability density function.

@param statistic Test statistic.
@return p-value.
*/
double chiSquareOneDOF(double);

std::vector<double> randomSample(std::vector<double> &, int);
VectorXd CovariateRegression(VectorXd &Y, MatrixXd &Z);
int randomInt(int from, int to);
double randomDouble(double from, double to);
double randomNormal(double mean, double sd);
MatrixXd covariance(MatrixXd &M);
MatrixXd correlation(MatrixXd &M);
VectorXd logisticRegression(VectorXd &Y, MatrixXd &Z);

MatrixXd calculateHatMatrix(MatrixXd &Z);
MatrixXd calculateHatMatrix(MatrixXd &Z, MatrixXd &W, MatrixXd &sqrtW);

/*
Approximates the p-value from the pdf of the normal distribution where x is a Z-score

@param x Z-score.
@return p-value.
*/
double pnorm(double x);

double qfc(std::vector<double> lambda, double evalpoint, int ncoef);


