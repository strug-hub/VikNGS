#pragma once
#include <vector>
#include <random>
#include "../vikNGS.h"
#include "../Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::VectorXi;

//RandomHelper.cpp
int randomInt(int from, int to);
double randomDouble(double from, double to);
double randomNormal(double mean, double sd);
int randomBinomial(int trials, double success);
std::vector<VectorXd> groupwiseShuffleWithoutReplacement(std::vector<VectorXd> &v);
std::vector<MatrixXd> groupwiseShuffleWithReplacement(std::vector<MatrixXd> &m);
std::vector<VectorXd> groupwiseShuffleWithReplacement(std::vector<VectorXd> &v);
std::vector<MatrixXd> shuffleWithoutReplacement(std::vector<MatrixXd> &m);
std::vector<VectorXd> shuffleWithoutReplacement(std::vector<VectorXd> &v);

//VectorHelper.cpp
double average(std::vector<VectorXd> v);
double variance(VectorXd &v);
double variance(std::vector<VectorXd> v);
VectorXd concatenate(std::vector<VectorXd> &v);
MatrixXd concatenate(std::vector<MatrixXd> &m);
MatrixXd replaceNAN(MatrixXd & M, double value);
VectorXd replaceNAN(VectorXd & V, double value);
inline MatrixXd nanToZero(MatrixXd &M) { return replaceNAN(M, 0); }
inline VectorXd nanToZero(VectorXd &V) { return replaceNAN(V, 0); }
VectorXd extractRows(VectorXd &v, VectorXi &where, int equals);
MatrixXd extractRows(MatrixXd &m, VectorXi &where, int equals);
VectorXi whereNAN(VectorXd &V1, VectorXd &V2, MatrixXd &M);
VectorXi whereNAN(VectorXd &V, MatrixXd &M);
VectorXi whereNAN(VectorXd &V);
VectorXi whereNAN(VectorXd &V1, VectorXd &V2);

//GeneticsHelper.cpp
double calcRobustVar(Vector3d &P);
VectorXd calculateGenotypeCalls(std::vector<Vector3d>& gl, Vector3d &P);
VectorXd calculateExpectedGenotypes(std::vector<Vector3d> &likelihood, Vector3d &P);
Vector3d calculateGenotypeFrequencies(std::vector<Vector3d> &likelihood);
Vector3d calculateGenotypeFrequencies(VectorXd & gt);

//StatisticsHelper.cpp
int maxValue(double zero, double one, double two);
double pnorm(double x);
double chiSquareOneDOF(double);
double qfc(std::vector<double> lambda, double evalpoint, int ncoef);
MatrixXd covariance(MatrixXd &M);
MatrixXd correlation(MatrixXd &M);
MatrixXd calculateHatMatrix(MatrixXd &Z);
MatrixXd calculateHatMatrix(MatrixXd &Z, MatrixXd &W, MatrixXd &sqrtW);
VectorXd getBeta(VectorXd &X, VectorXd &Y, MatrixXd &Z, Family family);
VectorXd getBeta(VectorXd &Y, MatrixXd &Z, Family family);
std::vector<VectorXd> fitModel(VectorXd &beta, std::vector<VectorXd> &y, std::vector<MatrixXd> &z, Family family);

inline double sigmoid(VectorXd &x, VectorXd &beta){
    return (1.0 / (1.0 + exp(-(beta.dot(x)))));
}
VectorXd logisticRegression(VectorXd &Y, MatrixXd &Z);
inline VectorXd CovariateNormalRegression(VectorXd &Y, MatrixXd &Z) {
    return Z.householderQr().solve(Y);
}
inline VectorXd CovariateBinomialRegression(VectorXd &Y, MatrixXd &Z) {
    return logisticRegression(Y,Z);
}
