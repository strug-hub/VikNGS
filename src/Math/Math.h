#pragma once
#include "../Enums.h"
#include "../Log.h"
#include <vector>
#include <map>
#include <random>
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
MatrixXd groupwiseShuffleWithReplacement(MatrixXd& M, VectorXi& G, std::map<int, std::vector<int>>& group);
VectorXd groupwiseShuffleWithReplacement(VectorXd& V, VectorXi& G, std::map<int, std::vector<int>>& group);
VectorXd groupwiseShuffleWithoutReplacement(VectorXd& V, VectorXi& G, std::map<int, std::vector<int>>& group);
MatrixXd shuffleColumnwiseWithoutReplacement(MatrixXd &M);

//VectorHelper.cpp
bool hasVariance(VectorXd &v);
double variance(VectorXd &v);
double variance(MatrixXd& M, int column);
double variance(VectorXd& X, VectorXi& G, int group);
MatrixXd subtractGroupMean(MatrixXd& M, VectorXi& G);
std::vector<VectorXd> splitIntoGroups(VectorXd& v, VectorXi& g, int gSize=-1);
std::vector<MatrixXd> splitIntoGroups(MatrixXd& m, VectorXi& g, int gSize=-1);
MatrixXd replaceNAN(MatrixXd& M, double value);
VectorXd replaceNAN(VectorXd& V, double value);
inline MatrixXd nanToZero(MatrixXd &M) { return replaceNAN(M, 0); }
inline VectorXd nanToZero(VectorXd &V) { return replaceNAN(V, 0); }
VectorXd extractRows(VectorXd& v, VectorXi& where, int equals);
VectorXi extractRows(VectorXi& v, VectorXi& where, int equals);
MatrixXd extractRows(MatrixXd& m, VectorXi& where, int equals);
VectorXi whereNAN(VectorXd& V);
VectorXi whereNAN(MatrixXd& M);

//GeneticsHelper.cpp
double calcRobustVar(Vector3d &P);
double calcRobustVar(double P1, double P2);
VectorXd calculateGenotypeCalls(std::vector<Vector3d>& gl, Vector3d &P);
VectorXd calculateExpectedGenotypes(std::vector<Vector3d> &likelihood, Vector3d &P);
Vector3d calculateGenotypeFrequencies(std::vector<Vector3d> &likelihood);
Vector3d calculateGenotypeFrequencies(VectorXd & gt);

//StatisticsHelper.cpp
int maxValue(double zero, double one, double two);
double pnorm(double x);
double chiSquareOneDOF(double);
MatrixXd covariance(MatrixXd &M);
MatrixXd correlation(MatrixXd &M);
MatrixXd calculateHatMatrix(MatrixXd &Z);
MatrixXd calculateHatMatrix(MatrixXd &Z, MatrixXd &W, MatrixXd &sqrtW);
VectorXd getBeta(VectorXd &X, VectorXd &Y, MatrixXd &Z, Family family);
VectorXd getBeta(VectorXd &Y, MatrixXd &Z, Family family);
VectorXd fitModel(VectorXd &beta, MatrixXd &Z, Family family);

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
