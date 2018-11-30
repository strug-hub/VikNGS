#pragma once
#include "../Math/Math.h"

MatrixXd getRobustVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                             std::map<int, Depth>& d, VectorXd robustVar, bool rvs);

MatrixXd getRobustVarianceNormal(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                           std::map<int, Depth>& d, VectorXd robustVar, bool rvs);

MatrixXd getRegularVariance(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& MU, Family family);

VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X);
