#pragma once
#include "../Math/Math.h"

class Group;

VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X);

MatrixXd getRobustVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs);
MatrixXd getRobustVarianceNormal(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs);
MatrixXd getRegularVariance(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& MU, Family family);
