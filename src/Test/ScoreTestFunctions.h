#pragma once
#include "../Math/EigenStructures.h"

class Group;
enum class Family;

VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X);

MatrixXd getRobustVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs);
MatrixXd getRobustVarianceNormal(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs);
MatrixXd getRegularVariance(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& MU, Family family);
