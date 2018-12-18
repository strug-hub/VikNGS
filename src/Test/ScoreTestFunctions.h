#pragma once
#include "../Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::VectorXi;

class Group;
#include "../Enums.h"

VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X);

MatrixXd getRobustVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs);
MatrixXd getRobustVarianceNormal(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs);
MatrixXd getRegularVariance(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& MU, Family family);
