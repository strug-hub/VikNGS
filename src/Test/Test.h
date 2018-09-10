#pragma once
#include "../Math/Math.h"
#include "../vikNGS.h"
#include "TestObject.h"

inline double getScore(VectorXd& Ycenter, MatrixXd& X, int col=0) {
    VectorXd x = X.col(col);
    return (Ycenter.array() * X.col(col).array()).sum();
}

//Test.cpp
double runTest(SampleInfo* sampleInfo, VariantSet* variant, Test test, int nboot, bool stopEarly);

//TestCommonHelper.cpp
double getVariance(TestObject& o, Test& test, Family family);
double getVarianceBinomial(VectorXd& Ycenter, VectorXd& X, VectorXi& G,
                           std::map<int, Depth>& d, double robustVar, bool rvs);
double getVarianceNormal(VectorXd& Y, VectorXd& X, VectorXi& G,
                          std::map<int, Depth>& d, double robustVar, bool rvs);
double getVarianceNormalCovariates(VectorXd& Ycenter, VectorXd& X, VectorXi& G,
                          std::map<int, Depth>& d, double robustVar, bool rvs);
double getVarianceRegular(VectorXd& Y, VectorXd& X, VectorXd& MU,  Family family);
double getVarianceRegularCovariates(VectorXd& Ycenter, VectorXd& X, MatrixXd& Z, VectorXd& MU, Family family);

double getVarianceMatrix(VectorXd* Y, MatrixXd* X, int col, VectorXd* MU,  Family family);


//TestRareHelper.cpp
VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X);
MatrixXd getVarianceMatrix(TestObject& o, Test& test, Family family);
MatrixXd getVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                             std::map<int, Depth>& d, VectorXd robustVar, bool rvs);
MatrixXd getVarianceNormal(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                           std::map<int, Depth>& d, VectorXd robustVar, bool rvs);
MatrixXd getVarianceRegular(VectorXd& Ycenter, MatrixXd& X, VectorXd& MU, Family family);
MatrixXd getVarianceRegularCovariates(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& MU, Family family);


