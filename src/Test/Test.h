#pragma once
#include "../Math/Math.h"
#include "../vikNGS.h"
#include "TestObject.h"

inline double getScore(VectorXd& Ycenter, MatrixXd& X, int col=0) {
    VectorXd y = Ycenter;
    MatrixXd x = X;

    double yy = y.sum();
    double xx = x.sum()/x.rows();

    double v = 0;
    double v2 = (Ycenter.array() * X.col(col).array()).sum();

    for(int i = 0; i < Ycenter.size(); i++){

        v += Ycenter[i] * X(i,col);

    }


    return (Ycenter.array() * X.col(col).array()).sum();
}

//Test.cpp
double runTest(SampleInfo* sampleInfo, VariantSet* variant, Test test, int nboot, bool stopEarly);

//TestCommonHelper.cpp
double getVariance(TestObject& o, Test& test, Family family);
double getVarianceBinomial(VectorXd& Ycenter, VectorXd& X, VectorXi& G,
                           std::map<int, Depth>& d, double robustVar, bool rvs);
double getVarianceNormal(VectorXd& Ycenter, VectorXd& X, VectorXi& G,
                          std::map<int, Depth>& d, double robustVar, bool rvs);
double getVarianceRegular(VectorXd& Y, VectorXd& X, VectorXd& MU,  Family family);
double getVarianceRegularCovariates(VectorXd& Ycenter, VectorXd& X, MatrixXd& Z, VectorXd& MU, Family family);

//TestRareHelper.cpp
VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X);
MatrixXd getVarianceMatrix(TestObject& o, Test& test, Family family);
MatrixXd getVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                             std::map<int, Depth>& d, VectorXd robustVar, bool rvs);
MatrixXd getVarianceNormal(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                           std::map<int, Depth>& d, VectorXd robustVar, bool rvs);
MatrixXd getVarianceRegular(VectorXd& Ycenter, MatrixXd& X, VectorXd& MU, Family family);
MatrixXd getVarianceRegularCovariates(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& MU, Family family);


