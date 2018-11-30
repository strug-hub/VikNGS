#include "ScoreTestFunctions.h"
#include "TestObject.h"
#include "Group.h"

VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X) {
    int nsnp = X.cols();
    VectorXd score(nsnp);
    for(int i = 0; i < nsnp; i++)
        score[i] = (Ycenter.array() * X.col(i).array()).sum();
    return score;
}

MatrixXd getVarianceMatrix2(VectorXd& Ycenter, VectorXd& Mu, MatrixXd& X, MatrixXd& Z, MatrixXd& P, Test& test, Family family){

    if(test.getVariance() == Variance::RVS){

  //      if(family == Family::BINOMIAL)
    //        return getRobustVarianceBinomial(Ycenter, X, *o.getG(), *o.getDepths(), o.robustVarVector(), test.isRVS());
   //     if(family == Family::NORMAL)
   //         return getRobustVarianceNormal(Ycenter, X, *o.getG(),
    //                                 *o.getDepths(), o.robustVarVector(), test.isRVS());
    }
    if(test.getVariance() == Variance::REGULAR)
        return getRegularVariance(Ycenter, X, Z, Mu, family);

    throwError("ScoreTestFunctions", "Unsure how to calculate variance in score test. This should not happen.");
}

MatrixXd getVarianceMatrix(TestObject& o, Test& test, Family family){
    if(test.isExpectedGenotypes()){

        if(family == Family::BINOMIAL)
            return getRobustVarianceBinomial(*o.getYcenter(), *o.getX(), *o.getGroup(), o.robustVarVector(), test.isRVS());
        if(family == Family::NORMAL)
            return getRobustVarianceNormal(*o.getYcenter(), *o.getX(), *o.getGroup(), o.robustVarVector(), test.isRVS());
    }
    else
            return getRegularVariance(*o.getYcenter(), *o.getX(), *o.getZ(), *o.getMU(), family);


    throwError("TEST", "Don't know how to compute variance during test, this error should not happen.");
}

MatrixXd getRobustVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs){
    int nsnp = X.cols();

    MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);
    MatrixXd diagRobustVar = robustVar.asDiagonal();

    std::vector<MatrixXd> x = splitIntoGroups(X, group);
    std::vector<VectorXd> y = splitIntoGroups(Ycenter, group);
    double minus1Factor = X.rows()*1.0/(X.rows()-1);

    for (size_t i = 0; i < x.size(); i++) {

        double ym = y[i].array().pow(2).sum();
        ym = ym * minus1Factor;

        MatrixXd diagYm = VectorXd::Constant(nsnp, sqrt(ym)).asDiagonal();

        if (group.depth(i) == Depth::HIGH && rvs) {
            MatrixXd var = diagRobustVar.transpose() * correlation(x[i]) * diagRobustVar;
            diagS += diagYm * var * diagYm;
        }
        else
            diagS += diagYm * covariance(x[i]) * diagYm;
    }

    return diagS;
}

MatrixXd getRobustVarianceNormal(VectorXd& Ycenter, MatrixXd& X, Group& group, VectorXd robustVar, bool rvs){
    int nsnp = X.cols();

    MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);
    MatrixXd diagRobustVar = robustVar.asDiagonal();
    double ym = Ycenter.array().pow(2).sum();

    MatrixXd diagYm = VectorXd::Constant(nsnp, sqrt(ym)).asDiagonal();
    MatrixXd diagVar = MatrixXd::Constant(nsnp, nsnp, 0);
    double n = 0;

    std::vector<MatrixXd> x = splitIntoGroups(X, group);

    for (size_t i = 0; i < x.size(); i++) {
        n += x[i].rows();
        if (group.depth(i) == Depth::HIGH && rvs){
            MatrixXd v = x[i].rows() * (diagRobustVar.transpose() * correlation(x[i]) * diagRobustVar).array();
            diagVar = diagVar + v;
        }
        else{
            MatrixXd v = x[i].rows() * covariance(x[i]).array();
            diagVar = diagVar + v;
        }
    }

    diagS = 1.0/n * (diagYm * diagVar * diagYm).array();

    return diagS;
}

MatrixXd getRegularVariance(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& Mu, Family family){

    VectorXd var1;

    if(family == Family::BINOMIAL)
        var1 = Mu.array() * (1 - Mu.array());
    else{
       double var = (Ycenter.array().pow(2).sum())/Mu.rows();
       var1 = VectorXd::Constant(Mu.rows(), var);
    }

    int nsnp = X.cols();
    int ncov = Z.cols();

    MatrixXd sum1 = MatrixXd::Constant(nsnp, nsnp, 0);
    MatrixXd sum2 = MatrixXd::Constant(nsnp, ncov, 0);
    MatrixXd sum3 = MatrixXd::Constant(ncov, ncov, 0);

    for(int i = 0; i < X.rows(); i++){
        MatrixXd one = var1[i] * (X.row(i).transpose() * X.row(i)).array();
        sum1 += one;

        MatrixXd two = var1[i] * (X.row(i).transpose() * Z.row(i)).array();
        sum2 += two;

        MatrixXd three = var1[i] * (Z.row(i).transpose() * Z.row(i)).array();
        sum3 += three;
    }

    MatrixXd var = sum2 * sum3.inverse() * sum2.transpose();
    return sum1 - var;
}

