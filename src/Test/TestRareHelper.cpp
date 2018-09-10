#include "Test.h"
#include "TestObject.h"

VectorXd getScoreVector(VectorXd& Ycenter, MatrixXd& X) {
    int nsnp = X.cols();
    VectorXd score(nsnp);
    for(int i = 0; i < nsnp; i++)
        score[i] = getScore(Ycenter, X, i);
    return score;
}

MatrixXd getVarianceMatrix(TestObject& o, Test& test, Family family){
    if(test.isExpectedGenotypes()){

        if(family == Family::BINOMIAL)
            return getVarianceBinomial(*o.getYcenter(), *o.getX(), *o.getG(),
                                       *o.getDepths(), o.robustVarVector(), test.isRVS());
        if(family == Family::NORMAL)
            return getVarianceNormal(*o.getYcenter(), *o.getX(), *o.getG(),
                                     *o.getDepths(), o.robustVarVector(), test.isRVS());
    }
    else{

        if(o.hasCovariates())
            return getVarianceRegularCovariates(*o.getYcenter(), *o.getX(), *o.getZ(), *o.getMU(), family);
        else
            return getVarianceRegular(*o.getYcenter(), *o.getX(), *o.getMU(), family);

    }
    throwError("TEST", "Don't know how to compute variance during test, this error should not happen.");
}

MatrixXd getVarianceBinomial(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                             std::map<int, Depth>& d, VectorXd robustVar, bool rvs){
    int nsnp = X.cols();

    MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

    double ym_hrd = 0; double ym_lrd = 0;

    for (int i = 0; i < Ycenter.rows(); i++)
        if (d[G[i]] == Depth::HIGH)
            ym_hrd += Ycenter[i] * Ycenter[i];
        else
            ym_lrd += Ycenter[i] * Ycenter[i];

    MatrixXd diagYm_hrd = VectorXd::Constant(nsnp, sqrt(ym_hrd)).asDiagonal();
    MatrixXd diagYm_lrd = VectorXd::Constant(nsnp, sqrt(ym_lrd)).asDiagonal();
    MatrixXd diagRobustVar = robustVar.asDiagonal();

    std::vector<MatrixXd> x = splitIntoGroups(X, G);

    for (size_t i = 0; i < x.size(); i++) {

        if (d[i] == Depth::HIGH) {
            if (rvs) {
                MatrixXd cor = correlation(x[i]);
                MatrixXd var_hrd = diagRobustVar.transpose() * cor * diagRobustVar;
                diagS += diagYm_hrd * var_hrd * diagYm_hrd;
            }
            else
                diagS += diagYm_hrd * covariance(x[i]) * diagYm_hrd;
        }
        else
            diagS += diagYm_lrd * covariance(x[i]) * diagYm_lrd;
    }

    return diagS;
}


MatrixXd getVarianceNormal(VectorXd& Ycenter, MatrixXd& X, VectorXi& G,
                           std::map<int, Depth>& d, VectorXd robustVar, bool rvs){
    int nsnp = X.cols();

    MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

    MatrixXd diagRobustVar = robustVar.asDiagonal();
    double ym = Ycenter.array().pow(2).sum();

    MatrixXd diagYm = VectorXd::Constant(nsnp, sqrt(ym)).asDiagonal();
    MatrixXd diagVar = MatrixXd::Constant(nsnp, nsnp, 0);
    double n = 0;

    std::vector<MatrixXd> x = splitIntoGroups(X, G);

    for (size_t i = 0; i < x.size(); i++) {
        n += x[i].rows();
        if (d[i] == Depth::HIGH && rvs){
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

MatrixXd getVarianceRegular(VectorXd& Ycenter, MatrixXd& X, VectorXd& MU, Family family){

    VectorXd var1;

    if(family == Family::BINOMIAL)
        var1 = MU.array() * (1 - MU.array());
    else{
       double var = (Ycenter.array().pow(2).sum())/MU.rows();
       var1 = VectorXd::Constant(MU.rows(), var);
    }

    int nsnp = X.cols();
    MatrixXd sum1 = MatrixXd::Constant(nsnp, nsnp, 0);
    VectorXd sum2 = VectorXd::Constant(nsnp, 0);
    MatrixXd sum3 = MatrixXd::Constant(1, 1, 0);
    for(int i = 0; i < X.rows(); i++){

        MatrixXd one = var1[i] * (X.row(i).transpose() * X.row(i)).array();
        sum1 += one;

        VectorXd two = var1[i] * X.row(i).transpose().array();
        sum2 += two;

        sum3(0,0) += var1[i];
    }

    MatrixXd var = sum2 * sum3.inverse() * sum2.transpose();
    return sum1 - var;
}

MatrixXd getVarianceRegularCovariates(VectorXd& Ycenter, MatrixXd& X, MatrixXd& Z, VectorXd& MU, Family family){

    VectorXd var1;

    if(family == Family::BINOMIAL)
        var1 = MU.array() * (1 - MU.array());
    else{
       double var = (Ycenter.array().pow(2).sum())/MU.rows();
       var1 = VectorXd::Constant(MU.rows(), var);
    }

    int nsnp = X.cols();
    MatrixXd sum1 = MatrixXd::Constant(nsnp, nsnp, 0);
    int ncov = Z.cols();
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
