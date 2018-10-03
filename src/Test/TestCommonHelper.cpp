#include "Test.h"
#include "TestObject.h"


double getVariance(TestObject& o, Test& test, Family family){

    VectorXd X = o.getX(0);
    if(test.isExpectedGenotypes()){

        if(o.hasCovariates()){
            //todo: getVarianceBinomialCovariates
            if(family == Family::BINOMIAL)
                return getVarianceBinomial(*o.getYcenter(), X, *o.getG(),
                                           *o.getDepths(), o.getRobustVar(), test.isRVS());
            if(family == Family::NORMAL)
                return getVarianceNormal(*o.getYcenter(), X, *o.getG(),
                                                   *o.getDepths(), o.getRobustVar(), test.isRVS());
        }
        else{

            if(family == Family::BINOMIAL)
                return getVarianceBinomial(*o.getYcenter(), X, *o.getG(),
                                           *o.getDepths(), o.getRobustVar(), test.isRVS());
            if(family == Family::NORMAL)
                return getVarianceNormal(*o.getYcenter(), X, *o.getG(),
                                         *o.getDepths(), o.getRobustVar(), test.isRVS());
        }
    }
    else{

        if(o.hasCovariates())
            return getVarianceRegularCovariates(*o.getYcenter(), X, *o.getZ(),
                                                    *o.getMU(), family);
        else
            return getVarianceRegular(*o.getY(), X, *o.getMU(), family);

    }
    throwError("TEST", "Don't know how to compute variance during test, this error should not happen.");
}

double getVarianceBinomial(VectorXd& Ycenter, VectorXd& X, VectorXi& G,
                           std::map<int, Depth>& d, double robustVar, bool rvs) {
    double var = 0;
    std::map<int, double> var_y;
    int gmax = 0;

    for (int i = 0; i < Ycenter.rows(); i++){
        if(var_y.count(G[i]) > 0)
            var_y[G[i]] += Ycenter[i] * Ycenter[i];
        else{
            var_y[G[i]] = Ycenter[i] * Ycenter[i];
            if(G[i] > gmax)
                gmax = G[i];
        }
    }

    for(int i = 0; i <= gmax; i++){

        if(var_y.count(i) < 1)
            continue;

        if (rvs && d[i] == Depth::HIGH)
            var += var_y[i] * robustVar;
        else
            var += var_y[i] * variance(X, G, i);
    }

    return var;
}

double normalVariance(double var_y, VectorXd& X, VectorXi& G,
                          std::map<int, Depth>& d, double robustVar, bool rvs) {
    std::map<int, double> n;

    double var = 0;
    int gmax = 0;

    for (int i = 0; i < G.rows(); i++)
        if(n.count(G[i]) > 0)
            n[G[i]]++;
        else{
            n[G[i]] = 1;
            if(G[i] > gmax)
                gmax = G[i];
        }

    for(int i = 0; i <= gmax; i++){

        if (rvs && d[i] == Depth::HIGH)
            var += n[i] * var_y * robustVar;
        else
            var += n[i] * var_y * variance(X, G, i);
    }

    return var / X.rows();
}

double getVarianceNormal(VectorXd& Ycenter, VectorXd& X, VectorXi& G,
                          std::map<int, Depth>& d, double robustVar, bool rvs) {
        double var_y = Ycenter.array().pow(2).sum();
        return normalVariance(var_y, X, G, d, robustVar, rvs);
}

double getVarianceRegular(VectorXd& Ycenter, VectorXd& X, VectorXd& MU, Family family) {

    VectorXd var1;

    if(family == Family::BINOMIAL)
        var1 = MU.array() * (1 - MU.array());
    else{
       double var = (Ycenter.array().pow(2).sum())/MU.rows();
       var1 = VectorXd::Constant(MU.rows(), var);
    }

    double sum1 = 0;

        double sum2 = 0;
        double sum3 = 0;
        for(int i = 0; i < X.rows(); i++){
            sum1 += var1[i] * X[i] * X[i];
            sum2 += var1[i] * X[i];
            sum3 += var1[i];
        }

    double var = (sum2 * (1.0/sum3) * sum2);
    return sum1 - var;
}


double getVarianceRegularCovariates(VectorXd& Ycenter, VectorXd& X, MatrixXd& Z, VectorXd& MU, Family family) {

    VectorXd var1;

    if(family == Family::BINOMIAL)
        var1 = MU.array() * (1 - MU.array());
    else{
       double var = (Ycenter.array().pow(2).sum())/MU.rows();
       var1 = VectorXd::Constant(MU.rows(), var);
    }

    double sum1 = 0;

    int ncov = Z.cols();
    VectorXd sum2 = VectorXd::Constant(ncov, 0);
    MatrixXd sum3 = MatrixXd::Constant(ncov, ncov, 0);

    for(int i = 0; i < X.rows(); i++){

        sum1 += var1[i] * X[i] * X[i];

        VectorXd two = var1[i] * X[i] * Z.row(i).array();
        sum2 += two;

        MatrixXd three = var1[i] * (Z.row(i).transpose() * Z.row(i)).array();
        sum3 += three;
    }

    double var = (sum2.transpose() * sum3.inverse() * sum2)(0);
    return sum1 - var;
}

//todo
/*
double  covariateCovarianceCalculation(MatrixXd* H, std::vector<MatrixXd>* x, int col){

    double cov1 = 0;
    for(int j = 0; j < H->cols() - 1; j++)
        for(int k = j+1; k < H->cols(); k++)
            cov1 += (H->col(j).array() * H->col(k).array() * Yvar.array()).sum();

    double cov2 = (*H*Yvar).sum() - H->diagonal().dot(Yvar);

    double meanx = average(*x, col);
    return meanx * meanx * (cov1 - cov2);
}

double covariateCovarianceCalculation(){

    double cov1 = 0;
    for(int j = 0; j < H.cols() - 1; j++)
        for(int k = j+1; k < H.cols(); k++)
            cov1 += (H.col(j).array() * H.col(k).array() * Yvar.array()).sum();

    double cov2 = (H*Yvar).sum() - H.diagonal().dot(Yvar);

    double meanx = average(this->x);
    return meanx * meanx * (cov1 - cov2);
}


double getVarianceBinomialCovariates(std::vector<VectorXd>& ycenter, std::vector<VectorXd>& x,
                           std::vector<Depth> readDepth, double robustVar, bool rvs) {

    if(useHatMatrix && covariates){
        double var = covariateVarianceCalculation(rvs);
        double cov = covariateCovarianceCalculation();
        return var + 2*cov;
    }
}
*/
