#include "RareTestCollapseObject.h"


void RareTestCollapseObject::setFamily(std::string fam){
    if(fam == "binomial"){
        this->normal = false;
        this->binomial = true;
        return;
    }
    if(fam == "normal"){
        this->normal = true;
        this->binomial = false;
        return;
    }
}

MatrixXd RareTestCollapseObject::getRobustVar(){
    VectorXd robustVar(size());
    for (int i = 0; i < size(); i++)
        robustVar[i] = sqrt(t[i].getRobustVar());

    return robustVar.asDiagonal();
}


MatrixXd RareTestCollapseObject::getYmHigh(){
    VectorXd ym(size());
    for (int i = 0; i < size(); i++)
        ym[i] = sqrt(t[i].getYmHigh(ycenter));

    return ym.asDiagonal();

}
MatrixXd RareTestCollapseObject::getYmLow(){
    VectorXd ym(size());
    for (int i = 0; i < size(); i++)
        ym[i] = sqrt(t[i].getYmLow(ycenter));

    return ym.asDiagonal();


}
MatrixXd RareTestCollapseObject::getYm(){
    VectorXd ym(size());
    for (int i = 0; i < size(); i++)
        ym[i] = sqrt(t[i].getYm(ycenter));

    return ym.asDiagonal();
}

VectorXd RareTestCollapseObject::getScore(){
    VectorXd score(size());
    for (int i = 0; i < size(); i++)
        score[i] = t[i].getScore(ycenter);
    return score;
}

MatrixXd RareTestCollapseObject::getVarianceRegular(){

    std::vector<MatrixXd> x = getX();
    VectorXd MU = concatenate(mu);

    VectorXd var1;

    if(binomial){
        var1 = MU.array() * (1 - MU.array());
    }
    else{
       VectorXd Y = concatenate(y);
       double var = ((Y - MU).array().pow(2).sum())/MU.rows();
       var1 = VectorXd::Constant(MU.rows(), var);
    }

    int nsnp = x[0].cols();
    MatrixXd sum1 = MatrixXd::Constant(nsnp, nsnp, 0);
    int index = 0;

    if(covariates){

        int ncov = z[0].cols();
        MatrixXd sum2 = MatrixXd::Constant(nsnp, ncov, 0);
        MatrixXd sum3 = MatrixXd::Constant(ncov, ncov, 0);

        for(int i = 0; i < x.size(); i++){
            for(int j = 0; j < x[i].rows(); j++){

                MatrixXd one = var1[index] * (x[i].row(j).transpose() * x[i].row(j)).array();
                sum1 += one;

                MatrixXd two = var1[index] * (x[i].row(j).transpose() * z[i].row(j)).array();
                sum2 += two;

                MatrixXd three = var1[index] * (z[i].row(j).transpose() * z[i].row(j)).array();
                sum3 += three;

                index++;
            }
        }

        MatrixXd var = sum2 * sum3.inverse() * sum2.transpose();
        return sum1 - var;
    }
    else{

        VectorXd sum2 = VectorXd::Constant(nsnp, 0);
        MatrixXd sum3 = MatrixXd::Constant(1, 1, 0);
        for(int i = 0; i < x.size(); i++){
            for(int j = 0; j < x[i].rows(); j++){


                MatrixXd one = var1[index] * (x[i].row(j).transpose() * x[i].row(j)).array();
                sum1 += one;

                VectorXd two = var1[index] * x[i].row(j).transpose().array();
                sum2 += two;

                sum3(0,0) += var1[index];

                index++;
            }
        }

        MatrixXd var = sum2 * sum3.inverse() * sum2.transpose();
        return sum1 - var;
    }
}


MatrixXd RareTestCollapseObject::getVarianceBinomial(bool rvs){
    int i, j;
    int nsnp = size();

    MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

    MatrixXd diagRobustVar = getRobustVar();
    MatrixXd diagYm_hrd = getYmHigh();
    MatrixXd diagYm_lrd = getYmLow();

    for (i = 0; i < t[0].size(); i++) {
        VectorXd x_0 = t[0].getX(i);
        MatrixXd x(x_0.rows(), nsnp);
        x.col(0) = x_0;

        for (j = 1; j < nsnp; j++)
            x.col(j) = t[j].getX(i);

        if (t[0].isHRG(i)) {
            if (rvs) {
                MatrixXd var_hrd = diagRobustVar.transpose() * correlation(x) * diagRobustVar;
                diagS += diagYm_hrd * var_hrd * diagYm_hrd;
            }
            else
                diagS += diagYm_hrd * covariance(x) * diagYm_hrd;
        }
        else
            diagS += diagYm_lrd * covariance(x) * diagYm_lrd;
    }

    return diagS;
}

MatrixXd RareTestCollapseObject::getVarianceNormal(bool rvs){
    int i, j;
    int nsnp = t.size();

    MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

    MatrixXd diagRobustVar = getRobustVar();
    MatrixXd diagYm = getYm();

    MatrixXd diagVar = MatrixXd::Constant(diagRobustVar.rows(), diagRobustVar.cols(), 0);
    double n = 0;

    for (i = 0; i < t[0].size(); i++) {
        VectorXd x_0 = t[0].getX(i);
        int groupSize = x_0.rows();
        n += groupSize;
        MatrixXd x(groupSize, nsnp);
        x.col(0) = x_0;

        for (j = 1; j < nsnp; j++)
            x.col(j) = t[j].getX(i);

        if (t[0].isHRG(i) && rvs){
            MatrixXd v = groupSize * (diagRobustVar.transpose() * correlation(x) * diagRobustVar).array();
            diagVar = diagVar + v;
        }
        else{
            MatrixXd v = groupSize * covariance(x).array();
            diagVar = diagVar + v;
        }
    }

    diagS = 1.0/n * (diagYm * diagVar * diagYm).array();

    return diagS;

}

void RareTestCollapseObject::bootstrap() {

    if(regular)
        regularBootstrap();
    else if(normal)
        normalBootstrap();
    else if (binomial)
        binomialBootstrap();

    return;
}

void RareTestCollapseObject::regularBootstrap() {

    for(int i = 0; i< size(); i++)
        t[i].regularBootstrap();

    if(covariates){
        z = shuffleWithoutReplacement(z_original);
        MatrixXd shuffleZ = concatenate(z);

        VectorXd beta = getBeta(Y, shuffleZ, "binomial");
        ycenter = fitModel(beta, y, z, "binomial");
    }
}

void RareTestCollapseObject::binomialBootstrap() {

    for(int i = 0; i< size(); i++)
        t[i].bootstrap();

    if(covariates){
        z = groupwiseShuffleWithReplacement(z_original);
        MatrixXd shuffleZ = concatenate(z);

        VectorXd beta = getBeta(Y, shuffleZ, "binomial");
        ycenter = fitModel(beta, y, z, "binomial");
    }
}


void RareTestCollapseObject::normalBootstrap() {

    if(covariates){
        y.clear();
        ycenter.clear();

        std::vector<VectorXd> residuals = groupwiseShuffleWithoutReplacement(ycenter_original);

        for (int i = 0; i < residuals.size(); i++)
            y.push_back(y_original[i] - ycenter_original[i] + residuals[i]);

        VectorXd shuffleY = concatenate(y);

        VectorXd beta = getBeta(shuffleY, Z, "normal");
        std::vector<VectorXd> ybar = fitModel(beta, y, z, "normal");
        for (int i = 0; i < y.size(); i++)
            ycenter.push_back(y[i].array() - ybar[i].array());

    }
    else{

        y = groupwiseShuffleWithoutReplacement(y_original);

        double ybar = average(y);
        ycenter.clear();

        for (int i = 0; i < y.size(); i++)
            ycenter.push_back(y[i].array() - ybar);
    }
}

