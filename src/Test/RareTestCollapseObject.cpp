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

    if(normal)
        normalBootstrap();
    else if (binomial)
        binomialBoostrap();

    return;
}

void RareTestCollapseObject::binomialBootstrap() {

    for(int i = 0; i< size(); i++)
        t[i].bootstrap();

    z = shuffleWithoutReplacement(z);

}


void RareTestCollapseObject::normalBootstrap() {

    if(covariates){
        std::vector<VectorXd> residuals = shuffleWithoutReplacement(ycenter_original);
        y.clear();

        for (int i = 0; i < residuals.size(); i++)
            y.push_back(y_original[i] - ycenter_original[i] + residuals[i]);


        VectorXd Y = concatenate(y);

        VectorXd beta = getBeta(Y, Z, "normal");
        ycenter = fitModel(beta, y, z, "normal");
    }
    else{

        y = shuffleWithoutReplacement(y_original);

        double ybar = average(y);
        ycenter.clear();

        for (int i = 0; i < y.size(); i++)
            ycenter.push_back(y[i].array() - ybar);
    }
}

