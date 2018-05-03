#include "CommonTestObject.h"

//todo: what if all observations in a group are filtered out?????
void CommonTestObject::filterNAN() {
    if(covariates){
        for (int i = 0; i < size(); i++) {
            VectorXd toRemove = whereNAN(x[i], y[i], z[i]);
            x[i] = extractRows(x[i], toRemove, 0);
            y[i] = extractRows(y[i], toRemove, 0);
            z[i] = extractRows(z[i], toRemove, 0);
        }
    }
    else{
        for (int i = 0; i < size(); i++) {
            VectorXd toRemove = whereNAN(x[i], y[i]);
            x[i] = extractRows(x[i], toRemove, 0);
            y[i] = extractRows(y[i], toRemove, 0);
        }
    }
}

void CommonTestObject::setFamily(std::string fam){
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

double CommonTestObject::getVariance(bool rvs) {

    if(covariates){
        double var = covariateVarianceCalculation(rvs);
        double cov = covariateCovarianceCalculation();
        return var + 2*cov;
    }

    double var = 0;

    if(normal){
        double var_y = variance(y);

        for (int i = 0; i < size(); i++) {
            if (rvs && readDepth[i] == 1)
                var += this->y[i].rows() * var_y * robustVar;
            else
                var += y[i].rows() * var_y * variance(x[i]);
        }
        return var;
    }

    //don't replace this later, use for comparsion
    if(binomial){

        for (int i = 0; i < size(); i++) {

            double var_y = ycenter[i].array().pow(2).sum();
            if (rvs && readDepth[i] == 1)
                var += var_y * robustVar;
            else
                var += var_y * variance(x[i]);
        }
        return var;
    }

    throwError("common test", "Don't know how to compute variance, this error should not happen.");
}

void CommonTestObject::bootstrap() {
    std::vector<VectorXd> newxboot;
    int length;

    for (int i = 0; i < size(); i++) {
        length = x[i].rows();
        VectorXd xrand(length);
        for (int j = 0; j < length; j++)
            xrand[j] = x[i][randomInt(0, length - 1)];

        newxboot.push_back(xrand);

    }
    this->xboot = newxboot;
}

double CommonTestObject::bootVariance(int i) { return ycenter[i].array().pow(2).sum() * variance(xboot[i]); }

void CommonTestObject::constructHatMatrix(VectorXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::vector<VectorXd> &yhat){
    VectorXd toRemove = whereNAN(X, Y, Z);
    MatrixXd Z_filtered = extractRows(Z, toRemove, 0);
    VectorXd G_filtered = extractRows(G, toRemove, 0);
    VectorXd Y_filtered = extractRows(Y, toRemove, 0);

    if(this->normal){

        this->H = calculateHatMatrix(Z_filtered);

        double constant = variance(Y_filtered);

        for (int i = 0; i < size(); i++)
                this->yvar.push_back(MatrixXd::Constant(x[i].rows(), 1, constant));
    }
    else if(this->binomial){
        VectorXd mu = concat(yhat);
        VectorXd mu_ = mu.array() * (VectorXd::Ones(mu.rows()) - mu).array();
        VectorXd mu_sqrt = mu_.array().sqrt();

        MatrixXd W = mu_.asDiagonal();
        MatrixXd sqrtW = mu_sqrt.asDiagonal();

        this->H = calculateHatMatrix(Z_filtered, W, sqrtW);

        for (int i = 0; i < size(); i++)
            this->yvar.push_back(extractRows(mu_, G_filtered, i).transpose());
    }

    VectorXd Hdiag = H.diagonal();

    for (int i = 0; i < size(); i++){
        this->h.push_back(extractRows(H, G_filtered, i).transpose());
        this->hdiag.push_back(extractRows(Hdiag, G_filtered, i));
    }

    this->Yvar = concat(this->yvar);
}

double CommonTestObject::covariateCovarianceCalculation(){

    double cov1 = 0;
    for(int j = 0; j < H.cols() - 1; j++)
        for(int k = j+1; k < H.cols(); k++)
            cov1 += (H.col(j).array() * H.col(k).array() * Yvar.array()).sum();

    double cov2 = (H*Yvar).sum() - H.diagonal().dot(Yvar);

    double meanx = average(this->x);
    return meanx * meanx * (cov1 - cov2);
}


double CommonTestObject::covariateVarianceCalculation(bool rvs){

    double var = 0;
    double meanx = average(x);
    double meanxsquared = meanx * meanx;

    for(int i = 0; i < size(); i++){

        double group_var = 0;

        for(int j = 0; j < h[i].cols(); j++){
            VectorXd row = h[i].col(j).array().pow(2);
            group_var += row.dot(Yvar);
        }

        group_var += yvar[i].sum() - 2*(hdiag[i].dot(yvar[i]));

        if (rvs && readDepth[i] == 1)
            var += group_var * (robustVar + meanxsquared);
        else
            var += group_var * (variance(x[i]) + meanxsquared);
    }

    return var;
}

