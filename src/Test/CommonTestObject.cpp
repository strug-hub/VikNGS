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

double CommonTestObject::getVarianceRegular() {

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

    double sum1 = 0;
    int index = 0;

    if(covariates){

        int ncov = z[0].cols();
        VectorXd sum2 = VectorXd::Constant(ncov, 0);
        MatrixXd sum3 = MatrixXd::Constant(ncov, ncov, 0);

        for(int i = 0; i < size(); i++){
            for(int j = 0; j < x[i].rows(); j++){

                sum1 += var1[index] * x[i][j] * x[i][j];

                VectorXd two = (var1[index] * x[i][j]) * z[i].row(j).array();
                sum2 += two;

                MatrixXd three = var1[index] * (z[i].row(j).transpose() * z[i].row(j)).array();
                sum3 += three;

                index++;
            }
        }

        double var = (sum2.transpose() * sum3.inverse() * sum2)(0);
        return sum1 - var;
    }
    else{

        double sum2 = 0;
        double sum3 = 0;
        for(int i = 0; i < size(); i++){
            for(int j = 0; j < x[i].rows(); j++){

                sum1 += var1[index] * (x[i].row(j) * x[i].row(j))[0];
                sum2 += var1[index] * x[i].row(j)[0];
                sum3 += var1[index];

                index++;
            }
        }

        double var = (sum2 * (1.0/sum3) * sum2);
        return sum1 - var;
    }

    //should never happen
    return -1;

    /* todo for later

    std::vector<VectorXd> hess1;
    for (int i = 0; i < size(); i++) {
        VectorXd h;
        h = (this->mu[i].array() * (1 - mu[i].array()) );
        hess1.push_back(h);
    }

    if(covariates){

        int ncov = z[0].cols();
        MatrixXd aaRegular = MatrixXd::Constant(ncov, ncov, 0);
        VectorXd abRegular = VectorXd::Constant(ncov, 0);
        double bbRegular=0;

        for (int i = 0; i < size(); i++) {
            for (int j = 0; j < hess1[i].rows(); j++) {

                VectorXd zrow = z[i].row(j);
                MatrixXd aa =  zrow * zrow.transpose();
                aa = aa.array() * hess1[i][j];
                aaRegular = aaRegular + aa;

                VectorXd ab = zrow.array() * (hess1[i][j] * x[i][j]);
                abRegular = abRegular + ab;

                bbRegular += hess1[i][j] * x[i][j] * x[i][j];
            }
        }

        double var = (abRegular.transpose() * aaRegular.inverse() * abRegular)(0);
        return bbRegular - var;

    }
    else{
        double aaRegular = 0;
        double abRegular = 0;
        double bbRegular = 0;

        for (int i = 0; i < size(); i++) {
            for (int j = 0; j < hess1[i].rows(); j++) {
                aaRegular += hess1[i][j];
                abRegular += hess1[i][j] * x[i][j];
                bbRegular += hess1[i][j] * x[i][j] * x[i][j];
            }
        }
        double var = bbRegular - abRegular * (1.0/aaRegular) * abRegular;
        return var;
    }

    */

    return 0;

}

double CommonTestObject::getVarianceBinomial(bool rvs) {

    if(useHatMatrix && covariates){
        double var = covariateVarianceCalculation(rvs);
        double cov = covariateCovarianceCalculation();
        return var + 2*cov;
    }

    double var = 0;
    for (int i = 0; i < size(); i++) {

        double var_y = ycenter[i].array().pow(2).sum();
        if (rvs && readDepth[i] == 1)
            var += var_y * robustVar;
        else
            var += var_y * variance(x[i]);
    }
    return var;

}

double CommonTestObject::getVarianceNormal(bool rvs) {

    if(covariates){

        double var = 0;
        double var_y = 0;
        double n = 0;
        for (int i = 0; i < size(); i++) {
            var_y += ycenter[i].array().pow(2).sum();

            double groupSize = x[i].rows();
            n+=groupSize;

            if (rvs && readDepth[i] == 1)
                var += groupSize * robustVar;
            else
                var += groupSize * variance(x[i]);
        }

        return (1.0/n) * var_y * var;
    }

    double var_y = variance(y);
    double var = 0;
    for (int i = 0; i < size(); i++) {
        if (rvs && readDepth[i] == 1)
            var += this->y[i].rows() * var_y * robustVar;
        else
            var += y[i].rows() * var_y * variance(x[i]);
    }
    return var;

}

void CommonTestObject::bootstrap() {

    //todo: shuffle Z if has covariates? OR remove completely
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



// hat matrix logic
// v v v v v v v v

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
        VectorXd mu = concatenate(yhat);
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

    this->Yvar = concatenate(this->yvar);
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

