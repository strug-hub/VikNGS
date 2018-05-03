#include "MathHelper.h"

static const std::string STATISTICS_HELPER = "statistics helper";

inline VectorXd CovariateNormalRegression(VectorXd &Y, MatrixXd &Z) {
	return Z.householderQr().solve(Y);
}

inline VectorXd CovariateBinomialRegression(VectorXd &Y, MatrixXd &Z) {
    return logisticRegression(Y,Z);
}

VectorXd getBeta(VectorXd &X, VectorXd &Y, MatrixXd &Z, std::string family) {

	VectorXd toRemove = whereNAN(X, Y, Z);
	VectorXd y_filtered = extractRows(Y, toRemove, 0);
    MatrixXd z_filtered = extractRows(Z, toRemove, 0);

    if(family == "binomial")
        return CovariateBinomialRegression(y_filtered, z_filtered);
    if(family == "normal")
        return CovariateNormalRegression(y_filtered, z_filtered);

    throwError(STATISTICS_HELPER, "Family not recognized, could not compute beta values.", family);
}

std::vector<VectorXd> fitModel(VectorXd &beta, std::vector<VectorXd> &y, std::vector<MatrixXd> &z, std::string family) {
    std::vector<VectorXd> yhat;

    for (int i = 0; i < y.size(); i++) {
	
        VectorXd meanValue = z[i] * beta;

        if (family == "binomial")
            meanValue = 1 / (1 + exp(-meanValue.array()));

        //if (family == "normal")
            //do nothing

        yhat.push_back(meanValue);
	}

    return yhat;
}

double calcRobustVar(VectorXd &p) {
	return (4 * p[2] + p[1]) - pow(2 * p[2] + p[1], 2);
}

std::random_device rd;
std::mt19937 generate(random());

int randomInt(int from, int to) {
	std::uniform_int_distribution<> sample(from, to);
	return sample(generate);
}
double randomDouble(double from, double to) {
	std::uniform_real_distribution<double> sample(from, to);
	return sample(generate);
}
double randomNormal(double mean, double sd) {
	std::normal_distribution<> sample(mean, sd);
	return sample(generate);
}

//same as doing pairwise.complete.obs in R
MatrixXd covariance(MatrixXd &M) {
	int n = M.cols();
	int m = M.rows();

	MatrixXd cov(n, n);
	size_t i, j, k;

	double count;
	double meani;
	double meanj;
	double sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {

			sum = 0;
			count = 0;
			meani = 0;
			meanj = 0;

			for (k = 0; k < m; k++) {
				if (!std::isnan(M(k, i)) && !std::isnan(M(k, j))) {
					count++;
					meani += M(k, i);
					meanj += M(k, j);
				}
			}

			meani /= count;
			meanj /= count;

			for (k = 0; k < m; k++)
				if (!std::isnan(M(k, i)) && !std::isnan(M(k, j)))
					sum += (M(k, i) - meani) * (M(k, j) - meanj);

			sum /= count - 1;
			cov(i, j) = sum;
			cov(j, i) = sum;
		}
	}

	return cov;
}


MatrixXd correlation(MatrixXd &M) {
	int n = M.cols();
	int m = M.rows();

	MatrixXd cor(n, n);
	size_t i, j, k;

	double count;
	double vari;
	double varj;
	double meani;
	double meanj;
	double sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			if (i == j) {
				cor(i, j) = 1;
				cor(j, i) = 1;
			}
			else {
				sum = 0;
				count = 0;
				vari = 0;
				varj = 0;
				meani = 0;
				meanj = 0;

				for (k = 0; k < m; k++) {
					if (!std::isnan(M(k, i)) && !std::isnan(M(k, j))) {
						count++;
						meani += M(k, i);
						meanj += M(k, j);
					}
				}

				meani /= count;
				meanj /= count;

				for (k = 0; k < m; k++) {
					if (!std::isnan(M(k, i)) && !std::isnan(M(k, j))) {
						vari += pow(M(k, i) - meani, 2);
						varj += pow(M(k, j) - meanj, 2);
						sum += (M(k, i) - meani) * (M(k, j) - meanj);
					}
				}

				sum /= sqrt(vari * varj);
				cor(i, j) = sum;
				cor(j, i) = sum;
			}
		}
	}

	return cor;
}

double pnorm(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p*x);

	return (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
}

std::vector<double> randomSample(std::vector<double> &vec, int nsample) {
	std::vector<double> rvec;
	std::uniform_int_distribution<> sample(0, vec.size() - 1);

	for (size_t i = 0; i < nsample; i++)
		rvec.push_back(vec[sample(generate)]);

	return rvec;
}

double chiSquareOneDOF(double statistic) {
	double z = statistic * 0.5;
	double sc = 2 * sqrt(z) * exp(-z);

	double sum = 1;
	double prevSum = sum;
	double nom = 1;
	double dnom = 1;
	double s = 0.5;

	for (int i = 0; i < 200; i++)
	{
		nom *= z;
		s++;
		dnom *= s;
		sum += nom / dnom;
		if (prevSum == sum) break;
		prevSum = sum;
	}

	double p = sum * sc;
	if (std::isnan(p) || std::isinf(p) || p <= 1e-8)
		return 1e-14;

	p /= tgamma(0.5);

	return std::max(1 - p, 1e-14);
}


inline double sigmoid(VectorXd &x, VectorXd &beta){
    return (1.0 / (1.0 + exp(-(beta.dot(x)))));
}

MatrixXd calculateHatMatrix(MatrixXd &Z){
    try{
        return Z*(Z.transpose()*Z).inverse()*Z.transpose();
    }catch(...){
        throwError(STATISTICS_HELPER, "Error while calculate covariate hat matrix. Not invertable?");
    }
}

MatrixXd calculateHatMatrix(MatrixXd &Z, MatrixXd &W, MatrixXd &sqrtW){
    try{
        return sqrtW*Z*(Z.transpose()*W*Z).inverse()*Z.transpose()*sqrtW;
    }catch(...){
        throwError(STATISTICS_HELPER, "Error while calculate covariate hat matrix. Not invertable?");
    }
}

VectorXd logisticRegression(VectorXd &Y, MatrixXd &X) {

    VectorXd beta = VectorXd::Constant(X.cols(), 0);
    VectorXd lastBeta = VectorXd::Constant(X.cols(), 0);

    int nobs = X.rows();

    std::vector<VectorXd> x;
    std::vector<MatrixXd> xxT;

    for(int i = 0; i< nobs; i++){
        x.push_back(X.row(i));
        xxT.push_back(X.row(i)*X.row(i).transpose());
    }

    VectorXd first;
    MatrixXd hess;

    int iteration = 0;

    while (true) {

        iteration++;

        VectorXd p(nobs);
        VectorXd W(nobs);

        for(int i =0; i < nobs; i++){
            double s = sigmoid(x[i], beta);
            p[i] = s;
            W[i] = s * (1-s);
        }

        //calculate 1st derivative vector of likelihood function
        first = X.transpose()*(Y - p);

        //calculate Hessian matrix of likelihood function
        hess = -X.transpose()*W.asDiagonal()*X;

        lastBeta = beta;
        try{
        beta = lastBeta - hess.inverse()*first;
        }catch(...){
            printWarning("Error while trying to invert matrix in logistic regression");
            throw variant_exception( "Failed to invert" );
        }

        bool stop = true;
        VectorXd grad = (beta - lastBeta).array().abs();
        for(int i = 0; i< grad.rows(); i++){
            if(grad[i] > 1e-7){
                stop = false;
                break;
            }
        }

        if(iteration > 50 || stop)
            break;
    }

    if(iteration > 50){
        printWarning("Logistic regression failed to converge after 50 iterations.");
        throw variant_exception( "Regression failed to converge" );
    }

    return beta;
}
