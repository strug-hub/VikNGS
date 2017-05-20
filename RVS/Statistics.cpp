#include "stdafx.h"
#include "RVS.h"

#include <vector>
#include <random>

double meanX(SNP &snp, Group &group) {

	double sum = 0;
	double n = 0;
	size_t i, j;

	for (i = 0; i < group.index.size(); i++) {
		 j = group.index[i];
		if (snp.EG[j] != NULL) {
			sum += snp.EG[j];
			n++;
		}
	}

	return sum / n;
}


double meanY(std::vector<Sample> &sample, SNP &snp) {

	double sum = 0;
	double n = 0;

	for (size_t i = 0; i < sample.size(); i++) {
		if (snp.EG[i] != NULL) {
			sum += sample[i].y;
			n++;
		}
	}

	return sum / n;
}

double varX(SNP &snp, Group &group) {

	double mean = 0;
	double n = 0;
	size_t i;

	for (i = 0; i < group.index.size(); i++) {
		if (snp.EG[group.index[i]] != NULL) {
			mean += snp.EG[group.index[i]];
			n++;
		}
	}

	mean = mean / n;
	double var = 0;

	for (i = 0; i < group.index.size(); i++) {
		if (snp.EG[group.index[i]] != NULL) {
			var += pow(snp.EG[group.index[i]] - mean, 2);
		}
	}

	return var / (n - 1);
}

double variance(std::vector<double> &vec) {

	double mean = 0;
	size_t i;

	for (i = 0; i < vec.size(); i++)
		mean += vec[i];

	mean = mean / vec.size();
	double var = 0;

	for (i = 0; i < vec.size(); i++)
		var += std::pow(vec[i] - mean, 2);

	return var / (vec.size() - 1);
}

std::random_device rd;
std::mt19937 gen(rd());

std::vector<double> randomSample(std::vector<double> &vec, int nsample) {
	std::vector<double> rvec;
	std::uniform_int_distribution<> sample(0, vec.size()-1);

	for (size_t i = 0; i < nsample; i++)
		rvec.push_back(vec[sample(gen)]);
	
	return rvec;
}


/*
Finds p-value for test statistic using a chi-squared distribution with one degree of freedom
using chi-squared probability density function.

@param statistic Test statistic.
@return p-value.
*/
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
	if (isnan(p) || isinf(p) || p <= 1e-8)
		return 1e-14;

	p /= tgamma(0.5);

	return 1 - p;
}



#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

std::vector<double> CovariateRegression(std::vector<Sample> &sample) {

	//todo: filter out NAs
	/*
	
	Coefficients:
	Estimate Std. Error t value Pr(>|t|)
	(Intercept)  4.043e-01  9.930e-02   4.071 7.46e-05 ***
	Z1          -2.265e-05  1.163e-02  -0.002    0.998
	Z2          -1.767e-01  1.257e-01  -1.406    0.162
	Z3          -6.635e-02  7.245e-02  -0.916    0.361
	
	*/
	size_t i, j;

	VectorXd y(sample.size());
	MatrixXd z(sample.size(), sample[0].numeric_cov.size() + 1);

	for (i = 0; i < sample.size(); i++) {
		y(i) = sample[i].y;
		
		z(i, 0) = 1;

		for (j = 0; j < sample[i].numeric_cov.size(); j++) 
			z(i, j+1) = sample[i].numeric_cov[j];
	}

	VectorXd yy(5);
	MatrixXd zz(5, 3);
	yy(0) = 0;
	yy(1) = 1;
	yy(2) = 0;
	yy(3) = 1;
	yy(4) = 0;

	zz(0, 0) = 1;
	zz(1, 0) = 1;
	zz(2, 0) = 1;
	zz(3, 0) = 1;
	zz(4, 0) = 1;

	zz(0, 1) = 1;
	zz(1, 1) = 2;
	zz(2, 1) = 3;
	zz(3, 1) = 4;
	zz(4, 1) = 5;

	zz(0, 2) = 1;
	zz(1, 2) = 0;
	zz(2, 2) = 1;
	zz(3, 2) = 0;
	zz(4, 2) = 3;

	std::cout << zz;

	std::cout << zz.fullPivHouseholderQr().solve(yy);
	
	std::vector<double> betas;
	return betas;
	
}



std::vector<double> LogisticRegression(std::vector<bool> &y, std::vector<double> &x) {
	double w0 = 0;
	double w = 0;

	double w0last = 1;
	double wlast = 1;

	double temp;
	double sigmoid;
	double det;

	double gradw0;
	double gradw;
	double a;
	double bc;
	double d;

	int iteration = 0;

	while (true) {

		iteration++;

		gradw0 = 0;
		gradw = 0;
		a = 0;
		bc = 0;
		d = 0;

		for (size_t i = 0; i < y.size(); i++) {
			if (x[i] != NULL) {
				sigmoid = (1.0 / (1.0 + exp(-(w * x[i] + w0))));

				//calculate 1st derivative of likelihood function
				temp = sigmoid - y[i];
				gradw0 += temp;
				gradw += temp * x[i];

				//calculate Hessian matrix
				temp = sigmoid * (1 - sigmoid);
				a += temp;
				bc += temp * x[i];
				d += temp * x[i] * x[i];
			}
		}

		temp = a*d - bc*bc;
		if (temp == 0) {
			std::cout << "Warning: Hessian not invertible (logistic regression).\n";
			std::vector<double> out;
			out.push_back(w0);
			out.push_back(w);
		}

		det = 1 / temp;

		w0 -= det * ((d * gradw0) + (-bc * gradw));
		w -= det * ((-bc * gradw0) + (a * gradw));

		if (iteration > 100 || (abs(w - wlast) < 1e-7 && abs(w0 - w0last) < 1e-7 && iteration > 1)) {

			if (iteration > 100)
				std::cout << "Warning: Logistic regression failed to converge after 100 iterations.\n";

			std::vector<double> out;
			out.push_back(w0);
			out.push_back(w);
			return out;
		}

		wlast = w;
		w0last = w0;
	}
}


double LogisticRegressionInterceptOnly(std::vector<bool> &y, std::vector<double> &x) {
	double w0 = 0;
	double w = 0;

	double w0last = 1;

	double sigmoid;

	double gradw0;
	double hessw0;

	int iteration = 0;

	while (true) {

		iteration++;

		gradw0 = 0;
		hessw0 = 0;

		for (size_t i = 0; i < y.size(); i++) {
			if (x[i] != NULL) {

				sigmoid = (1.0 / (1.0 + exp(-w0)));

				//calculate 1st derivative of likelihood function
				gradw0 += sigmoid - y[i];

				//calculate 2nd derivative of likelihood function
				hessw0 += sigmoid * (1 - sigmoid);
			}
		}

		w0 -= gradw0 / hessw0;

		if (iteration > 100 || (abs(w0 - w0last) < 1e-7 && iteration > 1)) {

			if (iteration > 100)
				std::cout << "Warning: Logistic regression failed to converge after 100 iterations.\n";

			return w0;
		}

		w0last = w0;
	}
}
