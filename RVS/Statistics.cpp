#include "stdafx.h"
#include "RVS.h"

#include <vector>
#include <random>

std::random_device random;
std::mt19937 generate(random());


double meanX(SNP &snp, Group &group) {

	double sum = 0;
	double n = 0;
	size_t i, j;

	for (i = 0; i < group.index.size(); i++) {
		 j = group.index[i];
		if (!isnan(snp.EG[j])) {
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
		if (!isnan(snp.EG[i])) {
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
		if (!isnan(snp.EG[group.index[i]])) {
			mean += snp.EG[group.index[i]];
			n++;
		}
	}

	mean = mean / n;
	double var = 0;

	for (i = 0; i < group.index.size(); i++) {
		if (!isnan(snp.EG[group.index[i]])) {
			var += pow(snp.EG[group.index[i]] - mean, 2);
		}
	}

	return var / (n - 1);
}

double var(VectorXd &X) {

	double mean = 0;
	double n = X.rows();
	int i;

	for (i = 0; i < n; i++)
		mean += X[i];

	mean = mean / n;
	double var = 0;

	for (i = 0; i < X.rows(); i++)
		var += pow(X[i] - mean, 2);

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


std::vector<double> randomSample(std::vector<double> &vec, int nsample) {
	std::vector<double> rvec;
	std::uniform_int_distribution<> sample(0, vec.size()-1);

	for (size_t i = 0; i < nsample; i++)
		rvec.push_back(vec[sample(generate)]);
	
	return rvec;
}

int generateRandomNumber(int from, int to) {
	std::uniform_int_distribution<> sample(from, to);
	return sample(generate);
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

VectorXd CovariateRegression(VectorXd &Y, MatrixXd &Z) {
	return Z.householderQr().solve(Y);
}

VectorXd fitModel(VectorXd &beta, MatrixXd &Z, std::string distribution) {

	double fitted;
	VectorXd meanValue(Z.rows());

	for (size_t i = 0; i < Z.rows(); i++) {
		fitted = 0;
		for (size_t j = 0; j < beta.rows(); j++) {

			fitted += beta[j] * Z(i, j);

			if (distribution == "binom")
				fitted = 1 / (1 + exp(-fitted));
		}
		meanValue[i] = fitted;
	}

	return meanValue;
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
			if (!isnan(x[i])) {
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
			if (!isnan(x[i])) {

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
