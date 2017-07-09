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

int generateRandomInteger(int from, int to) {
	std::uniform_int_distribution<> sample(from, to);
	return sample(generate);
}
double generateRandomDouble(double from, double to) {
	std::uniform_real_distribution<double> sample(from, to);
	return sample(generate);
}
double randomNormal(double mean, double sd) {
	std::normal_distribution<> sample(mean, sd);
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

MatrixXd nanToZero(MatrixXd &M) {
	std::cout << "\n";	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";

	std::cout << M;

	for (size_t i = 0; i < M.cols(); i++) 
		for (size_t j = 0; j < M.rows(); j++) 
			if (isnan(M(j, i)))
				M(j, i) = 0;

	return  M;
}

VectorXd nanToZero(VectorXd &V) {

	for (size_t i = 0; i < V.rows(); i++)
			if (isnan(V[i]))
				V[i] = 0;

	return V;
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
				if (!isnan(M(k,i)) && !isnan(M(k,j))) {
					count++;
					meani += M(k, i);
					meanj += M(k, j);
				}
			}

			meani /= count;
			meanj /= count;

			for (k = 0; k < m; k++) 
				if (!isnan(M(k, i)) && !isnan(M(k, j)))
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
					if (!isnan(M(k,i)) && !isnan(M(k,j))) {
						count++;
						meani += M(k,i);
						meanj += M(k,j);
					}
				}

				meani /= count;
				meanj /= count;

				for (k = 0; k < m; k++) {
					if (!isnan(M(k, i)) && !isnan(M(k, j))) {
						vari += pow(M(k,i) - meani, 2);
						varj += pow(M(k,j) - meanj, 2);
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

/*
Approximates the p-value from the pdf of the normal distribution where x is a Z-score

@param x Z-score.
@return p-value.
*/
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