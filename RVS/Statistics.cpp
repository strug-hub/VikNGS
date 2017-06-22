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

MatrixXd nanToZero(MatrixXd &M) {

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

MatrixXd covariance(MatrixXd &M) {

	//TODO:something better??
	if (M.rows() <= 1) {
		std::cout << "Cannot compute covariate matrix with 1 row!!";
		return M;
	}

	MatrixXd centered = M.rowwise() - M.colwise().mean();
	MatrixXd cov = (centered.adjoint() * centered) / double(M.rows() - 1);

	return cov;
}

MatrixXd correlation(MatrixXd &M) {
	std::cout << M;

	MatrixXd cov = covariance(M);
	VectorXd D = cov.diagonal().array().sqrt();
	
	MatrixXd cor = D.asDiagonal().inverse() * cov * D.asDiagonal().inverse();

	return cor;
}