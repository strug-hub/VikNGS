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