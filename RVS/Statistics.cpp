#include "stdafx.h"
#include "RVS.h"

#include <vector>

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