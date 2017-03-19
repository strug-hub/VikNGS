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
			var += pow((snp.EG[group.index[i]] - mean), 2);
		}
	}

	return var / (n - 1);
}