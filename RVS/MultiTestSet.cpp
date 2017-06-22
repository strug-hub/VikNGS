#include "stdafx.h"
#include "RVS.h"
#include "MultiTestSet.h"

MultiTestSet::MultiTestSet(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group) {

	for (size_t i = 0; i < snps.size(); i++)
		this->testset.push_back(TestSet(snps[i], sample, group, true));

}

VectorXd MultiTestSet::getRobustVariance() {

	VectorXd var(length());
	for (size_t i = 0; i < length(); i++)
		var[i] = testset[i].getRobustVariance();

	return var;
}


VectorXd MultiTestSet::getYm(bool hrg) {
	VectorXd Ym(length());

	for (size_t i = 0; i < length(); i++)
		if (hrg)
			Ym[i] = testset[i].getYmHRG();
		else
			Ym[i] = testset[i].getYmLRG();

	return Ym;
}

std::vector<double> MultiTestSet::calculateYbar()
{
	size_t i, j, k;
	std::vector<double> yhrg;
	std::vector<double> ylrg;

	TestSet ts;
//	TestGroup tg;

	double temp;

	for (i = 0; i < length(); i++) {
		ts = testset[i];
		temp = 0;
		for (j = 0; j < ts.length(); j++) {
//			tg = ts.groups[i];

			//	for (k = 0; k < tg.length(); k++)
			//
			//	variance += pow(Y[i] - Ybar[i], 2);



		}
	}


	return std::vector<double>();
}


MatrixXd MultiTestSet::getCovariateMatrix(int index) {

	MatrixXd X(testset[0].groupLRG[index].length_filterz(), length());

	for (size_t i = 0; i < length(); i++) 
		X.col(i) = testset[i].groupLRG[index].getX_filterz();

	return covariance(X);
}

MatrixXd MultiTestSet::getCorrelationMatrix(int index) {

	MatrixXd X(testset[0].groupHRG[index].length_filterz(), length());

	for (size_t i = 0; i < length(); i++)
		X.col(i) = testset[i].groupHRG[index].getX_filterz();

	return correlation(X);
}