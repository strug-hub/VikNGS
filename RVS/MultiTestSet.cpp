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


VectorXd MultiTestSet::getYmHRG() {
	VectorXd Ym(length());

	for (size_t i = 0; i < length(); i++)
			Ym[i] = testset[i].getYmHRG();

	return Ym;
}

VectorXd MultiTestSet::getYmLRG() {
	VectorXd Ym(length());

	for (size_t i = 0; i < length(); i++)
		Ym[i] = testset[i].getYmLRG();

	return Ym;
}


VectorXd MultiTestSet::getScoreVector() {
	VectorXd score(length());

	for (size_t i = 0; i < length(); i++)
		score[i] = testset[i].getScore();

	return score;
}

VectorXd MultiTestSet::getBootstrapScoreVector() {
	VectorXd score(length());

	for (size_t i = 0; i < length(); i++)
		score[i] = testset[i].getBootstrapScore();

	return score;
}

MatrixXd MultiTestSet::getXMatrix(int index) {

	MatrixXd X(testset[0].get(index).length_filterz(), length());

	for (size_t i = 0; i < length(); i++)
		X.col(i) = testset[i].get(index).getX_filterz();

	return X;
}


MatrixXd MultiTestSet::getBootstrapXMatrix(int index) {
	MatrixXd X(testset[0].get(index).length_filterz(), length());

//	for (size_t i = 0; i < length(); i++)
//		X.col(i) = testset[i].get(index).getBootstrapX_filterz();

	return X;
}
