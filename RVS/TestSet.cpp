#include "stdafx.h"
#include "RVS.h"
#include "TestSet.h"

TestSet::TestSet(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group, bool rare) {
	size_t i;
	int ncov = sample[0].covariates.size();

	for (i = 0; i < group.size(); i++)
		if (group[i].hrg)
			this->groupHRG.push_back(TestHRG(snp, sample, group[i], rare));
		else 
			this->groupLRG.push_back(TestLRG(snp, sample, group[i], rare));

	this->nhrg = 0;
	this->nlrg = 0;
	for (i = 0; i < groupHRG.size(); i++)
		this->nhrg += groupHRG[i].length();
	for (i = 0; i < groupLRG.size(); i++)
		this->nlrg += groupLRG[i].length();

	if (rare) {
		this->nhrg_filterz = 0;
		this->nlrg_filterz = 0;
		for (i = 0; i < groupHRG.size(); i++)
			this->nhrg_filterz += groupHRG[i].length_filterz();
		for (i = 0; i < groupLRG.size(); i++)
			this->nlrg_filterz += groupLRG[i].length_filterz();
	}

	//this->robustVar = calcRobustVar(snp.p);

	VectorXd beta = getBeta(snp, sample, group);
	for (i = 0; i < groupHRG.size(); i++)
		groupHRG[i].fitModel(beta, "normal");
	for (i = 0; i < groupLRG.size(); i++) 
		groupLRG[i].fitModel(beta, "normal");
}

VectorXd TestSet::getBeta(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group) {
	size_t i, j, k;
	int ncov = sample[0].covariates.size();
	int nobs = snp.EG.size();


	VectorXd y(nobs);
	MatrixXd z(nobs, ncov + 1);

	std::vector<bool> toRemove(nobs, false);

	for (i = 0; i < nobs; i++) {
		if (isnan(snp.EG[i]) || isnan(sample[i].y))
			toRemove[i] = true;
		else {
			for (j = 0; j < ncov; j++)
				if (isnan(sample[i].covariates[j]))
					toRemove[i] = true;
		}
	}

	int c = 0;
	for (i = 0; i < nobs; i++) {
		if (!toRemove[i]) {
			y[c] = sample[i].y;
			for (j = 0; j < ncov; j++)
				z(c, j + 1) = sample[i].covariates[j];

			z(c, 0) = 1;
			c++;
		}
	}

	VectorXd yblock = y.block(0, 0, c, 1);
	MatrixXd zblock = z.block(0, 0, c, ncov + 1);

	return CovariateRegression(yblock, zblock);
}

//---------------------------------
//	For RareTest
//---------------------------------

double TestSet::getYmLRG() {
	double Ym = 0;

	for (size_t i = 0; i < groupLRG.size(); i++)
		Ym += groupLRG[i].Ycenter.array().pow(2).sum();

	Ym = Ym / nlrg * nlrg_filterz;
	return sqrt(Ym);
}

double TestSet::getYmHRG() {
	double Ym = 0;

	for (size_t i = 0; i < groupHRG.size(); i++)
		Ym += groupHRG[i].Ycenter.array().pow(2).sum();

	Ym = Ym / nhrg * nhrg_filterz;
	return sqrt(Ym);
}

double TestSet::getScore() {
	double s = 0;

	for (size_t i = 0; i < length(); i++)
		s += get(i).score();

	return s;
}

double TestSet::getBootstrapScore() {
	double s = 0;

//	for (size_t i = 0; i < length(); i++)
//		s += get(i).bootstrapScore();

	return s;

}
