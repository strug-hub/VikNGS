#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <math.h> 

/*
Computes test statistic for common association test.

@param snp A single of SNP.
@param sample Vector with sample information.
@param group Vector with group information.
@param rvs Indicates how to calculate variance
@return test statistic for snp
*/
double testStatistic(SNP snp, std::vector<Sample> &sample, std::vector<Group> &group, bool rvs) {
	
	size_t i, j, k;
	double v;
	double score;
	double ybar;
	double temp;

	score = 0;
	v = 0;
	ybar = meanY(sample, snp);

	for (i = 0; i < sample.size(); i++) {
		if (!isnan(snp.EG[i]))
			score += (sample[i].y - ybar) * snp.EG[i];
	}

	for (i = 0; i < group.size(); i++) {
		temp = 0;

		for (j = 0; j < group[i].index.size(); j++) {
			k = group[i].index[j];
			if (!isnan(snp.EG[k])) {
				temp += pow(sample[k].y - ybar, 2);
			}
		}

		if (rvs && group[i].hrg)
			v += temp * calcRobustVar(snp.p);
		else
			v += temp * varX(snp, group[i]);
	}

	return pow(score, 2) / v;
	
}


/*
RVS using asymptotic distribution for score test statistic. This functions includes
two association tests which differ in the estimation of the variance of score. It uses
var_case(E(G | D)) when rvs = false and var_case(G) when rvs = true.

Reference: http://www.ncbi.nlm.nih.gov/pubmed/22570057
Reference: http://www.ncbi.nlm.nih.gov/pubmed/24733292
@param snps Vector of SNPs.
@param sample Vector with sample information.
@param group Vector with group information.
@param rvs Indicates how to calculate variance
@return p-values for the SNPs
*/
std::vector<double> RVSasy(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group, bool rvs) {

	std::vector<double> pvals;
	SNP snp;
	double tobs;

	for (size_t j = 0; j < snps.size(); j++) {
		snp = snps[j];
		tobs = testStatistic(snp, sample, group, rvs);
		pvals.push_back(chiSquareOneDOF(tobs));
	}

	return pvals;
}


double btrapHelper(int nboot, std::vector<std::vector<double>> &X, std::vector<std::vector<double>> &Y, 
					std::vector<int> &n, double tobs, double ybar) {

	std::vector<double> rand;
	double score;
	double var;
	double totalVar;
	double temp;
	double tcount;
	double bootcount = 0;

	size_t i, j, k;


	tcount = 0;
	for (int k = 0; k < nboot; k++) {
		score = 0;
		totalVar = 0;

		for (i = 0; i < X.size(); i++) {

			rand = randomSample(X[i], n[i]);
			var = variance(rand);

			for (j = 0; j < rand.size(); j++) {
				temp = Y[i][j] - ybar;
				score += temp * rand[j];
				totalVar += std::pow(temp, 2) * var;
			}
		}

		if (std::pow(score, 2) / totalVar >= tobs)
			tcount++;

		bootcount++;
		/*
		if (earlyStop && bootcount > 1000) {
			pstar = a / ((bootcount + c)*(1 + c));

			if (bootScoreCount / bootcount > pstar) {
				std::cout << "early stop";
				break;
			}
		}
		*/
	}

	return (tcount + 1) / (nboot + 1);
}

/*
Uses RVS to test associaton by bootstrap, given phenotype, expected values of genotypes,
estimated genotype frequency and number of bootstrap iterations.

@param snps Vector of SNPs.
@param sample Vector with sample information.
@param nboot Number of bootstrap iterations.
@param earlyStop Stop bootstrapping with p-value is large.
@param rvs Indicates how to calculate variance of score statistic.
@return Vector of p-values for the SNPs.
*/
std::vector<double> RVSbtrap(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group, int nboot, bool earlyStop, bool rvs) {
	std::vector<double> pvals;
	SNP snp;

	double tobs;
	double xbar;
	double ybar;
	double pval;
	int count;

	size_t i, j, k, l;

	for (size_t k = 0; k < snps.size(); k++) {
		snp = snps[k];
		tobs = testStatistic(snp, sample, group, rvs);
		ybar = meanY(sample, snp);

		std::vector<std::vector<double>> X;
		std::vector<std::vector<double>> Y;
		std::vector<int> n;

		//create vectors to sample from (x - xbar)
		for (i = 0; i < group.size(); i++) {
			xbar = meanX(snp, group[i]);

			std::vector<double> x;
			std::vector<double> y;

			count = 0;
			for (l = 0; l < group[i].index.size(); l++) {
				j = group[i].index[l];
				if (!isnan(snp.EG[j])) {
					x.push_back(snp.EG[j] - xbar);
					y.push_back(sample[j].y);
					count++;
				}
			}

			X.push_back(x);
			Y.push_back(y);
			n.push_back(count);
		}

		//start bootstrapping
		pval = btrapHelper(nboot, X, Y, n, tobs, ybar);
		pvals.push_back(pval);

	}

	return pvals;
}

