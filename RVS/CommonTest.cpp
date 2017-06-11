#include "stdafx.h"
#include "RVS.h"
#include "TestSet.h"

#include <iostream>  
#include <math.h> 
#include "Eigen/Dense"

/*
RVS using asymptotic distribution for score test statistic. This functions includes
two association tests which differ in the estimation of the variance of score. It uses
var_case(E(G | D)) when rvs = false and var_case(G) when rvs = true.

Reference: http://www.ncbi.nlm.nih.gov/pubmed/22570057
Reference: http://www.ncbi.nlm.nih.gov/pubmed/24733292
@param CommonTestObject Test object.
@param rvs Indicates how to calculate variance
@return p-value
*/
double commonAsymptotic(TestSet &ts, bool rvs) {
	int i;
	double variance = 0;
	double score = 0;

	for (i = 0; i < ts.length(); i++) {
		score += ts[i].score();
		variance += ts[i].variance(rvs);
	}

	return chiSquareOneDOF(pow(score, 2) / variance);
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

double commonBootstrap(TestSet &ts, int nboot, bool earlyStop, bool rvs) {
	int i;
	double variance = 0;
	double score = 0;

	for (i = 0; i < ts.length(); i++) {
		score += ts[i].score();
		variance += ts[i].variance(rvs);
	}

	double tobs = pow(score, 2) / variance;

	double tcount = 0;
	double average = 0;

	for (size_t n = 1; n <= nboot; n++) {
		variance = 0;
		score = 0;

		for (i = 0; i < ts.length(); i++) {
			ts[i].bootstrapX();
			score += ts[i].bootScore();
			variance += ts[i].bootVariance();
		}

		if (tobs <= pow(score, 2) / variance)
			tcount++;
	
		average += score;
	}

	std::cout << average / nboot;
	std::cout << '\n';


	return (tcount + 1) / (nboot + 1);



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

std::vector<double> runCommonTest(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group,
	int nboot, bool rvs) {
	std::vector<double> pvals;

	for (size_t i = 0; i < snps.size(); i++) {
		TestSet ts(snps[i], sample, group);

		if (nboot == 0)
			pvals.push_back(commonAsymptotic(ts, rvs));
		else
			pvals.push_back(commonBootstrap(ts, nboot, false, true));
	}

	return pvals;

}
