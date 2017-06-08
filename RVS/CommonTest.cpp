#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <math.h> 

#include "Eigen/Dense"


class CommonTestObject {

	private:
		std::vector<bool> H;

	public:
		VectorXd Y;
		VectorXd X;
		MatrixXd Z;
		VectorXi G;
		int length;
		int ngroup;
		std::vector<std::vector<double>> P;

		//assumes group labels are valid!!
		CommonTestObject(SNP &snp, std::vector<Sample> &sample, std::vector<Group> &group) {
			
			VectorXd x(sample.size());
			VectorXd y(sample.size());
			VectorXi g(sample.size());
			MatrixXd z(sample.size(), sample[0].covariates.size() + 1);
			std::vector<bool> hrg(group.size());
			std::vector<std::vector<double>> p;

			int c = 0;
			bool flag;

			for(size_t i = 0; i< group.size(); i++)
				hrg[group[i].groupIndex] = group[i].hrg;

			for (size_t i = 0; i < sample.size(); i++) {
				if (!isnan(snp.EG[i]) && !isnan(sample[i].y)) {

					y(c) = sample[i].y;
					p.push_back(snp.p);
					x(c) = snp.EG[i];
					g(c) = sample[i].groupIndex;
					
					flag = false;

					for (size_t j = 0; j < sample[i].covariates.size(); j++) {
						z(c, j + 1) = sample[i].covariates[j];

						if (isnan(sample[i].covariates[j])) {
							flag = true;
							break;
						}
					}

					if (!flag) {
						z(c, 0) = 1;
						c++;
					}
				}
			}
			length = c;
			ngroup = group.size();

			Y = y.block(0, 0, c, 1);
			X = x.block(0, 0, c, 1);
			Z = z.block(0, 0, c, sample[0].covariates.size() + 1);
			G = g.block(0, 0, c, 1);
			H = hrg;
			P = p;
		}

		inline bool isHRG(int groupID) {
			return(H[groupID]);
		}
};

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
double commonAsymptotic(CommonTestObject &o, bool rvs) {
	int i;
	double v = 0;
	double score = 0;

	std::vector<double> meanValues = CovariateRegression(o.Y, o.Z);
	std::vector<double> xvar;
	std::vector<double> vs;

	for (i = 0; i < o.ngroup; i++){
		vs.push_back(0);
		
		if (rvs && o.isHRG(i))
			xvar.push_back(calcRobustVar(o.P[i]));
		else
			xvar.push_back(varX(o.X, o.G, i));
	}

	for (i = 0; i < o.length; i++) {
		vs[o.G(i)] += pow(o.Y(i) - meanValues[i], 2);
		score += (o.Y(i) - meanValues[i]) * o.X(i);
	}

	for (i = 0; i < o.ngroup; i++)
		v += vs[i] * xvar[i];
	

	return chiSquareOneDOF(pow(score, 2) / v);
}

std::vector<double> runCommonTest(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group,
									int nboot, bool rvs) {
	std::vector<double> pvals;

	for (size_t i = 0; i < snps.size(); i++) {
		CommonTestObject o(snps[i], sample, group);

		if (nboot == 0)
			pvals.push_back(commonAsymptotic(o, rvs));
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

	size_t i, j;


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
		tobs = 0;// testStatistic(snp, sample, group, rvs);
		ybar = 0;// meanY(sample, snp);

		std::vector<std::vector<double>> X;
		std::vector<std::vector<double>> Y;
		std::vector<int> n;

		//create vectors to sample from (x - xbar)
		for (i = 0; i < group.size(); i++) {
			xbar = 0; // meanX(snp, group[i]);

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

