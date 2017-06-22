#include "stdafx.h"
#include "RVS.h"
#include <iostream>
#include "MultiTestSet.h"
#include "TestSet.h"
#include "TestGroup.h"

/*
Approximates the p-value from the pdf of the normal distribution where x is a Z-score

@param x Z-score.
@return p-value.
*/
double pnorm(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p*x);
	return (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
}

std::vector<double> calcScoreVector(std::vector<SNP> &snps, std::vector<Sample> &sample) {
	SNP snp;
	std::vector<double> score;
	double s;
	double ybar;
	size_t i, j;

	//calculate score vector
	for (i = 0; i < snps.size(); i++) {
		snp = snps[i];
		ybar = meanY(sample, snp);
		s = 0;

		for (j = 0; j < sample.size(); j++)
			s += (sample[j].y - ybar) * snp.EG[j];

		score.push_back(s);
	}

	return score;
}

MatrixXd correlation(std::vector<SNP> &snps, Group group) {
	int n = snps.size();
	MatrixXd cor(n, n);
	size_t i, j, k, l;
	
	double count;
	double vari;
	double varj;
	double meani;
	double meanj;
	double sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			if (i == j) {
				cor(i, j) = 1;
				cor(j, i) = 1;
			}
			else {
				sum = 0;
				count = 0;
				vari = 0;
				varj = 0;
				meani = 0;
				meanj = 0;

				for (l = 0; l < group.index.size(); l++) {
					k = group.index[l];
					if (!isnan(snps[i].EG[k]) && !isnan(snps[j].EG[k])) {
						count++;
						meani += snps[i].EG[k];
						meanj += snps[j].EG[k];
					}
				}

				meani /= count;
				meanj /= count;

				for (l = 0; l < group.index.size(); l++) {
					k = group.index[l];
					if (!isnan(snps[i].EG[k]) && !isnan(snps[j].EG[k])) {
						vari += pow((snps[i].EG[k] - meani), 2);
						varj += pow((snps[j].EG[k] - meanj), 2);
						sum += (snps[i].EG[k] - meani) * (snps[j].EG[k] - meanj);
					}
				}

				sum /= sqrt(vari * varj);
				cor(i, j) = sum;
				cor(j, i) = sum;
			}
		}
	}
	
	return cor;
}

MatrixXd covariance(std::vector<SNP> &snps, Group group) {
	int n = snps.size();
	MatrixXd cov(n, n);
	size_t i, j, k, l;

	double count;
	double meani;
	double meanj;
	double sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {

			sum = 0;
			count = 0;
			meani = 0;
			meanj = 0;

			for (l = 0; l < group.index.size(); l++) {
				k = group.index[l];
				if (!isnan(snps[i].EG[k]) && !isnan(snps[j].EG[k])) {
					count++;
					meani += snps[i].EG[k];
					meanj += snps[j].EG[k];
				}
			}

			meani /= count;
			meanj /= count;

			for (l = 0; l < group.index.size(); l++) {
				k = group.index[l];
				if (!isnan(snps[i].EG[k]) && !isnan(snps[j].EG[k]))
					sum += (snps[i].EG[k] - meani) * (snps[j].EG[k] - meanj);
			}
			sum /= count-1;
			cov(i, j) = sum;
			cov(j, i) = sum;
			
		}
	}

	return cov;
}

std::vector<double> testStatistic(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group, bool rvs) {
	size_t i, j, k, l, m;
	int n = snps.size();

	MatrixXd sigma;
	MatrixXd diagS(n, n);

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			diagS(i, j) = 0;

	double ybar;
	double sum;
	double count;
	std::vector<double> var;

	//TODO: check_hom_matrix

	for (i = 0; i < n; i++)
		var.push_back(sqrt(calcRobustVar(snps[i].p)));

	for (m = 0; m < group.size(); m++) {

		std::vector<double> Ym;

		for (j = 0; j < snps.size(); j++) {
			ybar = meanY(sample, snps[j]);
			sum = 0;
			count = 0;

			for (l = 0; l < group[m].index.size(); l++) {
				k = group[m].index[l];
				if (!isnan(snps[j].EG[k])) {
					sum += std::pow(sample[k].y - ybar, 2);
					count++;
				}
			}
			Ym.push_back(std::sqrt(sum / count * group[m].index.size()));
		}

		if (group[m].hrg && rvs) {
			sigma = correlation(snps, group[m]);
			for (i = 0; i < n; i++) {
				for (j = i; j < n; j++) {
					diagS(i, j) = diagS(i, j) + sigma(i, j) * var[i] * var[j] * Ym[i] * Ym[j];
					diagS(j, i) = diagS(i, j);
				}
			}
		}
		else {
			sigma = covariance(snps, group[m]);
			for (i = 0; i < n; i++) {
				for (j = i; j < n; j++) {
					diagS(i, j) = diagS(i, j) + sigma(i, j) * Ym[i] * Ym[j];
					diagS(j, i) = diagS(i, j);
				}
			}
		}
	}

	VectorXd e = diagS.eigenvalues().real();
	std::vector<double> eigenvals(e.data(), e.data() + e.size());
	
	sum = 0;
	for (i = 0; i < n; i++) 
		for (j = 0; j < n; j++) 
			sum += diagS(i, j);

	double s = 0;
	double s2 = 0;
	std::vector<double> score = calcScoreVector(snps, sample);

	for (i = 0; i < n; i++) {
		s += score[i];
		s2 += score[i] * score[i];
	}
	
	//CAST -> linear test
	double cast = pnorm(s / sqrt(sum));
	//C-alpha -> quadratic test
	double calpha = qfc(eigenvals, s2, snps.size());

	return {cast, calpha};
}

/*
RVS analysis with rare variants for one group. Get p-values from RVS
with CAST and C-alpha (resampling: bootstrap; variance estimate : robust)

@param snps Vector of SNPs.
@param sample Vector with sample information.
@param group Vector with group information.
@param nboot Number of bootstrap iterations.
@param rvs Indicates how to calculate variance of score statistic.
@param njoint Number of SNPs grouping together for one test, default is 5.

checkHomMatrix parameters:
@param hom 1 or 2; 1 means making changes with the 1st non-NAN element, 2 means making changes with a random non-NAN element
@param multiplier Value 1 or 2; 1 is dividing by 2 and 2 is multiplying by 2.

@return Vector with two p-values. First element is linear p-value (CAST), second is quadratic p-value (C-alpha).
*/
std::vector<double> RVSrare(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group, int nboot, bool rvs, int njoint, int hom, int multiplier) {
	
	SNP snp;
	std::vector<double> tobs = testStatistic(snps, sample, group, rvs);
	std::vector<double> tsamp;
	std::vector<SNP> samples = cloneAll(snps);

	//SNPs > groups > EG values
	std::vector<std::vector<std::vector<double>>> X;
	std::vector<double> temp;

	double tcountLinear = 0;
	double tcountQuadratic = 0;
	double bootcount = 0;

	size_t i, j, k, l, h;

	//create vectors to sample from (x - xbar)
	for (i = 0; i < snps.size(); i++) {
		snp = snps[i];
		std::vector<std::vector<double>> xsnp;

		for (j = 0; j < group.size(); j++) {

			std::vector<double> xeg;
			double xbar = meanX(snp, group[j]);

			for (l = 0; l < group[j].index.size(); l++) {
				k = group[j].index[l];
				if (!isnan(snp.EG[k]))
					xeg.push_back(snp.EG[k] - xbar);
				//TODO: ask if NAs should be sampled?? (they aren't in common test)
				else
					xeg.push_back(snp.EG[k]);
			}
			xsnp.push_back(xeg);
		}
		X.push_back(xsnp);
	}
	
	//start bootstrapping
	for (h = 0; h < nboot; h++) {
		
		for (i = 0; i < snps.size(); i++) {
			for (j = 0; j < group.size(); j++) {
				temp = randomSample(X[i][j], X[i][j].size());
				
				for (l = 0; l < group[j].index.size(); l++) {
					k = group[j].index[l];
					samples[i].EG[k] = temp[l];
				}
			}
		}
		
		tsamp = testStatistic(samples, sample, group, false);

		if (tsamp[0] <= tobs[0])
			tcountLinear++;
		if (tsamp[1] <= tobs[1])
			tcountQuadratic++;

		bootcount++;

	}

	return{ (tcountLinear+1)/(bootcount+1), (tcountQuadratic+1)/(bootcount+1) };
}


double RVSrare(MultiTestSet ts, int nboot) {
	size_t i, j;
	int n = ts.length();

	VectorXd var = sqrt(ts.getRobustVariance().array());
	MatrixXd diagS = var.asDiagonal();

	VectorXd YmHRG = ts.getYm(true);
	VectorXd YmLRG = ts.getYm(false);

	ts.getCovariateMatrix(0);
	ts.getCorrelationMatrix(0);

	return 0;
}

double runRareTest(std::vector<SNP> &snps, std::vector<Sample> &sample, std::vector<Group> &group,
	int nboot, bool rvs) {

	MultiTestSet ts(snps, sample, group);
	return RVSrare(ts, nboot);
}

