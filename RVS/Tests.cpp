#include "stdafx.h"
#include "RVS.h"
#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>
#include <vector>
#include <random>
#include <time.h>

/*
Calculates the robust variance of E(G | D). var(x) = E(x^2) - E(x)^2

@param p Genotype frequency from EM algorithm
@return Robust variance of E(G | D)
*/
inline double calcRobustVar(std::vector<double> p) {
	return (4 * p[2] + p[1]) - pow(2 * p[2] + p[1], 2);
}


void calcMeanVar(std::vector<bool> &IDmap, std::vector<SNP> &snps) {
	double var;
	double mean;
	double n;
	double controlvar;
	double controlmean;
	double ncontrol;
	double casevar;
	double casemean;
	double ncase;
	double eg;

	for (size_t j = 0; j < snps.size(); j++) {

		var = 0;
		mean = 0;
		n = 0;
		controlvar = 0;
		controlmean = 0;
		ncontrol = 0;
		casevar = 0;
		casemean = 0;
		ncase = 0;

		for (size_t i = 0; i < snps[j].EG.size(); i++) {
			eg = snps[j].EG[i];
			if (eg != NULL) {
				if (!IDmap[i]) {
					ncontrol++;
					controlmean += eg;
				}
				else {
					ncase++;
					casemean += eg;
				}
			}
		}
		
		mean = casemean + controlmean;
		n = ncase + ncontrol;
	
		mean /= n;
		controlmean /= ncontrol;
		casemean /= ncase;

		for (size_t i = 0; i < snps[j].EG.size(); i++) {
			eg = snps[j].EG[i];

			if (eg != NULL) {
				var += pow((eg - mean), 2);

				if (!IDmap[i])
					controlvar += pow((eg - controlmean), 2);
				else {
					casevar += pow((eg - casemean), 2);
				}
			}
		}

		var /= n - 1;
		controlvar /= ncontrol - 1;
		casevar /= ncase - 1;

		snps[j].var = var;
		snps[j].mean = mean;
		snps[j].n = n;
		snps[j].controlvar = controlvar;
		snps[j].controlmean = controlmean;
		snps[j].ncontrol = ncontrol;
		snps[j].casevar = casevar;
		snps[j].casemean = casemean;
		snps[j].ncase = ncase;
	}
	
	return;
}


/*
RVS using asymptotic distribution for score test statistic. This functions includes
two association tests which differ in the estimation of the variance of score. It uses
var_case(E(G | D)) when rvs = false and var_case(G) when rvs = true.

Reference: http://www.ncbi.nlm.nih.gov/pubmed/22570057
Reference: http://www.ncbi.nlm.nih.gov/pubmed/24733292
@param snps Vector of SNPs.
@param IDmap Vector with phenotypes (case/control).
@param rvs Indicates how to calculate variance
@return p-values for the SNPs
*/
std::vector<double> RVSasy(std::vector<SNP> &snps, std::vector<bool> &IDmap, bool rvs) {
	
	std::vector<double> pvals;
	SNP snp;
	double v;
	double p;

	for (size_t j = 0; j < snps.size(); j++) {
		snp = snps[j];

		if (rvs) {
			p = snp.ncase / snp.n;

			v = (1-p) * calcRobustVar(snp.p) + p * snp.controlvar;
			v = snp.var / v;
		}
		else
			v = 1;

		pvals.push_back(chiSquareOneDOF(scoreTest(IDmap, snp.EG) * v));
	}

	return pvals;
}


double bstrapHelp1(SNP &snp, std::vector<bool> &IDmap, int nboot, double tobs) {
	srand(time(NULL));

	std::vector<double> x0;
	std::vector<double> x1;
	std::vector<double> counter0;
	std::vector<double> counter1;

	for (size_t i = 0; i < snp.EG.size(); i++) {
		if (snp.EG[i] != NULL) {
			if (IDmap[i]) {
				x1.push_back(snp.EG[i] - snp.casemean);
				counter1.push_back(0);
			}
			else {
				x0.push_back(snp.EG[i] - snp.controlmean);
				counter0.push_back(0);
			}
		}
	}

	double bootmean0;
	double bootvar0;
	double bootmean1;
	double bootvar1;
	double statistic;
	double bootScoreCount = 1;
	int ncase = int(snp.ncase);
	int ncontrol = int(snp.ncontrol);

	for (int k = 0; k < nboot; k++) {
		bootmean0 = 0;
		bootvar0 = 0;
		bootmean1 = 0;
		bootvar1 = 0;
		//reset counters
		for (int i = 0; i < ncontrol; i++)
			counter0[i] = 0;
		for (int i = 0; i < ncase; i++)
			counter1[i] = 0;

		//sample from observed values
		for (int i = 0; i < ncontrol; i++)
			counter0[rand() % ncontrol]++;
		for (int i = 0; i < ncase; i++)
			counter1[rand() % ncase]++;

		//calculate means
		for (int i = 0; i < ncontrol; i++)
			bootmean0 += counter0[i] * x0[i];
		for (int i = 0; i < ncase; i++)
			bootmean1 += counter1[i] * x1[i];

		bootmean0 /= snp.ncontrol;
		bootmean1 /= snp.ncase;

		//calculate variances
		for (int i = 0; i < ncontrol; i++)
			bootvar0 += counter0[i] * pow(x0[i] - bootmean0, 2);
		for (int i = 0; i < ncase; i++)
			bootvar1 += counter1[i] * pow(x1[i] - bootmean1, 2);

		bootvar0 /= snp.ncontrol - 1;
		bootvar1 /= snp.ncase - 1;

		statistic = (bootmean1 - bootmean0) / sqrt(bootvar1 / snp.ncase + bootvar0 / snp.ncontrol);
		if (abs(statistic) >= tobs)
			bootScoreCount++;
	}

	return 	bootScoreCount / nboot;
}

double bstrapHelp2(SNP &snp, std::vector<bool> &IDmap, int nboot, double tobs) {
	srand(time(NULL));

	std::vector<double> x;
	double bootScoreCount = 1;
	double p = snp.ncase / snp.n;
	double q = 1 - p;
	double vs = p * q * snp.n * snp.var;
	double statistic;
	double casesum;
	double controlsum;

	double xindex;

	for (size_t i = 0; i < snp.EG.size(); i++)
		if (snp.EG[i] != NULL)
			x.push_back(i);

	for (int k = 0; k < nboot; k++) {
		casesum = 0;
		controlsum = 0;
		xindex = 0;
		std::random_shuffle(x.begin(), x.end());

		for (size_t i = 0; i < snp.EG.size(); i++) {
			if (snp.EG[i] != NULL) {
				if (IDmap[i])
					casesum += snp.EG[x[xindex]];
				else
					controlsum += snp.EG[x[xindex]];
				xindex++;
			}
		}

		statistic = pow(q * casesum - p * controlsum, 2) / vs;

		if (abs(statistic) >= tobs)
			bootScoreCount++;
	}
	return 	bootScoreCount / nboot;
}

double RVSbtrap(std::vector<SNP> &snps, std::vector<bool> &IDmap, int nboot, bool rvs) {
	SNP snp = snps[0];

	double maf = 0;
	double p;
	double s;
	double tobs;

	//calculate observed score statistic
	if (rvs)
		tobs = (snp.casemean - snp.controlmean) / 
		sqrt(calcRobustVar(snp.p) / snp.ncase + snp.controlvar / snp.ncontrol);
	else {
		p = snp.ncase / snp.n;
		s = (1 - p) * snp.casemean * snp.ncase - p * snp.controlmean * snp.ncontrol;
		tobs = pow(s, 2) / (p * (1 - p) * snp.n * snp.var);
	}
	tobs = abs(tobs);

	//bootstrap test
	if (rvs)
		return bstrapHelp1(snp, IDmap, nboot, tobs);
	else
		return bstrapHelp2(snp, IDmap, nboot, tobs);
}


//TODO
//takes matrix
//finds columns which are all identical (except NULLs)
//picks first row without *any* NULL or random row  without *any* NULL (depends on hom)
//in columns with identical rows, takes above^^ row and if it is 0, sets to 10e-15
//if not zero, either divides cell by 2 or multiplies by 2 and takes min of 2 and itself
//eg.
//[, 1][, 2][, 3]			[, 1][, 2][, 3]
//[1, ]	1	1	0			[1, ]	NA	1	0.00E+00
//[2, ]	1	1	0			[2, ]	2	1	1.00E-15
//[3, ]	1	2	0			[3, ]	1	2	0.00E+00
//[4, ]	1	3	0	->		[4, ]	1	3	0.00E+00
//[5, ]	1	1	0			[5, ]	1	1	0.00E+00
//[6, ]	1	1	0			[6, ]	1	1	0.00E+00
//[7, ]	1	4	0			[7, ]	1	4	0.00E+00
//[8, ]	1	5	0			[8, ]	1	5	0.00E+00
//[9, ]	1	3	0			[9, ]	1	3	0.00E+00
//[10,]	1	2	0			[10, ]	1	2	0.00E+00
//void checkHomMatrix()
	
double pnorm(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p*x);
	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

	return 0.5*(1.0 + sign*y);
}



void RVSrare(std::vector<SNP> &snps, std::vector<bool> &IDmap, int nboot, bool rvs, int njoint, int hom, int multiplier) {
	SNP snp = snps[0];
	double temp;
	double s = 0;
	double s2 = 0;
	typedef std::vector<std::vector<double>> Matrix;
	typedef std::vector<double> Row;
	MatrixXd sigma(njoint, njoint);

	double sum;
	double SLobs;
	double SQobs;

	//calculate once
	double totalncase = 0;
	double totalncontrol = 0;
	double meany = 0;
	double vary = 0;

	for (size_t k = 0; k < IDmap.size(); k++)
		if (IDmap[k])
			totalncase++;
	totalncontrol = IDmap.size() - totalncase;
	meany = totalncase / IDmap.size();

	vary += totalncase * pow(1 - meany, 2);
	vary += totalncontrol * pow(meany, 2);
	vary /= (totalncase + totalncontrol - 1);

	double p = totalncontrol / totalncase;
	double q = 1 / p;

	//========== for loop here

	//calculate score vector
	for (int h = 0; h < njoint; h++) {
		snp = snps[h];
		temp = snp.ncase / snp.n;
		temp = snp.casemean * snp.ncase * (1 - temp) - temp * snp.controlmean * snp.ncontrol;
		s += temp;
		s2 += temp*temp;
	}

	//TODO: check_hom_matrix

	if (rvs) {
		double n;
		double vari;
		double varj;
		double meani;
		double meanj;
		std::vector<double> v;

		//calculate correlation matrix
		for (size_t i = 0; i < njoint; i++) {
			for (size_t j = 0; j < njoint; j++) {
				if (i == j) {
					sigma(i, j) = 1;
					sigma(j, i) = 1;
					v.push_back(sqrt(calcRobustVar(snps[i].p)));
				}
				else if (i < j) {
					sum = 0;
					n = 0;
					vari = 0;
					varj = 0;
					meani = 0;
					meanj = 0;

					for (size_t k = 0; k < IDmap.size(); k++) {
						if (IDmap[k] && snps[i].EG[k] != NULL && snps[j].EG[k] != NULL) {
							n++;
							meani += snps[i].EG[k];
							meanj += snps[j].EG[k];
						}
					}

					meani /= n;
					meanj /= n;

					for (size_t k = 0; k < IDmap.size(); k++) {
						if (IDmap[k] && snps[i].EG[k] != NULL && snps[j].EG[k] != NULL) {
							vari += pow((snps[i].EG[k] - meani), 2);
							varj += pow((snps[j].EG[k] - meanj), 2);
							sum += (snps[i].EG[k] - meani) * (snps[j].EG[k] - meanj);
						}
					}

					vari /= n - 1;
					varj /= n - 1;
					sum /= (n - 1) * sqrt(vari * varj);
					sigma(i, j) = sum;
					sigma(j, i) = sum;
				}
			}
		}

		//TODO: check for NULLs in sigma?
		//TODO: check_hom_matrix

		for (size_t i = 0; i < njoint; i++){
			for (size_t j = 0; j < njoint; j++) {
				if (i <= j) {
					sigma(i, j) *= v[i] * v[j] * totalncase * p;
					sigma(j, i) = sigma(i, j);
				}
			}
		}

		for (size_t i = 0; i < njoint; i++) {
			for (size_t j = 0; j < njoint; j++) {
				if (i == j) {
					sum = 0;
					for (size_t k = 0; k < IDmap.size(); k++)
						if (!IDmap[k] && snps[i].EG[k] != NULL)
							sum += pow(snps[i].EG[k] - snps[i].controlmean, 2);

					sigma(i, i) += q * sum;
					sigma(i, i) *= vary;
				}
				else if (i < j) {
					sum = 0;
					for (size_t k = 0; k < IDmap.size(); k++)
						if (!IDmap[k] && snps[i].EG[k] != NULL && snps[j].EG[k] != NULL)
							sum += (snps[i].EG[k] - snps[i].controlmean) * (snps[j].EG[k] - snps[j].controlmean);

					sigma(i, j) += q * sum;
					sigma(i, j) *= vary;
					sigma(j, i) = sigma(i, j);
				}
			}
		}
	}
	else {
		//TODO: check_hom_matrix?

		double sumcase;
		double sumcontrol;
		double ncase;
		double ncontrol;
		double meancontroli;
		double meancontrolj;
		double meancasei;
		double meancasej;

		std::vector<double> vs;

		//calculate covariance matrix		
		for (size_t i = 0; i < njoint; i++) {
			for (size_t j = 0; j < njoint; j++) {
				if (i <= j) {
					ncase = 0;
					ncontrol = 0;
					sumcase = 0;
					sumcontrol = 0;
					meancasei = 0;
					meancasej = 0;
					meancontroli = 0;
					meancontrolj = 0;

					for (size_t k = 0; k < IDmap.size(); k++) {
						if (snps[i].EG[k] != NULL && snps[j].EG[k] != NULL) {
							if (IDmap[k]) {
								ncase++;
								meancasei += snps[i].EG[k];
								meancasej += snps[j].EG[k];
							}
							else {
								ncontrol++;
								meancontroli += snps[i].EG[k];
								meancontrolj += snps[j].EG[k];
							}
						}
					}

					meancasei /= ncase;
					meancasej /= ncase;
					meancontroli /= ncontrol;
					meancontrolj /= ncontrol;

					for (size_t k = 0; k < IDmap.size(); k++) {
						if (snps[i].EG[k] != NULL && snps[j].EG[k] != NULL) {
							if (IDmap[k])
								sumcase += (snps[i].EG[k] - meancasei) * (snps[j].EG[k] - meancasej);
							else
								sumcontrol += (snps[i].EG[k] - meancontroli) * (snps[j].EG[k] - meancontrolj);
						}
					}

					sumcase /= ncase - 1;
					sumcontrol /= ncontrol - 1;
					sumcase = totalncontrol / totalncase * sumcase * (totalncase - 1) +
						totalncase / totalncontrol * sumcontrol * (totalncontrol - 1);

					sigma(i, j) = sumcase;
					sigma(j, i) = sumcase;

					if (i == j) 
						vs.push_back(sqrt(ncase*ncontrol) / (ncase + ncontrol));
				}
			}
		}

		for (size_t i = 0; i < njoint; i++) 
			for (size_t j = 0; j < njoint; j++) 
				sigma(i, j) *= vs[i] * vs[j];
			
	}

	sum = 0;
	for (size_t i = 0; i < njoint; i++)
		for (size_t j = 0; j < njoint; j++)
			sum += sigma(i, j);

	SLobs = 2 * pnorm(s / sqrt(sum));
	//SQobs


	VectorXd e = sigma.eigenvalues().real();
	std::vector<double> eigenvals(e.data(), e.data() + e.size());
	
	std::cout << qfc(eigenvals, s2, njoint);
}
