#include "stdafx.h"
#include "RVS.h"

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

#include "Eigen/Dense"
using Eigen::MatrixXd;

void RVSrare(std::vector<SNP> &snps, std::vector<bool> &IDmap, int nboot, bool rvs, int njoint, int hom, int multiplier) {
	SNP snp = snps[0];
	double temp;
	double s = 0;
	double s2 = 0;
	typedef std::vector<std::vector<double>> Matrix;
	typedef std::vector<double> Row;
	MatrixXd sigma(njoint, njoint);
	std::cout << sigma << std::endl;

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
	std::cout << SLobs;

	std::cout << '\n';

	for (size_t i = 0; i < njoint; i++) {
		for (size_t j = 0; j < njoint; j++) {
			std::cout << sigma(i, j);
			std::cout << '\t';

		}
		std::cout << '\n';

	}

	std::cout << sigma.eigenvalues();

}


/*
# RVS analysis with rare variants
# Robust Variance Estimate
# Get p-values from RVS with CAST and C-alpha (Resampling: Bootstrap, Variance Estimate: Robust)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values, it has dimensions n by J, J - number of variants, n - sample size
#' @param nboot - number of bootstrap
#' @param njoint - number of SNPs grouping together for one test, default is 5
#' @param P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic.
#' @return two p-values for CAST and C-alpha
#' @export
RVS_rare= function(Y,X,P,njoint=5,nboot,RVS='TRUE',hom=1,multiplier=1){
if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
cat('Wrong input for option RVS in RVS_rare, should be True or False.\n')
return(NULL)
}

nsnp=nrow(P)
if(nsnp%%njoint!=0) {nloop=(nsnp-nsnp%%njoint)/njoint ## the last couple of snps(<5) combined to the last group
}else {nloop=nsnp/njoint}

p.RVS.rare=NULL
for(i in 1:nloop){
t1<-(i-1)*njoint+1;
if(i < nloop) {t2<-i*njoint;
}else{t2 <- nsnp}
Xsub=X[,t1:t2]

tt=RVS_rare1(Y,X[,t1:t2],P[t1:t2,],njoint,nboot,RVS,hom=hom,multiplier=multiplier,snp_loop=i)
p.RVS.rare<-rbind(p.RVS.rare,tt)
}

row.names(p.RVS.rare)=paste0('loop',1:nloop)
p.RVS.rare=data.frame(p.RVS.rare)
colnames(p.RVS.rare)=c('p.CAST','p.Calpha')
return (p.RVS.rare)
}



#' Function to calculate the score test for given expected genotype probability for case and controls seperately (M1, M2)
#'
#' This function calcuates the score test \eqn{T=S^2/var(S)}, for a given \eqn{j, S_j=\sum_i (Y_i-bar(Y))E(G_{ij}|D_{ij})=\sum_{case}(1-bar(Y))E(G|D)-\sum_{cont}bar(Y)E(G|D)}. It is called in \code{\link{regScore_perm}}.
#' @param M1   a vector of the expected  genotype probability for case \eqn{E(G_{ij}|D_{ij})} on one snp
#' @param M2   a vector of the expected genotype probability for control \eqn{E(G_{ij}|D_{ij})} on one snp
#' @return the score test statistic calculated from M1 and M2
calc_score_test = function(M1,M2){
X = c(M1,M2)
Y = c(rep(1,length(M1)),rep(0,length(M2)))
p = length(M1)/length(X)
q = 1 - p
S = q*sum(M1)-p*sum(M2)
vs = p*q*(length(X))*var(X)
x = S^2/vs
return (x)
}










#' use RVS to test associaton by bootrtrap, given phenotype (Y), expected values of genotypes for case and controls (X)
#' estimated genotype frequency (P) and number of times of bootstrap (nboot)
#' @param Y - a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X - a vector of genotype for a snp, first ncase and then ncont
#' @param P (vector with length 3  (double)) - estimated genotype frequency for a variant P(G=0), P(G=1) and P(G=2)
#' @param nboot - number of bootstrap
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return   p-values for the snp
#' @export
RVS_btrap = function(Y, X, P, nboot, RVS = 'TRUE') {
	if (!RVS %in% c('TRUE', 'True', 'true', 'T', 't', 'FALSE', 'False', 'false', 'F', 'f')) {
		cat('Wrong input for option RVS in RVS_brap, should be True or False.\n')
			return(NULL)
	}
	Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		maf = sum(X) / (length(X) * 2); maf0 = min(maf, 1 - maf)
		if (maf0<0.05) { cat('Warning: MAF of the SNP is', maf0, ', try to group with other rare variants and use RVS_rare1.\n') }

	ncase1 = sum(Y == 1)
		ncont1 = sum(Y == 0)
		p = length(X[Y == 1]) / length(X)
		q = 1 - p
#  if(RVS %in% c('TRUE','True','true','T','t')) {
		vcase = as.numeric(calc_robust_Var(P))
#  } else{
		#    vcase = var(X[Y == 1])
#  }
		vcont = var(X[Y == 0])
		if (RVS %in% c('TRUE', 'True', 'true', 'T', 't')) {
			Tobs = (mean(X[Y == 1]) - mean(X[Y == 0])) / sqrt(vcase / ncase1 + vcont / ncont1)
		}
		else {
			X1 = X[Y == 1]; X2 = X[Y == 0]
				Tobs = calc_score_test(X1, X2)
		}
		X1 = X[Y == 1] - mean(X[Y == 1])
			X2 = X[Y == 0] - mean(X[Y == 0])
			C = NULL
			if (RVS %in% c('TRUE', 'True', 'true', 'T', 't')) {
				for (j in 1 : nboot) {
					Xca = sample(X1[], ncase1, replace = TRUE)
						Xco = sample(X2[], ncont1, replace = TRUE)
						vcase = var(Xca)
						vcont = var(Xco)
						C = c(C, (mean(Xca) - mean(Xco)) / sqrt(vcase / ncase1 + vcont / ncont1))
				}## enf for j
			}
			else {
				for (j in 1 : nboot) {
					k = sample(length(X));
					Xnew = X[k]
						X1 = Xnew[Y == 1]; X2 = Xnew[Y == 0]
						C = c(C, calc_score_test(X1, X2))
				} ## end for j
			}## else if RVS
				cc = (sum(abs(C) >= abs(Tobs)) + 1) / (nboot + 1)
						return(cc)
}
*/