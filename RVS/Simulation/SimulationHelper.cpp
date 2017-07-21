#include "stdafx.h"
#include "../RVS.h"
#include "Simulation.h"

VectorXd simulateMinorAlleleFrequency(int nsnp, double min, double max) {
	VectorXd maf(nsnp);
	for (int i = 0; i < nsnp; i++) 
		maf[i] = generateRandomDouble(min, max);

	return maf;
}

double inline generateGenotype(double prob_y, double prob_x0, double prob_x1) {

	double rand = generateRandomDouble(0, prob_y);

	if (rand > prob_x0)
		if (rand < prob_x0 + prob_x1)
			return 1;
		else
			return 2;
	else
		return 0;
}

MatrixXd simulatePopulationX(int npop, int ncase, double oddsRatio, VectorXd maf) {
	int nsnp = maf.size();
	MatrixXd X(npop, nsnp);
	
	int ncont = npop - ncase;

	double p_y0 = ncont / (double)npop;
	double p_y1 = 1 - p_y0;

	for (int i = 0; i < nsnp; i++) {

		double beta0 = ncase / (ncont * pow(maf[i] * (1 - oddsRatio) - 1, 2));

		double odds_y1_x0 = beta0;
		double odds_y1_x1 = beta0 * oddsRatio;
		double odds_y1_x2 = beta0 * oddsRatio * oddsRatio;


		double p_x0_y0 = ncont / (double) npop * pow(1 - maf[i], 2);
		double p_x1_y0 = 2 * maf[i] * (1 - maf[i]) * ncont / (double) npop;
		double p_x0_y1 = odds_y1_x0 * p_x0_y0;
		double p_x1_y1 = odds_y1_x0 * oddsRatio * p_x1_y0;

		for (int j = 0; j < npop; j++) 
			X(j, i) = generateGenotype(p_y1, p_x0_y1, p_x1_y1);
	}

	return X;
}

VectorXd simulatePopulationY(int npop, int ncase) {
	VectorXd Y(npop);

	for (int i = 0; i < npop; i++) {
		if (i < ncase)
			Y[i] = 1;
		else
			Y[i] = 0;
	}

	return Y;
}


VectorXd sampleY(VectorXd Y, int nsamp, int ncase, int ncase_pop) {
	return simulatePopulationY(nsamp, ncase);
}

MatrixXd sampleX(MatrixXd X, int nsamp, int ncase, int ncase_pop) {
	MatrixXd x(nsamp, X.cols());

	//sample without replacement
	//use std::shuffle

	return x;
}


void generateSeqData(MatrixXd x, VectorXd y, VectorXd g) {



}




/*
generate_seqdata_OR1_pop<-function(pop.data, ncase, ncont, rd.group, nhrd1, nhrd2, nlrd, mdhrd1, sdhrd1, mdhrd2, sdhrd2, mdlrd, sdlrd, me, sde) {
	Ntotal = ncase + ncont
		pheno.geno = choose.sample(pop.data, ncase, ncont)
		v = NULL # variant in the reads concatenated person by person
		erv = NULL # error rate for each reads concatenated person by person
		rdv = NULL # vector of read depth for each indvidual(sum(rdv) = length(v) = length(erv))

		for (i in 1 : Ntotal)
		{
			if (rd.group[i] == 0) {
				rd = round(sdhrd1*rnorm(1) + mdhrd1)
			}
			else if (rd.group[i] == 1) {
				rd = round(sdhrd2*rnorm(1) + mdhrd2)
			}
			else {
				rd = round(sdlrd*rnorm(1) + mdlrd)
			}
			if (rd <= 0) { rd = 1 }
			error = sde*rnorm(rd) + rep(me, rd)
				k = pheno.geno$x[i]
				gen_vect = c('TT', 'CT', 'CC')
				if (k == 0) { genotype = c('T', 'T') }
			if (k == 1) { genotype = c('C', 'T') }
			if (k == 2) { genotype = c('C', 'C') }

			a = seq_call(genotype, error, ndepth = rd) ## read from one individual
				v = c(v, a)
				erv = c(erv, error)
				rdv = c(rdv, rd)
		}

	#  write.table(v, paste0(ncase, 'case_', ncont, 'control_OR', OR, '_reads_for_simulation.txt', col.names = F, sep = ' ', quote = F, append = T))
		#  write.table(rdv, 'read_depth_for_reads.txt', col.names = F, sep = ' ', quote = F, append = T)
		# final = list(pheno = pheno.geno$y, geno<-pheno.geno$x, exp_geno = exp_geno$exp_geno, geno_freq = exp_geno$geno_freq)

		exp_geno<-calc_simu_exp(v, erv, rdv)
## reads have different length for each simulation, so rbind will have warning.
		final = list(geno = pheno.geno$x, pheno = pheno.geno$y, rdepth_vector = rdv, exp_geno = exp_geno$exp_geno, pop_freq = exp_geno$geno_freq)
		return(final)
}

*/



VectorXd baseCall(std::vector<char> trueGenotype, VectorXd error, int readDepth) {
	
	if (error.size() != readDepth) {
		std::cout << "dimension error in seq_call.\n";
	}

	VectorXd bases(readDepth);
	double e3;
	double e23;
	double e;
	double r;
	char trueBase;
	std::vector<char> errorBase;

	for (size_t i = 0; i < readDepth; i++) {
		trueBase = trueGenotype[generateRandomInteger(0, 1)];
		if (trueBase == 'C')
			errorBase = {'A', 'G', 'T'};
		else if (trueBase == 'T')
			errorBase = { 'A', 'G', 'C' };
		else if (trueBase == 'G')
			errorBase = { 'A', 'T', 'C' };
		else if (trueBase == 'A')
			errorBase = { 'G', 'T', 'C' };

		e3 = error[i] / 3;
		e23 = e3 * 2;
		e = e3 * 3;
		r = generateRandomDouble(0, 1);

		if (r <= e3)
			bases[i] = errorBase[0];
		else if (r > e3 && r <= e23)
			bases[i] = errorBase[1];
		else if (r > e23 & r <= e)
			bases[i] = errorBase[2];
		else
			bases[i] = trueBase;
	}

	return bases;
}


