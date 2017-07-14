#include "stdafx.h"
#include "../RVS.h"
#include "Simulation.h"

void simulate() {
	
	int npop = 10000; //The number of population
	double preval = 0.2; //A decimal between[0, 1], prevalence rate of the disease.
	double nsnp = 1000;  //Integer.The number of variants or bases.

	double me = 0.01; //The mean error rate of sequencing.
	double sde = 0.025;  //The standard deviation for the error rate.

	int nsamp = 2000;
	int ncase = 500;
	int ncont = nsamp - ncase;

	std::vector<SimulationGroup> groups;
	groups.push_back(makeSimulationGroup(200, true, 100, 10));
	groups.push_back(makeSimulationGroup(300, true, 80, 5));
	groups.push_back(makeSimulationGroup(1500, false, 4, 1));



	//Generate population data

	//We can randomly generate MAF for each SNP
		
	//mafco = runif(nsnp, min = 0.001, max = 0.05) # For common variant choose bigger parameters for uniform.e.g.min = 0.1, max = 0.5


	//or we can determine the minor allele frequency fixed for each collapeed SNPs(5 SNPs in the current setting)
	int njoint = 5;
	int loopno = nsnp / njoint;

//		mafco_5 = c(0.040856639, 0.021548447, 0.015053696, 0.005911941, 0.022559459)
	//	mafco = rep(mafco_5, loopno)

//		pop.data = list()

	double oddsRatio = 1.0;  //Under H0

	double maf = 0.040856639;
/*
	MatrixXd populationCase = simulateCasePopulation(npop, preval, maf, oddsRatio);
	MatrixXd populationControl = simulateControlPopulation(npop, preval, maf, oddsRatio);

	MatrixXd sampleCase = chooseSample(populationCase, ncase);
	MatrixXd sampleControl = chooseSample(populationControl, ncont);

	double readDepth;
	SimulationGroup g;

	for (size_t i = 0; i < groups.size(); i++) {
		g = groups[i];

		for (size_t j = 0; j < g.n; j++) {

			readDepth = round(randomNormal(g.mean, g.sd));
			if (readDepth < 1)
				readDepth = 1;

			VectorXd error(readDepth);

			for (size_t k = 0; k < readDepth; k++) 
				error[k] = randomNormal(me, sde);
		
			std::vector<char> genotype;
			int k = sampleCase() //todo: add group index to sample case/control vector???

			k = pheno.geno$x[i]
				gen_vect = c('TT', 'CT', 'CC')
				if (k == 0) { genotype = c('T', 'T') }
			if (k == 1) { genotype = c('C', 'T') }
			if (k == 2) { genotype = c('C', 'C') }

			
	}




		a = seq_call(genotype, error, ndepth = rd) // read from one individual
			v = c(v, a)
			erv = c(erv, error)
			rdv = c(rdv, rd)

			*/
	}

	




	/*			
	v = NULL # variant in the reads concatenated person by person
	erv = NULL # error rate for each reads concatenated person by person
	rdv = NULL # vector of read depth for each indvidual(sum(rdv) = length(v) = length(erv))



	#  write.table(v, paste0(ncase, 'case_', ncont, 'control_OR', OR, '_reads_for_simulation.txt', col.names = F, sep = ' ', quote = F, append = T))
	#  write.table(rdv, 'read_depth_for_reads.txt', col.names = F, sep = ' ', quote = F, append = T)
	# final = list(pheno = pheno.geno$y, geno<-pheno.geno$x, exp_geno = exp_geno$exp_geno, geno_freq = exp_geno$geno_freq)

	exp_geno<-calc_simu_exp(v, erv, rdv)
	## reads have different length for each simulation, so rbind will have warning.
	final = list(geno = pheno.geno$x, pheno = pheno.geno$y, rdepth_vector = rdv, exp_geno = exp_geno$exp_geno, pop_freq = exp_geno$geno_freq)
	return(final)
	}
	*/




	//while(true){}
//}

