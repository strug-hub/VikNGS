#include "stdafx.h"
#include "../RVS.h"
#include "Simulation.h"

void simulate() {

	int npop = 10000; //The number of population
	double prevalence = 0.2; //A decimal between[0, 1], prevalence rate of the disease.
	int ncase_pop = floor(npop * prevalence);
	int ncont_pop = npop - ncase_pop;


	double nsnp = 1000;  //Integer.The number of variants or bases.

	double me = 0.01; //The mean error rate of sequencing.
	double sde = 0.025;  //The standard deviation for the error rate.

	int nsamp = 2000;
	int ncase = 500;
	int ncont = nsamp - ncase;

	double oddsRatio = 1.0;  //Under H0

	//todo: function to create groups
	//take in # groups, high/low status vector and number per group (must sum to nsamp) and mean, sd
	//--------------------------------------------------------

	//note: put cases first!!
	std::map<int, SimulationGroup> group;
	group[0] = makeSimulationGroup(200, true, 100, 10);
	group[1] = makeSimulationGroup(300, true, 80, 5); 
	group[2] = makeSimulationGroup(1500, false, 4, 1); 

	VectorXd g(nsamp);
	for (int i = 0; i < nsamp; i++) {
		if (i < 200)
			g[i] = 0;
		else if (i >= 200 && i < 500)
			g[i] = 1;
		else
			g[i] = 2;
	}
	//--------------------------------------------------------

	VectorXd maf = simulateMinorAlleleFrequency(nsnp, 0.001, 0.05);
	//todo:?
	/* or we can determine the minor allele frequency fixed for each collapeed SNPs(5 SNPs in the current setting)
		njoint = 5
		loopno = nsnp / njoint

		mafco_5 = c(0.040856639, 0.021548447, 0.015053696, 0.005911941, 0.022559459)
		mafco = rep(mafco_5, loopno)
	*/

	MatrixXd X = simulatePopulationX(npop, ncase_pop, oddsRatio, maf);
	VectorXd Y = simulatePopulationY(npop, ncase_pop);





/*

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

