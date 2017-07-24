#include "stdafx.h"
#include "../RVS.h"
#include "Simulation.h"

std::vector<VectorXd> generateSeqData(VectorXd x, VectorXd y, VectorXd g, std::map<int, SimulationGroup> group, double me, double sde) {
	int i, j, k;
	double mean, sd;
	int rd;
	double error;
	std::vector<char> trueGenotype;

	MatrixXd M(3, x.rows());
	std::vector<MatrixXd> mM;

	for (int i = 0; i < x.rows(); i++) {

		mean = group[i].mean;
		sd = group[i].sd;

		rd = std::round(randomNormal(mean, sd));
		rd = std::max(rd, 1);

		VectorXd error(rd);
		for (int j = 0; j <= error.rows(); i++)
			error[j] = randomNormal(me, sde);

		k = x[i];
		if (k == 0) {
			trueGenotype.push_back('T');
			trueGenotype.push_back('T');
		}
		if (k == 1) {
			trueGenotype.push_back('C');
			trueGenotype.push_back('T');
		}
		if (k == 2) {
			trueGenotype.push_back('C');
			trueGenotype.push_back('C');
		}

		std::vector<char> bases = baseCall(trueGenotype, error, rd);
		M.row(i) = calculateLikelihood(bases, error);
		mM.push_back(calculateLikelihood2(bases, error));
	}

	VectorXd p = calcEM(M);
	VectorXd EG = calculateExpectedGenotypes(mM, p);
	return{ EG, p };
}

void simulate() {

	int npop = 100; //The number of population
	double prevalence = 0.2; //A decimal between[0, 1], prevalence rate of the disease.
	int ncase_pop = floor(npop * prevalence);
	int ncont_pop = npop - ncase_pop;


	double nsnp = 10;  //Integer.The number of variants or bases.

	double me = 0.01; //The mean error rate of sequencing.
	double sde = 0.025;  //The standard deviation for the error rate.

	int nsamp = 20;
	int ncase = 5;
	int ncont = nsamp - ncase;

	double oddsRatio = 1.0;  //Under H0

	//todo: function to create groups
	//take in # groups, high/low status vector and number per group (must sum to nsamp) and mean, sd
	//--------------------------------------------------------

	//note: put cases first!!
	std::map<int, SimulationGroup> group;
	group[0] = makeSimulationGroup(2, true, 100, 10);
	group[1] = makeSimulationGroup(3, true, 80, 5);
	group[2] = makeSimulationGroup(15, false, 4, 1);

	VectorXd g(nsamp);
	for (int i = 0; i < nsamp; i++) {
		if (i < 2)
			g[i] = 0;
		else if (i >= 2 && i < 5)
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
	MatrixXd x = sampleX(X, nsamp, ncase, ncase_pop);
	MatrixXd y = sampleY(Y, nsamp, ncase, ncase_pop);


	for (int i = 0; i < x.cols(); i++) {
		std::cout << generateSeqData(x.col(i), y, g, group, me, sde)[0];
	}
}


