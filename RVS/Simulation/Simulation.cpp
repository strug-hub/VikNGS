#include "stdafx.h"
#include "../RVS.h"
#include "Simulation.h"

std::vector<VectorXd> generateSeqData(VectorXd x, VectorXd y, VectorXd g, std::map<int, SimulationGroup> group, double me, double sde) {
	int i, j, k;
	double mean, sd;
	int rd;

	MatrixXd M(x.rows(), 3);
	std::vector<MatrixXd> mM;

	for (int i = 0; i < x.rows(); i++) {

		mean = group[g[i]].mean;
		sd = group[g[i]].sd;

		rd = std::round(randomNormal(mean, sd));
		rd = std::max(rd, 1);

		VectorXd error(rd);
		for (int j = 0; j < error.rows(); j++)
			error[j] = randomNormal(me, sde);

		k = x[i];

		std::vector<char> trueGenotype;

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

void simulate(MatrixXd &X, VectorXd &Y, VectorXd &G, std::map<int, int> &readGroup, MatrixXd &P) {

	std::cout << "Simulating population data\n";

	int npop = 10000; //The number of population
	double prevalence = 0.2; //A decimal between[0, 1], prevalence rate of the disease.
	int ncase_pop = floor(npop * prevalence);
	int ncont_pop = npop - ncase_pop;


	int nsnp = 101;  //Integer.The number of variants or bases.

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

	VectorXd maf = simulateMinorAlleleFrequency(nsnp, 0.1, 0.5);
	//todo:?
	/* or we can determine the minor allele frequency fixed for each collapeed SNPs(5 SNPs in the current setting)
		njoint = 5
		loopno = nsnp / njoint

		mafco_5 = c(0.040856639, 0.021548447, 0.015053696, 0.005911941, 0.022559459)
		mafco = rep(mafco_5, loopno)
	*/

	std::cout << "Simulating sample data\n";

	MatrixXd Xpop = simulatePopulationX(npop, ncase_pop, oddsRatio, maf);
	VectorXd Ypop = simulatePopulationY(npop, ncase_pop);
	MatrixXd x = sampleX(Xpop, nsamp, ncase, ncase_pop);
	VectorXd y = sampleY(Ypop, nsamp, ncase, ncase_pop);

	MatrixXd EG(nsamp, nsnp);
	MatrixXd p(nsnp, 3);

	for (int i = 0; i < EG.cols(); i++) {

		if (i % 100 == 0) {
			std::cout << "Simulating SNP ";
			std::cout << i;
			std::cout << "/";
			std::cout << EG.cols();
			std::cout << "\n";
		}

		std::vector<VectorXd> tmp = generateSeqData(x.col(i), y, g, group, me, sde);
		EG.col(i) = tmp[0];
		p.row(i) = tmp[1];

	}
	
	X = EG;
	Y = y;
	G = g;
	readGroup = simulationToReadGroup(group);
	P = p;
}


