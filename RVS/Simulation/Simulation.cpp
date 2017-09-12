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

/*
Simulates a dataset which can be used for an association test.

@param req Parameters required to set up dataset.
@param X Matrix of explanatory variable.
@param Y Vector of response variable.
@param G Vector of group ID.
@param readGroup Mapping of group ID to high or low read group.
@param P Vector with probability of 0, 1 or 2 minor alleles.

@return None.
@effect Fills Y, Z, G, P and readGroup with information calculated using SimulationRequest.
*/
void simulate(SimulationRequest req, MatrixXd &X, VectorXd &Y, VectorXd &G, std::map<int, int> &readGroup, MatrixXd &P) {

	printInfo("Setting up simulation parameters.");

    int npop = req.npop;
    double prevalence = req.prevalence;

	int ncase_pop = floor(npop * prevalence);
	int ncont_pop = npop - ncase_pop;

    int nsnp = req.nsnp;

    double me = req.me;
    double sde = req.sde;

    double oddsRatio = req.oddsRatio;
    double upperMAF = req.upperMAF;
    double lowerMAF = req.lowerMAF;


	//take in # groups, high/low status vector and number per group (must sum to nsamp) and mean, sd
	//--------------------------------------------------------
	
	int nsamp;
	int ncase = 0;
	int ncont = 0;

	for (SimulationRequestGroup srg : req.groups) {
		if (srg.isCase)
			ncase += srg.n;
		else
			ncont += srg.n;
	}

	nsamp = ncase + ncont;

	std::map<int, SimulationGroup> group;
	VectorXd g(nsamp);

	int groupIndex = 0;
	int gIndex = 0;

	//note: put cases first!!
	for (SimulationRequestGroup srg : req.groups) {
		if (srg.isCase) {
			group[groupIndex] = makeSimulationGroup(srg.n, srg.isHrg, srg.meanDepth, srg.sdDepth);

			for (int i = gIndex; i < srg.n + gIndex; i++)
				g[i] = groupIndex;

			gIndex += srg.n;
			groupIndex++;
		}
	}
	for (SimulationRequestGroup srg : req.groups) {
		if (!srg.isCase) {
			group[groupIndex] = makeSimulationGroup(srg.n, srg.isHrg, srg.meanDepth, srg.sdDepth);

			for (int i = gIndex; i < srg.n + gIndex; i++)
				g[i] = groupIndex;

			gIndex += srg.n;
			groupIndex++;
		}
	}


	//--------------------------------------------------------

    VectorXd maf = simulateMinorAlleleFrequency(nsnp, lowerMAF, upperMAF);
	//todo:?
	/* or we can determine the minor allele frequency fixed for each collapeed SNPs(5 SNPs in the current setting)
		njoint = 5
		loopno = nsnp / njoint

		mafco_5 = c(0.040856639, 0.021548447, 0.015053696, 0.005911941, 0.022559459)
		mafco = rep(mafco_5, loopno)
	*/

	printInfo("Simulating population data.");

	MatrixXd Xpop = simulatePopulationX(npop, ncase_pop, oddsRatio, maf);
	VectorXd Ypop = simulatePopulationY(npop, ncase_pop);

	printInfo("Simulating sample data.");

	MatrixXd x = sampleX(Xpop, nsamp, ncase, ncase_pop);
	VectorXd y = sampleY(Ypop, nsamp, ncase, ncase_pop);

	MatrixXd EG(nsamp, nsnp);
	MatrixXd p(nsnp, 3);

	int tenth = std::max(nsnp / 10, 1);

	for (int i = 0; i < EG.cols(); i++) {

		if (i % tenth == 0)
			printInfo("Simulating SNP " + std::to_string(i) + "/" + std::to_string(EG.cols()) + ".");

		std::vector<VectorXd> results = generateSeqData(x.col(i), y, g, group, me, sde);
		EG.col(i) = results[0];
		p.row(i) = results[1];
	}
	
	X = EG;
	Y = y;
	G = g;
	readGroup = simulationToReadGroup(group);
	P = p;
}


