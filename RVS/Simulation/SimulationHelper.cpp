#include "../RVS.h"
#include "Simulation.h"

/*
Produces a vector of simulated minor allele frequencies (uniformly random between min and max).

@param nsnp Number of MAF to produce.
@param min Minimum value for MAF.
@param max Maximum value for MAF.

@return Vector of minor allele frequencies.
*/
VectorXd simulateMinorAlleleFrequency(int nsnp, double min, double max) {
	VectorXd maf(nsnp);

	for (int i = 0; i < nsnp; i++) 
		maf[i] = randomDouble(min, max);

	return maf;
}

double inline generateGenotype(double prob_y, double prob_x0, double prob_x1) {

	double rand = randomDouble(0, prob_y);

	if (rand > prob_x0)
		if (rand < prob_x0 + prob_x1)
			return 1;
		else
			return 2;
	else
		return 0;
}

/*
Produces a matrix of expected genotypes for the population using minor allele frequency and odds ratio.

@param npop Number of individuals in the population.
@param ncase Number of affected individuals in the population.
@param oddsRatio Odds ratio of being affected given harmful variant.
@param Vector of minor allele frequencies for each variant.

@return Matrix of population expected genotypes.
*/
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

/*
Produces a vector of case/control status for the population.

@param npop Number of individuals in the population.
@param ncase Number of affected individuals in the population.

@return Vector of population case/control status.
*/
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

/*
Produces a matrix of expected genotypes for samples.

@param X Matrix of expected genotypes from the population.
@param nsamp Number of individuals in the sample.
@param ncase Number of affected individuals in the sample.
@param ncase_pop Number of affected individuals in the population.

@return Matrix of sample expected genotypes.
*/
MatrixXd sampleX(MatrixXd X, int nsamp, int ncase, int ncase_pop) {
	VectorXi indices;
	int rows = X.rows();
	int cols = X.cols();

	//randomize cases
	MatrixXd Xcase = X.block(0, 0, ncase_pop, cols);
	indices = VectorXi::LinSpaced(ncase_pop, 0, ncase_pop);
	std::random_shuffle(indices.data(), indices.data() + Xcase.rows());
	MatrixXd xcase = indices.asPermutation() * Xcase;
	MatrixXd x1 = xcase.block(0, 0, ncase, xcase.cols());

	//randomize controls
	int ncont_pop = rows - ncase_pop;
	int ncont = nsamp - ncase;
	MatrixXd Xcont = X.block(ncase_pop, 0, ncont_pop, cols);
	indices = VectorXi::LinSpaced(ncont_pop, 0, ncont_pop);
	std::random_shuffle(indices.data(), indices.data() + Xcont.rows());
	MatrixXd xcont = indices.asPermutation() * Xcont;
	MatrixXd x0 = xcont.block(0, 0, ncont, xcont.cols());

	MatrixXd x(nsamp, cols);
	x << x1, x0;

	return x;
}

/*
Produces a vector of characters representing the called bases given a
true genotype, read depth and sequencing platform error rate.

@param trueGenotype The base corresponding to the true genotype.
@param error Vector holding probability of error for each read.
@param readDepth Number of reads.

@return Vector of called bases.
*/
std::vector<char> baseCall(std::vector<char> trueGenotype, VectorXd error, int readDepth) {
	
	if (error.size() != readDepth) {
		std::cout << "dimension error in seq_call.\n";
	}

	std::vector<char> bases;
	double e3;
	double e23;
	double e;
	double r;
	char trueBase;
	std::vector<char> errorBase;

	for (int i = 0; i < readDepth; i++) {
		trueBase = trueGenotype[randomInt(0, 1)];
		if (trueBase == 'C')
			errorBase = {'A', 'G', 'T'};
		else if (trueBase == 'T')
			errorBase = { 'A', 'G', 'C' };
		else if (trueBase == 'G')
			errorBase = { 'A', 'T', 'C' };
		else if (trueBase == 'A')
			errorBase = { 'G', 'T', 'C' };

		e3 = error[i] / 3.0;
		e23 = e3 * 2.0;
		e = e3 * 3.0;
		r = randomDouble(0, 1);

		if (r <= e3)
			bases.push_back(errorBase[0]);
		else if (r > e3 && r <= e23)
			bases.push_back(errorBase[1]);
		else if (r > e23 & r <= e)
			bases.push_back(errorBase[2]);
		else
			bases.push_back(trueBase);
	}

	return bases;
}

/*
Produces a vector of 3 probabilities as follows:
0: probability of calling every base in bases in a homozygote individual (2 major alleles).
1: probability of calling every base in bases in a heterzygote individual (1 minor allele).
2: probability of calling every base in bases in a homozygote individual (2 minor alleles).

@param bases Vector of called bases.
@param error Vector holding probability of error for each read.

@return Vector of 3 probabilities.
*/
VectorXd calculateLikelihood(std::vector<char> &bases, VectorXd &error) {
	
	VectorXd ll(3);
	ll[0] = 1;
	ll[1] = 1;
	ll[2] = 1;

	for (int i = 0; i < bases.size(); i++) {
		ll[0] = ll[0] * pSingle(bases[i], 'T', 'T', error[i]);
		ll[1] = ll[1] * pSingle(bases[i], 'C', 'T', error[i]);
		ll[2] = ll[2] * pSingle(bases[i], 'C', 'C', error[i]);
	}
	return ll;
}

/*
Produces a matrix where each column has 3 probabilities as follows:
0: probability of calling base in a homozygote individual (2 major alleles).
1: probability of calling base in a heterzygote individual (1 minor allele).
2: probability of calling base in a homozygote individual (2 minor alleles).

Each row of the vector is a base call cooresponding to a different read.

@param bases Vector of called bases.
@param error Vector holding probability of error for each read.

@return nx3 matrix of probabilities (n = read depth).
*/
MatrixXd calculateLikelihood2(std::vector<char> &bases, VectorXd &error) {
	MatrixXd ll(bases.size(), 3);

	for (int i = 0; i < bases.size(); i++) {
		std::vector<double> row;
		ll(i, 0) = pSingle(bases[i], 'T', 'T', error[i]);
		ll(i, 1) = pSingle(bases[i], 'C', 'T', error[i]);
		ll(i, 2) = pSingle(bases[i], 'C', 'C', error[i]);
	}
	return ll;
}

/*
Calculates the probability of calling a base given the true genotype and base call error rate.

@param base Called base.
@param true1 True base 1.
@param true2 True base 2 (same as true1 = homozygote, different = heterozygote).
@param error Probability of a base call error.

@return Probability of calling base.
*/
double pSingle(char base, char true1, char true2, double error) {
	
	double p1, p2;

	if (base == true1)
		p1 = 1 - error;
	else
		p1 = error / 3.0;
	
	if (base == true2)
		p2 = 1 - error;
	else
		p2 = error / 3.0;

	return 0.5 * p1 + 0.5 * p2;
}

/*
Uses EM algorithm to estimate the genotype frequencies in the sample.

@param  M A vector with genotype likelihoods (calculated by function calculateLikelihood).
@return A vector with probability of 0, 1 or 2 minor alleles.
*/
VectorXd calcEM(MatrixXd M) {

	double p = 0.15;
	double q = 0.15;
	double qn = 1;
	double pn = 0;
	double dn = 0;
	double d = 0;

	int k = 0;

	while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

		d = 1 - p - q;
		Vector3d v = { p, q, d };
		VectorXd pD = M * v;
		VectorXd Ep = M.col(0).array() * (p / pD.array()).array();
		VectorXd Eq = M.col(1).array() * (q / pD.array()).array();
		pn = p;
		qn = q;
		dn = 1 - q - p;
		p = Ep.sum() / Ep.rows();
		q = Eq.sum() / Eq.rows();

		k++;
		if (k == 1000)
			break;
	}

	VectorXd freq(3);
	freq[0] = std::max(0.0, p);
	freq[1] = std::max(0.0, q);
	freq[2] = std::max(0.0, 1 - p - q);

	return freq;
}

/*
Generates the expected probabilities of the genotypes E(G_ij | D_ij).

@param  M A vector with genotype likelihoods (calculated by function calculateLikelihood2).
@param  p Population frequency (calculated by function calcEM).

@return Vector of expected genotypes.
*/
VectorXd calculateExpectedGenotypes(std::vector<MatrixXd> M, VectorXd p){
		
	VectorXd m(3);
	Vector3d g = { 0, 1, 2 };
	VectorXd EG( M.size());
	double pm;

	for (int i = 0; i < M.size(); i++) {
		for (int j = 0; j < 3; j++) {
			double L = 1;

			for (int k = 0; k < M[i].rows(); k++)
				L = L * M[i](k, j);

			m[j] = L * p[j];			
		}

		pm = (m.array() / m.sum() * g.array()).sum();
		pm = std::max(0.0, pm);
		pm = std::min(2.0, pm);
		EG[i] = pm;
	}

	return EG;
}

std::map<int, int> simulationToReadGroup(std::map<int, SimulationGroup> &group) {
	std::map<int, int> readGroup;

	for (int i = 0; i < group.size(); i++)
		readGroup[i] = group[i].hrg;

	return readGroup;
}