#include "stdafx.h"
#include "RVS.h"
#include "Simulation.h"

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

MatrixXd simulateCasePopulation(int n, double preval, double maf, double oddsRatio) {

	int ncont = n - floor(n*preval);
	int ncase = n - ncont;

	double prob_y0 = ncont / (n * 1.0);
	double prob_y1 = 1 - prob_y0;

	double odds_y1_x0 = (n - ncont) / (ncont* pow(maf*(1 - oddsRatio) - 1, 2));

	double prob_x0_y0 = ncont / (1.0*n)*pow(1 - maf, 2);
	double prob_x1_y0 = 2 * maf*(1 - maf) * ncont / (1.0*n);
	double prob_x0_y1 = odds_y1_x0 * prob_x0_y0;
	double prob_x1_y1 = odds_y1_x0 * oddsRatio * prob_x1_y0;

	MatrixXd cases(n, 2);

	for (size_t i = 0; i < ncase; i++) {
		cases(i, 0) = 1;
		cases(i, 1) = generateGenotype(prob_y1, prob_x0_y1, prob_x1_y1);
	}

	return cases;
}

MatrixXd simulateControlPopulation(int n, double preval, double maf, double oddsRatio) {

	int ncont = n - floor(n*preval);

	double prob_y0 = ncont / (n * 1.0);
	double prob_x0_y0 = ncont / (1.0*n)*pow(1 - maf, 2);
	double prob_x1_y0 = 2 * maf*(1 - maf) * ncont / (1.0*n);

	MatrixXd controls(n, 2);

	for (size_t i = 0; i < ncont; i++) {
		controls(i, 0) = 0;
		controls(i, 1) = generateGenotype(prob_y0, prob_x0_y0, prob_x1_y0);
	}

	return controls;
}


MatrixXd chooseSample(MatrixXd population, int nsamp) {
	MatrixXd sample(nsamp, 2);
	
	int temp = 0;
	for (size_t i = 0; i < nsamp; i++) {
		temp = generateRandomInteger(0, population.rows());
		sample(i, 0) = population(temp, 0);
		sample(i, 1) = population(temp, 1);
	}

	return sample;
}