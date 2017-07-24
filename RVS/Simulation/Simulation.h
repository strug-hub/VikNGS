#pragma once

struct SimulationGroup {
	int n;
	bool hrg;
	double mean;
	double sd;
};

inline SimulationGroup makeSimulationGroup(int n, bool hrg, double mean, double sd) {
	SimulationGroup g;
	g.n = n;
	g.hrg = hrg;
	g.mean = mean;
	g.sd = sd;
	return g;
}

std::map<int, int> simulationToReadGroup(std::map<int, SimulationGroup> &group);

VectorXd simulateMinorAlleleFrequency(int nsnp, double min, double max);
VectorXd simulatePopulationY(int npop, int ncase);
MatrixXd simulatePopulationX(int npop, int ncase, double oddsRatio, VectorXd maf);
VectorXd sampleY(VectorXd Y, int nsamp, int ncase, int ncase_pop);
MatrixXd sampleX(MatrixXd X, int nsamp, int ncase, int ncase_pop);

std::vector<char> baseCall(std::vector<char> trueGenotype, VectorXd error, int readDepth);

double inline generateGenotype(double prob_y, double prob_x0, double prob_x1);
VectorXd calculateLikelihood(std::vector<char> &bases, VectorXd &error);
MatrixXd calculateLikelihood2(std::vector<char> &bases, VectorXd &error);
double pSingle(char base, char true1, char true2, double error);
VectorXd calcEM(MatrixXd M);
VectorXd calculateExpectedGenotypes(std::vector<MatrixXd> M, VectorXd p);