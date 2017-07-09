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

double inline generateGenotype(double prob_y, double prob_x0, double prob_x1);
MatrixXd simulateControlPopulation(int n, double preval, double maf, double oddsRatio);
MatrixXd simulateCasePopulation(int n, double preval, double maf, double oddsRatio);
MatrixXd chooseSample(MatrixXd population, int nsamp);
VectorXd baseCall(std::vector<char> trueGenotype, VectorXd error, int readDepth);
