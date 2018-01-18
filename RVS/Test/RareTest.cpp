#include "../RVS.h"
#include "RareTestObject.h"

double getTestStatistic(std::vector<RareTestObject> &t, bool rvs, std::string test) {
	int nsnp = t.size();
	int i, j;

	MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

	VectorXd var(nsnp);
	VectorXd ym_hrd(nsnp);
	VectorXd ym_lrd(nsnp);

	for (i = 0; i < nsnp; i++) {
		var[i] = sqrt(t[i].getRobustVar());
		ym_hrd[i] = t[i].getYm(1);
		ym_lrd[i] = t[i].getYm(0);
	}

	MatrixXd diagVar = var.asDiagonal();
	MatrixXd diagYm_hrd = ym_hrd.asDiagonal();
	MatrixXd diagYm_lrd = ym_lrd.asDiagonal();

	for (i = 0; i < t[0].size(); i++) {
		VectorXd x_0 = t[0].getX(i);
		MatrixXd x(x_0.rows(), nsnp);
		x.col(0) = x_0;

		for (j = 1; j < nsnp; j++)
			x.col(j) = t[j].getX(i);

		if (t[0].isHRG(i)) {
			if (rvs) {
				MatrixXd var_hrd = diagVar.transpose() * correlation(x) * diagVar;
				diagS += diagYm_hrd * var_hrd * diagYm_hrd;
			}
			else
				diagS += diagYm_hrd * covariance(x) * diagYm_hrd;
		}
		else
			diagS += diagYm_lrd * covariance(x) * diagYm_lrd;
	}

	VectorXd score(nsnp);
	for (i = 0; i < nsnp; i++)
		score[i] = t[i].getScore();

	VectorXd e = diagS.eigenvalues().real();
	std::vector<double> eigenvals(e.data(), e.data() + e.size());

	double tstat = 0;

	//CAST -> linear test
	if(test == "cast")
		tstat = pnorm(score.sum() / sqrt(diagS.sum()));
	//C-alpha -> quadratic test
	else if(test == "calpha")
		tstat = qfc(eigenvals, score.array().pow(2).sum(), score.rows());

	return tstat;
}

double rareTest(std::vector<RareTestObject> &t, int nboot, bool rvs, std::string test) {
	int i, j, h;
	int nsnp = t.size();

	double tobs = getTestStatistic(t, rvs, test);

	//start bootstrapping!

	double bootCount = 0;
	double count = 0;
	double tsamp;

	for (h = 0; h < nboot; h++) {
		for (i = 0; i < nsnp; i++)
			t[i].bootstrap();

		tsamp = getTestStatistic(t, false, test);

		if (tsamp <= tobs)
			count++;

		bootCount++;
	}

	return{ (count + 1) / (bootCount + 1)};
}

//0.00153997      0.0362593
//first 5
//0.667967        0.940141
//next 5
std::vector<double> runRareTestWithCovariates(Request req, TestInput input) {

	MatrixXd X = input.X;
	VectorXd Y = input.Y;
	MatrixXd Z = input.Z;

	VectorXd G = input.G;
	std::map<int, int> readGroup = input.readGroup;
	MatrixXd P = input.P;

	//todo
	int collapseNumber = 5;
	std::string test = "cast";

	int i, j, h;
	std::vector<double> pvals;

	std::vector<MatrixXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;
	std::vector<int> rd;

	int ngroups = 1 + (int)G.maxCoeff();

	for (i = 0; i < ngroups; i++) {
		x.push_back(extractRows(X, G, i));
		y.push_back(extractRows(Y, G, i));
		z.push_back(extractRows(Z, G, i));
		rd.push_back(readGroup[i]);
	}

	std::vector<RareTestObject> t;

	for (h = 0; h < X.cols(); h+=collapseNumber) {
		for (i = 0; i < collapseNumber; i++) {
			std::vector<VectorXd> x_i;
			for (j = 0; j < ngroups; j++)
				x_i.push_back(x[j].col(i));

			VectorXd X_i = X.col(i);
			VectorXd P_i = P.row(i);
			RareTestObject t_i(X_i, Y, Z, x_i, y, z, rd, P_i);
			t.push_back(t_i);
		}

        pvals.push_back(rareTest(t, req.nboot, req.rvs, test));
	}

	return pvals;
}

std::vector<double> runRareTestWithoutCovariates(Request req, TestInput input) {

	MatrixXd X = input.X;
	VectorXd Y = input.Y;

	VectorXd G = input.G;
	std::map<int, int> readGroup = input.readGroup;
	MatrixXd P = input.P;

	//todo
	int collapseNumber = 5;
	std::string test = "cast";


	int i, j, h;
	std::vector<double> pvals;

	std::vector<MatrixXd> x;
	std::vector<VectorXd> y;
	std::vector<int> rd;

	int ngroups = 1 + (int)G.maxCoeff();

	for (i = 0; i < ngroups; i++) {
		x.push_back(extractRows(X, G, i));
		y.push_back(extractRows(Y, G, i));
		rd.push_back(readGroup[i]);
	}

	std::vector<RareTestObject> t;

	for (h = 0; h < X.cols(); h += collapseNumber) {

		std::cout << "Collapsing variants ";
		std::cout << h;
		std::cout << " through ";
		std::cout << h + collapseNumber;
		std::cout << "\n";

		for (i = 0; i < collapseNumber; i++) {
			std::vector<VectorXd> x_i;
			for (j = 0; j < ngroups; j++)
				x_i.push_back(x[j].col(i));

			VectorXd X_i = X.col(i);
			VectorXd P_i = P.row(i);
			RareTestObject t_i(x_i, y, rd, P_i);
			t.push_back(t_i);
		}

		pvals.push_back(rareTest(t, req.nboot, req.rvs, test));

	}
	return pvals;
}


std::vector<double> runRareTest(Request req, TestInput input) {
	if (input.hasCovariates())
		return runRareTestWithCovariates(req, input);
	else
		return runRareTestWithoutCovariates(req, input);
}