#pragma once
#include "../Math/MathHelper.h"

class CommonTestObject {
private:
	//NANs have been removed
	std::vector<VectorXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;

    MatrixXd H;
    std::vector<MatrixXd> h;
    std::vector<VectorXd> hdiag;
    VectorXd Yvar;
    std::vector<VectorXd> yvar;

    //----------------------
	std::vector<int> readDepth;
	std::vector<VectorXd> ycenter;
	double robustVar;
    bool covariates;

    bool normal;
    bool binomial;

    void filterNAN();
    void constructHatMatrix(VectorXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::vector<VectorXd> &yhat);
    void setFamily(std::string fam);

public:
    CommonTestObject(VectorXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G,
		std::vector<VectorXd> &x, std::vector<VectorXd> &y, std::vector<MatrixXd> &z,
        std::vector<int> &readDepth, VectorXd &P, std::string family) {

        covariates = true;
        setFamily(family);

		this->x = x;
		this->y = y;
		this->z = z;
		this->readDepth = readDepth;

        filterNAN();

        VectorXd beta = getBeta(X, Y, Z, family);
        std::vector<VectorXd> yhat = fitModel(beta, this->y, this->z, family);
        for (int i = 0; i < y.size(); i++)
            ycenter.push_back(this->y[i] - yhat[i]);

        constructHatMatrix(X, Y, Z, G, yhat);

        robustVar = calcRobustVar(P);
	}

	CommonTestObject(std::vector<VectorXd> &x, std::vector<VectorXd> &y,
        std::vector<int> &readDepth, VectorXd &P, std::string family) {

        covariates = false;
        setFamily(family);

		this->x = x;
		this->y = y;
		this->readDepth = readDepth;


        filterNAN();

		double ybar = average(this->y);

        for (int i = 0; i < size(); i++)
			ycenter.push_back(this->y[i].array() - ybar);

		robustVar = calcRobustVar(P);
	}

    inline double getScore() {
        double score = 0;
        for(int i = 0; i < size(); i++)
            score += (ycenter[i].array() * x[i].array()).sum();

        return score;
    }

    double getVariance(bool rvs);
    double covariateVarianceCalculation(bool rvs);
    double covariateCovarianceCalculation();
    inline int size() { return x.size(); }

	//bootstrapping ------------------------------
	std::vector<VectorXd> xboot;

    void bootstrap();
	inline double bootScore(int i) { return (ycenter[i].array() * xboot[i].array()).sum(); }
    double bootVariance(int i);
	//bootstrapping ------------------------------

};
