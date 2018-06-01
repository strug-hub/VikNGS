#pragma once
#include "../Math/MathHelper.h"

class CommonTestObject {
private:
    double n;
	//NANs have been removed
	std::vector<VectorXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;

    //----------------------
	std::vector<int> readDepth;
    std::vector<VectorXd> mu;
	std::vector<VectorXd> ycenter;
	double robustVar;

    bool covariates;
    bool regular;
    bool normal;
    bool binomial;

    void filterNAN();
    void setFamily(std::string fam);

    //hat matrix logic, not used
    //--------------------------
    bool useHatMatrix = false;
    //--------------------------
    void constructHatMatrix(VectorXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::vector<VectorXd> &yhat);
    double covariateVarianceCalculation(bool rvs);
    double covariateCovarianceCalculation();
    MatrixXd H;
    std::vector<MatrixXd> h;
    std::vector<VectorXd> hdiag;
    VectorXd Yvar;
    std::vector<VectorXd> yvar;
    //--------------------------

public:
    CommonTestObject(VectorXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G,
		std::vector<VectorXd> &x, std::vector<VectorXd> &y, std::vector<MatrixXd> &z,
        std::vector<int> &readDepth, VectorXd &P, std::string family, bool useRegular) {

        covariates = true;
        regular = useRegular;
        setFamily(family);

		this->x = x;
		this->y = y;
		this->z = z;
		this->readDepth = readDepth;

        filterNAN();

        VectorXd beta = getBeta(X, Y, Z, family);
        this->mu = fitModel(beta, this->y, this->z, family);
        for (int i = 0; i < y.size(); i++)
            ycenter.push_back(this->y[i] - mu[i]);

        if(useHatMatrix)
            constructHatMatrix(X, Y, Z, G, mu);

        robustVar = calcRobustVar(P);
	}

	CommonTestObject(std::vector<VectorXd> &x, std::vector<VectorXd> &y,
        std::vector<int> &readDepth, VectorXd &P, std::string family, bool useRegular) {

        covariates = false;
        regular = useRegular;
        setFamily(family);

		this->x = x;
		this->y = y;
		this->readDepth = readDepth;

        filterNAN();

        double ybar = average(this->y);

        for (int i = 0; i < size(); i++){
            VectorXd temp = VectorXd::Constant(y[i].rows(), ybar);
            mu.push_back(temp);
            ycenter.push_back(this->y[i].array() - ybar);
        }

		robustVar = calcRobustVar(P);
	}

    inline double getScore() {
        double score = 0;
        for(int i = 0; i < size(); i++)
            score += (ycenter[i].array() * x[i].array()).sum();

        return score;
    }

    double getVarianceBinomial(bool rvs);
    double getVarianceNormal(bool rvs);

    inline double getVariance(bool rvs){
        if(binomial)
            return getVarianceBinomial(rvs);
        if(normal)
            return getVarianceNormal(rvs);

        throwError("COMMON_TEST_OBJECT", "Don't know how to compute variance during test, this error should not happen.");
    }

    inline int size() { return x.size(); }

	//bootstrapping ------------------------------
	std::vector<VectorXd> xboot;

    void bootstrap();
	inline double bootScore(int i) { return (ycenter[i].array() * xboot[i].array()).sum(); }
    double bootVariance(int i);
	//bootstrapping ------------------------------

};
