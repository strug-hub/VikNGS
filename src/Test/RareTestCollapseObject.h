#include "../Math/MathHelper.h"
#include "RareTestObject.h"

class RareTestCollapseObject {
private:

    std::vector<VectorXd> y_original;
    std::vector<VectorXd> ycenter_original;
    MatrixXd Z;

    std::vector<VectorXd> y;
    std::vector<VectorXd> ycenter;
    std::vector<MatrixXd> z;

    bool covariates;
    bool normal;
    bool binomial;

    std::vector<RareTestObject> t;
    void setFamily(std::string fam);

    MatrixXd getVarianceNormal(bool rvs);
    MatrixXd getVarianceBinomial(bool rvs);
    void normalBootstrap();
    void binomialBootstrap();

public:
    RareTestCollapseObject(VectorXd Y, MatrixXd Z, std::vector<VectorXd> y, std::vector<MatrixXd> z,
                           std::vector<VectorXd> &ycenter, bool covariates, std::string family) {

        this->y = y;
        y_original = y;
        this->z = z;
        Z = concatenate(z);
        this->ycenter = ycenter;
        ycenter_original = ycenter;
        this->covariates = covariates;
        setFamily(family);
	}

    void addVariant(std::vector<VectorXd> &x, std::vector<int> &readDepth, VectorXd &p) {

      RareTestObject newTest(x, readDepth, p);
      t.push_back(newTest);

    }

    MatrixXd getRobustVar();
    MatrixXd getYmHigh();
    MatrixXd getYmLow();
    MatrixXd getYm();

    VectorXd getScore();

    inline MatrixXd getVariance(bool rvs){
        if (normal)
            return getVarianceNormal(rvs);
        else if (binomial)
            return getVarianceBinomial(rvs);

    }

    inline int size() { return t.size(); }

	//bootstrapping ------------------------------

    void bootstrap();

	//bootstrapping ------------------------------

};
