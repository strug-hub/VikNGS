#include "../Math/MathHelper.h"
#include "RareTestObject.h"

class RareTestCollapseObject {
private:

    std::vector<VectorXd> y_original;
    std::vector<VectorXd> ycenter_original;
    std::vector<MatrixXd> z_original;

    MatrixXd Z;
    VectorXd Y;

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
    RareTestCollapseObject(std::vector<VectorXd> y, std::vector<MatrixXd> z,
                           std::vector<VectorXd> &ycenter, bool covariates, std::string family) {

        this->covariates = covariates;
        setFamily(family);

        this->y = y;
        y_original = y;

        if(covariates){
            this->z = z;
            z_original = z;
            Z = concatenate(z);
        }

        Y = concatenate(y);

        this->ycenter = ycenter;
        ycenter_original = ycenter;
	}

    void addVariant(std::vector<VectorXd> &x, std::vector<int> &readDepth, VectorXd &p) {

        RareTestObject newObj(x, readDepth, p);
        t.push_back(newObj);
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
    void bootstrap();
};
