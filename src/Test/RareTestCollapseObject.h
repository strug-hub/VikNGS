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
    std::vector<VectorXd> mu;
    std::vector<VectorXd> ycenter;
    std::vector<MatrixXd> z;

    bool covariates;
    bool normal;
    bool binomial;
    bool regular;

    std::vector<RareTestObject> t;
    void setFamily(std::string fam);

    MatrixXd getVarianceRegular();
    MatrixXd getVarianceNormal(bool rvs);
    MatrixXd getVarianceBinomial(bool rvs);
    void normalBootstrap();
    void binomialBootstrap();

public:
    RareTestCollapseObject(std::vector<VectorXd> &y, std::vector<MatrixXd> &z,
                           std::vector<VectorXd> &mu,
                           bool covariates, std::string family, bool useRegular) {

        this->covariates = covariates;
        setFamily(family);
        this->regular = useRegular;

        this->y = y;
        y_original = y;

        if(covariates){
            this->z = z;
            z_original = z;
            Z = concatenate(z);
        }

        Y = concatenate(y);

        this->mu = mu;

        for (int i = 0; i < y.size(); i++)
            ycenter.push_back(y[i].array() - mu[i].array());

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

        if(regular)
            return getVarianceRegular();
        if (normal)
            return getVarianceNormal(rvs);
        if (binomial)
            return getVarianceBinomial(rvs);

    }

    inline std::vector<MatrixXd> getX(){

        std::vector<MatrixXd> x;

        for(int i = 0; i < t[0].size(); i++){
            MatrixXd x_i(t[0].xSize(i), size());

            for(int j = 0; j < size(); j++)
                x_i.col(j) = t[j].getX(i);

            x.push_back(x_i);
        }
        return x;
    }

    inline int size() { return t.size(); }
    void bootstrap();
};
