#pragma once
#include "../Math/Math.h"

class Genotype {
    MatrixXd X;

    MatrixXd Xboot;
    bool bootstrapped;
    MatrixXd P;
    VectorXd robustVar;
    GenotypeSource source;

    void robustVarVector(){
        VectorXd var(P.rows());
        for (int i = 0; i < P.rows(); i++)
            calcRobustVar(P(i,1), P(i,2));

        this->robustVar = var;
    }

public:

    Genotype(MatrixXd genotypes, MatrixXd alleleFreq, GenotypeSource gtSource) : X(genotypes), P(alleleFreq), source(gtSource) {
        calculateRobustVar();
        bootstrapped = false;
    }

    //toRemove: 0 = keep, 1 = remove
    inline void filterX(VectorXi& toRemove){
        this->X = extractRows(this->X, toRemove, 0);
    }

    inline MatrixXd* getX() { return (bootstrapped) ? &Xboot : &X; }
    inline int size() { return this->X.cols(); }

};
