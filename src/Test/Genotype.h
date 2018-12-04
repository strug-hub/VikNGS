#pragma once
#include "../Math/Math.h"

class Genotype {
    MatrixXd X;
    MatrixXd P;
    GenotypeSource source;

public:

    Genotype(MatrixXd genotypes, MatrixXd alleleFreq, GenotypeSource gtSource) : X(genotypes), P(alleleFreq), source(gtSource) {    }

    //toRemove: 0 = keep, 1 = remove
    inline void filterX(VectorXi& toRemove){
        this->X = extractRows(this->X, toRemove, 0);
    }

    inline void replaceNA(double value){
        this->X = replaceNAN(X, value);
    }

    inline MatrixXd* getX() { return &X; }
    inline int size() { return this->X.cols(); }

    inline VectorXd robustVarVector(){
        VectorXd robustVar(P.rows());
        for (int i = 0; i < P.rows(); i++)
            robustVar[i] = sqrt(calcRobustVar(P(i,1), P(i,2)));

        return robustVar;
    }

    inline VectorXd mafWeightVector(){

        VectorXd maf(P.rows());
        for (int i = 0; i < P.rows(); i++){
            double m = P(i,0) + 0.5*P(i,1);
            if(m >= 1){
                double n = this->X.rows();
                m = (n + 0.5)/(n+1);
            }
            maf[i] = 1/sqrt(m * (1-m));
        }
        return maf;
    }

};
