#pragma once
#include "../Math/Math.h"

class Phenotype {
    VectorXd Y;
    VectorXd Mu;

    Family family;
    MatrixXd Z;

    inline void calculateMu() {
        if(hasCovariates()){
            VectorXd beta = getBeta(Y, Z, family);
            Mu = fitModel(beta, Z, family);
        }
        else
            Mu = VectorXd::Constant(Y.rows(), Y.mean());
    }

    bool isMuCalculated = false;
public:

    Phenotype(MatrixXd phenotype, MatrixXd covariates, Family distribution) : Y(phenotype), Z(covariates), family(distribution) {
        isMuCalculated = false;
        if(Z.rows() < 1 || Z.cols() < 1)
            Z = MatrixXd::Constant(Y.rows(), 1, 1);
    }
    Phenotype(MatrixXd phenotype, Family distribution) : Y(phenotype), family(distribution) {
        isMuCalculated = false;
        Z = MatrixXd::Constant(Y.rows(), 1, 1);
    }

    //toRemove: 0 = keep, 1 = remove
    inline void filterY(VectorXi& toRemove){
        this->Y = extractRows(this->Y, toRemove, 0);
    }
    inline void filterZ(VectorXi& toRemove){
        this->Z = extractRows(this->Z, toRemove, 0);
    }

    inline VectorXd* getY() { return &Y; }
    inline VectorXd* getMu() { if(!isMuCalculated) calculateMu(); return &Mu; }
    inline VectorXd getYCenter() { if(!isMuCalculated) calculateMu(); return (Y-Mu); }
    inline Family getFamily() { return family; }

    inline MatrixXd* getZ() { return &Z; }

    inline int ncov() { return this->Z.cols(); }
    inline bool hasCovariates() {
        int cols = Z.cols();
        return cols > 1;
    }
    inline int size() { return this->Y.rows(); }

};
