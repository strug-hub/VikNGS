#pragma once
#include "../vikNGS.h"
#include "../Math/Math.h"
#include "Group.h"

class TestObject {

    MatrixXd X;
    VectorXd Y;
    MatrixXd Z;
    MatrixXd P;
    Family family;
    Group& group;

    VectorXd Ycenter;
    VectorXd Ycenter_original;

    VectorXd MU;

    MatrixXd Xboot;
    VectorXd Yboot;
    MatrixXd Zboot;

    inline void calculateYCenter() {

        if(hasCovariates()){
            VectorXd beta = getBeta(Y, Z, family);
            this->MU = fitModel(beta, Z, family);
        }
        else{
            double ybar = Y.mean();
            this->MU = VectorXd::Constant(this->Y.rows(), ybar);
        }

        Ycenter = Y - MU;
    }

public:

    TestObject(MatrixXd& genotypes, VectorXd& phenotypes, MatrixXd& covariates,
               MatrixXd& frequency, Family distribution, Group& groups,  bool rareVariant) :
        X(genotypes), Y(phenotypes), Z(covariates), P(frequency), family(distribution), group(groups) {

        //Filter NAN
        VectorXi toRemove = whereNAN(Y);
        if(!rareVariant)
            toRemove = toRemove + whereNAN(X);
        if(hasCovariates())
            toRemove = toRemove + whereNAN(Z);

        X = extractRows(X, toRemove, 0);
        Y = extractRows(Y, toRemove, 0);
        if(hasCovariates())
            Z = extractRows(Z, toRemove, 0);
        else
            Z = MatrixXd::Constant(X.rows(),1,1);

        groups.filterG(toRemove);

        //Replace NAN with 0 if rare
        if(rareVariant)
            X = replaceNAN(X, 0);

        calculateYCenter();

        bootstrapped = false;
        groupVectorCache = false;
        XcenterCache = false;
    }

    inline bool hasCovariates() { return this->Z.cols() > 1; }
    inline double getRobustVar(int i=0){
        return calcRobustVar(P(i,1), P(i,2));
    }
    inline VectorXd robustVarVector(){
        VectorXd robustVar(P.rows());
        for (int i = 0; i < P.rows(); i++)
            robustVar[i] = sqrt(getRobustVar(i));

        return robustVar;
    }
    inline VectorXd mafWeightVector(){

        VectorXd maf(P.rows());
        for (int i = 0; i < P.rows(); i++){
            double m = P(i,0) + 0.5*P(i,1);
            if(m >= 1){
                double n = this->Y.rows();
                m = (n + 0.5)/(n+1);
            }
            maf[i] = 1/sqrt(m * (1-m));
        }


        return maf;
    }

    inline MatrixXd* getX(){ return (bootstrapped) ? &Xboot : &X; }
    inline VectorXd getX(int i){ return (bootstrapped) ? Xboot.col(i) : X.col(i); }
    inline VectorXd* getY(){ return (bootstrapped) ? &Yboot : &Y; }
    inline MatrixXd* getZ(){ return (bootstrapped) ? &Zboot : &Z; }
    inline Group* getGroup(){ return &group; }
    inline MatrixXd* getP(){ return &P; }
    inline VectorXd* getMU(){ return &MU; }
    inline VectorXd* getYcenter(){ return &Ycenter; }

    inline void bootstrap(Test& test, Family family) {

       if(test.isExpectedGenotypes()){
           calculateXcenter();
           calculateGroupVector();

           if(family == Family::NORMAL)
               normalBootstrap();
           if(family == Family::BINOMIAL)
               binomialBootstrap();
       }
       else
           permute();

       bootstrapped = true;
       calculateYCenterBoot();
   }


//bootstrap functions
private:

    bool bootstrapped;

    inline void permute() {

        if(!bootstrapped)
            Xboot = X;

        Yboot = shuffleWithoutReplacement(Y);
        //Xboot = shuffleColumnwiseWithoutReplacement(Xcenter);

        if(hasCovariates())
           Zboot = shuffleColumnwiseWithoutReplacement(Z);
    }

    void binomialBootstrap() {
        if(!bootstrapped)
            Yboot = Y;
        Xboot = groupwiseShuffleWithReplacement(Xcenter, *group.getG(), groupVector);

        //todo?
        if(hasCovariates())
            Zboot = groupwiseShuffleWithReplacement(Z, *group.getG(), groupVector);
    }

    void normalBootstrap() {

       if(hasCovariates()){

           if(!bootstrapped){
               Xboot = X;
               Zboot = Z;
               Ycenter_original = Ycenter;
           }

           VectorXd residuals = groupwiseShuffleWithoutReplacement(Ycenter_original, *group.getG(), groupVector);
           Yboot = Y - Ycenter_original + residuals;

       }
       else{
           Yboot = groupwiseShuffleWithoutReplacement(Y, *group.getG(), groupVector);
           if(!bootstrapped)
               Xboot = X;
       }
    }

    bool groupVectorCache;
    std::map<int, std::vector<int>> groupVector;

    inline void calculateGroupVector() {

        if(groupVectorCache)
            return;

        for(int i = 0; i < group.size(); i++){

            if(groupVector.count(group[i]) < 1){
                std::vector<int> v;
                v.push_back(i);
                groupVector[group[i]] = v;
            }
            else
                groupVector[group[i]].push_back(i);
        }

        groupVectorCache = true;
    }

    bool XcenterCache;
    MatrixXd Xcenter;
    inline void calculateXcenter() {
        if(XcenterCache)
            return;

        Xcenter = subtractGroupMean(X, *group.getG());
        XcenterCache = true;
    }

    inline void calculateYCenterBoot() {

        if(hasCovariates()){
            VectorXd beta = getBeta(Yboot, Zboot, family);
            this->MU = fitModel(beta, Zboot, family);
        }
        else{
            double ybar = Yboot.mean();
            this->MU = VectorXd::Constant(Yboot.rows(), ybar);
        }

        Ycenter = Yboot - MU;
    }
};
