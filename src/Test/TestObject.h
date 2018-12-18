#pragma once
#include "../vikNGS.h"
#include "../Math/Math.h"
#include "Group.h"
#include "Genotype.h"
#include "Phenotype.h"

class TestObject {

    Genotype& geno;
    Phenotype& pheno;
    Group& group;

    VectorXd Ycenter;
    VectorXd Ycenter_original;

    VectorXd MU;

    MatrixXd Xboot;
    VectorXd Yboot;
    MatrixXd Zboot;



public:

    TestObject(Genotype& genotype, Phenotype& phenotype, Group& groups, bool rareVariant) :
        geno(genotype), pheno(phenotype), group(groups) {

        //Filter NAN
        VectorXi toRemove = whereNAN(*pheno.getY());
        if(!rareVariant)
            toRemove = toRemove + whereNAN(*geno.getX());
        if(pheno.hasCovariates())
            toRemove = toRemove + whereNAN(*pheno.getZ());

        geno.filterX(toRemove);
        group.filterG(toRemove);
        pheno.filterY(toRemove);
        if(pheno.hasCovariates())
            pheno.filterZ(toRemove);

        Zboot = *pheno.getZ();
        Ycenter = pheno.getYCenter();

        //Replace NAN with 0 if rare
        if(rareVariant)
            geno.replaceNA(0);

        bootstrapped = false;
        groupVectorCache = false;
        XcenterCache = false;
    }

    inline VectorXd robustVarVector(){
        return geno.robustVarVector();
    }
    inline VectorXd mafWeightVector(){

       return geno.mafWeightVector();
    }

    inline MatrixXd* getX(){ return (bootstrapped) ? &Xboot : geno.getX(); }
    inline VectorXd* getY(){ return (bootstrapped) ? &Yboot : pheno.getY(); }
    inline MatrixXd* getZ(){ return (bootstrapped) ? &Zboot : pheno.getZ(); }
    inline Group* getGroup(){ return &group; }
    inline VectorXd* getMU(){ return pheno.getMu(); }
    inline VectorXd* getYcenter(){ return &Ycenter; }

    inline void bootstrap(TestSettings& test, Family family) {

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
            Xboot = *geno.getX();

        Yboot = shuffleWithoutReplacement(*pheno.getY());
        //Xboot = shuffleColumnwiseWithoutReplacement(*geno.getX());

        if(pheno.hasCovariates())
           Zboot = shuffleColumnwiseWithoutReplacement(*pheno.getZ());
    }

    void binomialBootstrap() {
        if(!bootstrapped)
            Yboot = *pheno.getY();
        Xboot = groupwiseShuffleWithReplacement(Xcenter, *group.getG(), groupVector);


        //todo?
        if(pheno.hasCovariates())
            Zboot = groupwiseShuffleWithReplacement(*pheno.getZ(), *group.getG(), groupVector);
    }

    void normalBootstrap() {


       if(pheno.hasCovariates()){

           if(!bootstrapped){
               Xboot = *geno.getX();
               Zboot = *pheno.getZ();
               Ycenter_original = Ycenter;
           }

           VectorXd residuals = groupwiseShuffleWithoutReplacement(Ycenter_original, *group.getG(), groupVector);
           Yboot = *pheno.getY() - Ycenter_original + residuals;

       }
       else{
           Yboot = groupwiseShuffleWithoutReplacement(*pheno.getY(), *group.getG(), groupVector);
           if(!bootstrapped)
               Xboot = *geno.getX();
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

        Xcenter = subtractGroupMean(*geno.getX(), *group.getG());
        XcenterCache = true;
    }

    inline void calculateYCenterBoot() {

        if(pheno.hasCovariates()){
            VectorXd beta = getBeta(Yboot, Zboot, pheno.getFamily());
            this->MU = fitModel(beta, Zboot, pheno.getFamily());
        }
        else{
            double ybar = Yboot.mean();
            this->MU = VectorXd::Constant(Yboot.rows(), ybar);
        }

        Ycenter = Yboot - MU;
    }
};
