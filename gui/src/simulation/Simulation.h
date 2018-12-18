#pragma once
#include "../src/vikNGS.h"
#include "../src/Log.h"
#include <string>
#include <vector>
#include <map>

static const std::string ERROR_SOURCE = "SIMULATION_REQUEST";

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Vector3d;
using Eigen::DiagonalMatrix;

struct SimulationRequestGroup {

    Depth readDepth;

    int index;
    int n; //The number of samples in this group
    int n_increment = 0;

    Family family;
    bool isCase;
    double normalMean;
    double normalSd;

    double meanDepth;
    double sdDepth;
    double errorRate;

    inline std::string getCohort(){
        if(family == Family::BINOMIAL)
            return isCase ? "case" : "control";
        else if (family == Family::NORMAL)
            return "normal";

        return "none";
    }

    inline int getSampleSize(int step){ return n + n_increment*step; }

    inline int getIncreaseSize(int step){
        if(step < 1e-12)
            return n;
        if(step > 0)
            return n_increment;


        return 0;
    }

    inline double generatePhenotype(){
        if(family == Family::BINOMIAL)
             return (isCase ? 1 : 0);
        if(family == Family::NORMAL)
            return normalMean;
         return 0;
    }

    //for debugging
    void print() {
        std::string s;
        printInfo("n = " + std::to_string(n));
        s = "read depth = ";
        s = s + (readDepth == Depth::HIGH ? "high" : "low");
        printInfo(s);
        s = "group = ";
        s = s + (isCase ? "case" : "control");
        printInfo("mean depth = " + std::to_string(meanDepth));
        printInfo("depth sd = " + std::to_string(sdDepth));
        printInfo("error rate = " + std::to_string(errorRate));
    }
};

struct SimulationRequest {

    int nsnp;
    double effectSize;

    double mafMin;
    double mafMax;

    int steps;

    //for covariates...
    double covariate = -1;
    bool corX = true; //otherwise corY = true
    //...

    std::vector<SimulationRequestGroup> groups;

    Family family;
    Statistic testStatistic;
    bool useBootstrap;
    int nboot;
    int collapse = 1;
    bool stopEarly;
    int nthreads = 1;

    inline bool underNull(){
        double epsilon = 1e-8;

        if(family == Family::BINOMIAL)
            return (effectSize + epsilon > 1 && effectSize - epsilon < 1);
        else if (family == Family::NORMAL)
            return (effectSize + epsilon > 0 && effectSize - epsilon < 0);

        else return false;
    }

    int mxSize;
    bool maxSizeCache = false;
    inline int maxSize(){
        if(!maxSizeCache)
            mxSize = nsamp(steps-1);
        maxSizeCache = true;
        return mxSize;
    }

    VectorXi getGroupVector(){
        VectorXi G(maxSize());
        int index = 0;

        for(int i = 0; i < steps; i++)
            for(size_t j = 0; j < groups.size(); j++){
                int n = groups[j].getIncreaseSize(i);
                for (int k = 0; k < n; k++){
                    G[index] = groups[j].index;
                    index++;
                }
            }

        return G;
    }

    std::map<int, Depth> getReadDepthMap(){
        std::map<int, Depth> groupDepth;
        for(size_t j = 0; j < groups.size(); j++)
            groupDepth[groups[j].index] = groups[j].readDepth;
        return groupDepth;
    }

    void validate(){

        int nsample = 0;
        for (SimulationRequestGroup g : groups)
            nsample += g.n;

        if(isRare(testStatistic)){
            if (collapse < 2)
                throwError(ERROR_SOURCE, "Number of variants per collapsed region should be a number greater than 1.", std::to_string(collapse));
            if (nsnp % collapse != 0)
                throwError(ERROR_SOURCE, "Number of variants should be divisible by the collapse size (" + std::to_string(collapse) + ").", std::to_string(nsnp));
        }

        if (nsnp <= 0)
            throwError(ERROR_SOURCE, "Number of variants should be greater than zero.", std::to_string(nsnp));

        if (mafMin > 0.5 || mafMin <= 0)
            throwError(ERROR_SOURCE, "Minimum minor allele frequency should be a value between 0 and 0.5.", std::to_string(mafMin));
        if (mafMax > 0.5 || mafMax <= 0)
            throwError(ERROR_SOURCE, "Maximum minor allele frequency should be a value between 0 and 0.5.", std::to_string(mafMax));
        if (mafMax < mafMin)
            throwError(ERROR_SOURCE, "Maximum minor allele frequency should greater than or equal to minimum minor allele frequency.");
        if (nthreads < 1)
            throwError(ERROR_SOURCE,"Number of threads should be at least 1.");

        //if (!request.useCommonTest && request.nsnp < 5)
            //throw std::domain_error("Rare test requires number of variants to be at least 5 (value given: " +
                //std::to_string(request.nsnp) + ").");

        if (useBootstrap && nboot < 1)
            throwError(ERROR_SOURCE, "Number of bootstrap iterations should be at least 1.", std::to_string(nboot));


        if(family == Family::BINOMIAL){
            if (effectSize <= 0)
                throwError(ERROR_SOURCE, "Odds ratio should be a value greater than zero.", std::to_string(effectSize));

            int ncase = 0;
            int ncontrol = 0;
            if(groups.size() < 2)
                throwError(ERROR_SOURCE, "At least one group should be specified in the group table in order to simulate case-control data.");

            for (SimulationRequestGroup g : groups) {

                if (g.isCase)
                    ncase++;

                if (!g.isCase)
                    ncontrol++;

                if(g.n < 1)
                    throwError(ERROR_SOURCE, "Each group should have a sample size of at least 1.");
                if (g.meanDepth <= 0)
                    throwError(ERROR_SOURCE, "Mean depth should be a value greater than zero for all groups.");
                if (g.sdDepth <= 0)
                    throwError(ERROR_SOURCE, "Read depth standard deviation should be a value greater than zero for all groups.");
                if (g.errorRate < 0)
                    throwError(ERROR_SOURCE, "Mean error rate should be a value greater than or equal to zero for all groups.");
            }

            if (ncase < 1 || ncontrol < 1)
                throwError(ERROR_SOURCE, "Case/control simulation requires at least one case and one control group (" +
                    std::to_string(ncase) + " case and " +
                    std::to_string(ncontrol) + " control groups given).");
        }
        if(family == Family::NORMAL){
            if (effectSize < 0 || effectSize > 1)
                throwError(ERROR_SOURCE, "R² (variance explained) value should be between zero and one.", std::to_string(effectSize));

            if(groups.size() < 1)
                throwError(ERROR_SOURCE, "At least one group should be specified in the group table in order to simulate data.");

            for (SimulationRequestGroup g : groups) {

                if(g.n < 1)
                    throwError(ERROR_SOURCE, "Each group should have a sample size of at least 1.");
                if (g.meanDepth <= 0)
                    throwError(ERROR_SOURCE, "Mean depth should be a value greater than zero for all groups.");
                if (g.sdDepth <= 0)
                    throwError(ERROR_SOURCE, "Read depth standard deviation should be a value greater than zero for all groups.");
                if (g.errorRate < 0)
                    throwError(ERROR_SOURCE, "Mean error rate should be a value greater than or equal to zero for all groups.");
            }
        }
    }

    inline int stepIncrementSize(int step){
        int n = 0;
        for (size_t i = 0; i < groups.size(); i++)
            n+=groups[i].getIncreaseSize(step);

        return n;
    }

    inline int nsamp(int step){
        int n = 0;
        for (size_t i = 0; i < groups.size(); i++)
            n += groups[i].getSampleSize(step);
        return n;
    }

    inline int ncase(int step){
        int n = 0;
        for (size_t i = 0; i < groups.size(); i++)
            if (groups[i].isCase)
                n += groups[i].getSampleSize(step);
        return n;

    }

    inline int ncontrol(int step){
        int n = 0;
        for (size_t i = 0; i<groups.size(); i++)
            if (!groups[i].isCase)
                n += groups[i].getSampleSize(step);
        return n;
    }


    //for debugging
    void print() {
        printInfo("nsnp = " + std::to_string(nsnp));
        printInfo("effectSize = " + std::to_string(effectSize));
        printInfo("MAFmin = " + std::to_string(mafMin));
        printInfo("MAFmax = " + std::to_string(mafMax));

        for (size_t i = 0; i < groups.size(); i++) {
            printInfo( "group " + std::to_string(i) + ":");
            groups[i].print();
        }
    }
};

Data startSimulation(SimulationRequest& simReq);
std::vector<VariantSet> simulateVariants(SimulationRequest& simReq);
std::vector<VariantSet> simulateVariant(SimulationRequest& simReq);

Variant randomVariant();
MatrixXd simulateYCaseControl(SimulationRequest& simReq);
MatrixXd simulateYNormal(SimulationRequest& simReq, MatrixXd& X, VectorXd& mafs);
MatrixXd simulateZ(MatrixXd M, SimulationRequest& simReq);
VectorXd generateMafs(SimulationRequest& simReq);
MatrixXd simulateXCaseControl(SimulationRequest& simReq, VectorXd& mafs);
MatrixXd simulateXNormal(SimulationRequest& simReq, VectorXd& mafs);
std::vector<std::vector<Vector3d>> simulateSequencing(SimulationRequest& simReq, MatrixXd& Xtrue);

double generateGenotype(double prob_x0, double prob_x1);

VectorXi generateReadDepths(SimulationRequest& simReq);
std::vector<std::vector<int>> generateBaseCalls(SimulationRequest& simReq, VectorXd& X, VectorXi& readDepths);
std::vector<Vector3d> generateLikelihoods(SimulationRequest& simReq, std::vector<std::vector<int>>& baseCalls);

