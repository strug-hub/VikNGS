#pragma once
#include "../../../src/vikNGS.h"
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

    double meanDepth;
    double sdDepth;
    double errorRate;


    inline int getSampleSize(int step){ return n + n_increment*step; }

    inline int getIncreaseSize(int step){
        if(step == 0)
            return n;
        if(step > 0)
            return n_increment;

        return 0;
    }

    inline double generatePhenotype(){
        if(family == Family::BINOMIAL)
             return (isCase ? 1 : 0);

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
    double oddsRatio;
    double r2;

    double mafMin;
    double mafMax;
    int steps;

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
        return (oddsRatio + epsilon > 1 && oddsRatio - epsilon < 1);
    }

    int mxSize;
    bool maxSizeCache = false;
    inline int maxSize(){
        if(!maxSizeCache)
            mxSize = nsamp(steps-1);
        maxSizeCache = true;
        return mxSize;
    }

    VectorXi getGroups(int step){

        VectorXi G(nsamp(step));
        int index = 0;
        for (int i = 0; i <= step; i++){
            for (int j = 0; j < groups.size(); j++){
                for (int k = 0; k < groups[j].getIncreaseSize(i); k++){
                    G[index] = static_cast<int>(j);
                    index++;
                }
            }
        }
        return G;
    }

    void validate(){

        int nsample = 0;
        for (SimulationRequestGroup g : groups)
            nsample += g.n;

        if(isRare(testStatistic))
            if (collapse < 2)
                throwError(ERROR_SOURCE, "Number of variants per collapsed region should be a number greater than 1.", std::to_string(collapse));

        if (nsnp <= 0)
            throwError(ERROR_SOURCE, "Number of variants should be greater than zero.", std::to_string(nsnp));

        if (mafMin > 0.5 || mafMin <= 0)
            throwError(ERROR_SOURCE, "Minimum minor allele frequency should be a value between 0 and 0.5.", std::to_string(mafMin));
        if (mafMax > 0.5 || mafMax <= 0)
            throwError(ERROR_SOURCE, "Maximum minor allele frequency should be a value between 0 and 0.5.", std::to_string(mafMax));
        if (mafMax < mafMin)
            throwError(ERROR_SOURCE, "Maximum minor allele frequency should greater than or equal to minimum minor allele frequency.");

        if (groups.size() < 2)
            throwError(ERROR_SOURCE,"Simulation requires at least two groups.");

        if (nthreads < 1)
            throwError(ERROR_SOURCE,"Number of threads should be at least 1.");

        //if (!request.useCommonTest && request.nsnp < 5)
            //throw std::domain_error("Rare test requires number of variants to be at least 5 (value given: " +
                //std::to_string(request.nsnp) + ").");

        if (useBootstrap && nboot < 1)
            throwError(ERROR_SOURCE, "Number of bootstrap iterations should be at least 1.", std::to_string(nboot));


        if(family == Family::BINOMIAL){
            if (oddsRatio <= 0)
                throwError(ERROR_SOURCE, "Odds ratio should be a value greater than zero.", std::to_string(oddsRatio));

            int ncase = 0;
            int ncontrol = 0;

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
        printInfo("oddsRatio = " + std::to_string(oddsRatio));
        printInfo("MAFmin = " + std::to_string(mafMin));
        printInfo("MAFmax = " + std::to_string(mafMax));

        for (size_t i = 0; i < groups.size(); i++) {
            printInfo( "group " + std::to_string(i) + ":");
            groups[i].print();
        }
    }
};

Data startSimulation(SimulationRequest& simReq);
SampleInfo simulateSampleInfo(SimulationRequest& simReq);
std::vector<VariantSet> simulateVariants(SimulationRequest& simReq);

Variant randomVariant();
VectorXi simulateG(SimulationRequest& simReq);
VectorXd simulateY(SimulationRequest& simReq);
MatrixXd simulateX(SimulationRequest& simReq, double oddsRatio, VectorXd& maf);
double inline generateGenotype(double prob_x0, double prob_x1);

VectorXi generateReadDepths(SimulationRequest& simReq);
std::vector<std::vector<int>> generateBaseCalls(SimulationRequest& simReq, VectorXd& X, VectorXi& readDepths);
std::vector<Vector3d> generateLikelihoods(SimulationRequest& simReq, std::vector<std::vector<int>>& baseCalls);

