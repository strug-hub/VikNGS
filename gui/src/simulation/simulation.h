#pragma once
#include <string>
#include <vector>
#include <map>

#include "../../../src/Math/MathHelper.h"
#include "../../../src/Log.h"
#include "../../../src/RVS.h"
#include "../../../src/Variant.h"
#include "../../../src/Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Vector3d;
using Eigen::DiagonalMatrix;

struct SimulationRequestGroup {
    int index;
    int n; //The number of samples in this group
    int n_increment = 0;
    bool isHrg;

    bool isCaseControl;
    bool isCase;

    double meanDepth;
    double sdDepth;
    double errorRate;

    inline int getSampleSize(int step){ return n + n_increment*step; }

    inline int getStepSize(int step){

        if(step==0)
            return n;
        if(step>0)
            return n_increment;
    }

    inline double generatePhenotype(){

        if(isCaseControl)
             return (isCase ? 1 : 0);

         return 0;
    }

    //for debugging
    void print() {
        std::string s;
        printInfo("n = " + std::to_string(n));
        s = "read depth = ";
        s = s + (isHrg ? "high" : "low");
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
    double mafMin;
    double mafMax;
    int steps;

    std::vector<SimulationRequestGroup> groups;

    std::string test;
    bool useBootstrap;
    int nboot;
    int collapse = 1;
    bool stopEarly;
    bool rvs;
    bool regular;
    bool isCaseControl;


    inline bool underNull(){
        return oddsRatio == 1;
    }

    VectorXd getCaseControlStatus(int step){

        VectorXd y(nsamp(step));
        int index = 0;
        for (int i = 0; i <= step; i++){
            for (int j = 0; j < groups.size(); j++){
                for (int k = 0; k < groups[j].getStepSize(i); k++){
                    y[index] = groups[j].isCase;
                    index++;
                }
            }
        }


        return y;
    }

    VectorXd getGroups(int step){

        VectorXd g(nsamp(step));
        int index = 0;
        for (int i = 0; i <= step; i++){
            for (int j = 0; j < groups.size(); j++){
                for (int k = 0; k < groups[j].getStepSize(i); k++){
                    g[index] = j;
                    index++;
                }
            }
        }
        return g;
    }

    void validate(){

        std::string error_loc = "SIMULATION REQUEST";

        int nsample = 0;
        for (SimulationRequestGroup g : groups)
            nsample += g.n;

        if(isRare())
            if (collapse < 2)
                throwError(error_loc, "Number of variants per collapsed region should be a number greater than 1.", std::to_string(collapse));

        if (nsnp <= 0)
            throwError(error_loc, "Number of variants should be greater than zero.", std::to_string(nsnp));

        if (mafMin > 0.5 || mafMin <= 0)
            throwError(error_loc, "Minimum minor allele frequency should be a value between 0 and 0.5.", std::to_string(mafMin));
        if (mafMax > 0.5 || mafMax <= 0)
            throwError(error_loc, "Maximum minor allele frequency should be a value between 0 and 0.5.", std::to_string(mafMax));
        if (mafMax < mafMin)
            throwError(error_loc, "Maximum minor allele frequency should greater than or equal to minimum minor allele frequency.");

        if (groups.size() < 2)
            throwError(error_loc,"Simulation requires at least two groups.");

        if (oddsRatio <= 0)
            throwError(error_loc, "Odds ratio should be a value greater than zero.", std::to_string(oddsRatio));

        //if (!request.useCommonTest && request.nsnp < 5)
            //throw std::domain_error("Rare test requires number of variants to be at least 5 (value given: " +
                //std::to_string(request.nsnp) + ").");

        if (useBootstrap && nboot < 1)
            throwError(error_loc, "Number of bootstrap iterations should be at least 1.", std::to_string(nboot));

        int ncase = 0;
        int ncontrol = 0;

        for (SimulationRequestGroup g : groups) {

            if (g.isCase)
                ncase++;

            if (!g.isCase)
                ncontrol++;

            if(g.n < 1)
                throwError(error_loc, "Each group should have a sample size of at least 1.");
            if (g.meanDepth <= 0)
                throwError(error_loc, "Mean depth should be a value greater than zero for all groups.");
            if (g.sdDepth <= 0)
                throwError(error_loc, "Read depth standard deviation should be a value greater than zero for all groups.");
            if (g.errorRate < 0)
                throwError(error_loc, "Mean error rate should be a value greater than or equal to zero for all groups.");
        }

        if (ncase < 1 || ncontrol < 1)
            throwError(error_loc, "Case/control simulation requires at least one case and one control group (" +
                std::to_string(ncase) + " case and " +
                std::to_string(ncontrol) + " control groups given).");          

    }

    inline int stepIncrementSize(int step){
        int n = 0;
        for (int i = 0; i < groups.size(); i++)
            n+=groups[i].getStepSize(step);

        return n;
    }

    inline int nsamp(int step){
        int n = 0;
        for (int i = 0; i < groups.size(); i++)
            n += groups[i].getSampleSize(step);
        return n;
    }

    inline int ncase(int step){
        int n = 0;
        for (int i = 0; i < groups.size(); i++)
            if (groups[i].isCase)
                n += groups[i].getSampleSize(step);
        return n;

    }

    inline int ncontrol(int step){
        int n = 0;
        for (int i = 0; i<groups.size(); i++)
            if (!groups[i].isCase)
                n += groups[i].getSampleSize(step);
        return n;
    }

    int isRare(){
        return test == "calpha" || test == "cast";
    }

    //for debugging
    void print() {
        printInfo("nsnp = " + std::to_string(nsnp));
        printInfo("oddsRatio = " + std::to_string(oddsRatio));
        printInfo("MAFmin = " + std::to_string(mafMin));
        printInfo("MAFmax = " + std::to_string(mafMax));

        for (int i = 0; i < groups.size(); i++) {
            printInfo( "group " + std::to_string(i) + ":");
            groups[i].print();
        }
    }
};

std::vector<std::vector<Variant>> startSimulation (SimulationRequest& simReq);
std::vector<TestInput> simulate(SimulationRequest& simReq);

VectorXd simulateMinorAlleleFrequency(int nsnp, double min, double max);
std::vector<VectorXd> simulateG(SimulationRequest& simReq);
std::vector<VectorXd> simulateY(SimulationRequest& simReq);
std::vector<MatrixXd> simulateX(SimulationRequest& simReq, double oddsRatio, VectorXd& maf, int collapse=1);

std::vector<int> generateReadDepths(VectorXd& x, double meanDepth, double depthSd);
std::vector<std::vector<int>> generateBaseCalls(VectorXd& x, double errorRate, std::vector<int>& readDepths);
std::vector<GenotypeLikelihood> generateLikelihoods(VectorXd& x, double errorRate, std::vector<std::vector<int>>& baseCalls);

double inline generateGenotype(double prob_x0, double prob_x1);

VectorXd calculateExpectedGenotypes(std::vector<GenotypeLikelihood>& gl, VectorXd& p);
VectorXd calculateGenotypeCalls(std::vector<GenotypeLikelihood>& gl);
Variant randomVariant();
