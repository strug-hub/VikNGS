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
    bool isHrg;
    bool isCase;
    double meanDepth;
    double sdDepth;
    double errorRate;

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
    int intervalSize;

    std::vector<SimulationRequestGroup> groups;

    std::string test;
    bool useBootstrap;
    int nboot;
    int collapse;
    bool stopEarly;
    bool rvs;
    bool regular;

    int ncase_cache = -1;
    int ncontrol_cache = -1;

    inline bool underNull(){
        return oddsRatio == 1;
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

    inline int ncase(bool useCache=false){
        if(!useCache || this->ncase_cache < 0){
            int n = 0;
            for (int i = 0; i<groups.size(); i++)
                if (groups[i].isCase)
                    n += groups[i].n;
            if(useCache)
                this->ncase_cache = n;
            return n;
        }
        return ncase_cache;
    }

    inline int ncontrol(bool useCache=false){
        if(!useCache || this->ncontrol_cache < 0){
            int n = 0;
            for (int i = 0; i<groups.size(); i++)
                if (!groups[i].isCase)
                    n += groups[i].n;
            if(useCache)
                this->ncontrol_cache = n;
            return n;
        }
        return ncontrol_cache;
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

std::vector<std::vector<Variant>> startSimulation (std::vector<SimulationRequest> simReqs);
std::vector<TestInput> simulate(std::vector<SimulationRequest> simReqs);

VectorXd simulateMinorAlleleFrequency(int nsnp, double min, double max);
VectorXd simulateY(SimulationRequest simReq);
MatrixXd simulateX(SimulationRequest simReq, double oddsRatio, VectorXd maf);
MatrixXd simulateX(SimulationRequest simReq, double oddsRatio, VectorXd maf, int collapse);

Variant generateSeqData(VectorXd x, VectorXd g, std::map<int, SimulationRequestGroup> group, Variant &variant);

std::vector<char> baseCall(std::vector<char> trueGenotype, double error, int readDepth);

double inline generateGenotype(double prob_y, double prob_x0, double prob_x1);

Variant randomVariant();
double pSingle(char base, char true1, char true2, double error);
VectorXd calculateLikelihood(std::vector<char> &bases, double error);
MatrixXd calculateLikelihood2(std::vector<char> &bases, double error);
VectorXd calcEM(MatrixXd M);
VectorXd calculateExpectedGenotypes(std::vector<MatrixXd> M, VectorXd p);

