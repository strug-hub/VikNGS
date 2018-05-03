#pragma once
#include <string>
#include <vector>
#include <map>
#include <QString>

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
    }
};

struct SimulationRequest {
    int npop; //The number of population
    double prevalence; //A decimal between[0, 1], prevalence rate of the disease.

    int nsnp;  //Integer.The number of variants or bases.
    double me; //The mean error rate of sequencing.
    double sde;  //The standard deviation for the error rate.

    double oddsRatio;  //Under H0
    double maf;

    std::vector<SimulationRequestGroup> groups;

    std::string test;
    bool useBootstrap;
    int nboot;
    bool stopEarly;

    void validate(){

        std::string error_loc = "SIMULATION REQUEST";

        int nsample = 0;
        for (SimulationRequestGroup g : groups)
            nsample += g.n;

        if (nsample > npop)
            throwError(error_loc, "Population size should be greater than sample size (summed across groups: " +
                std::to_string(nsample) + ").", std::to_string(npop));

        if (npop <= 0)
            throwError(error_loc, "Population size should be greater than zero.", std::to_string(npop));

        if (nsnp <= 0)
            throwError(error_loc, "Number of variants should be greater than zero.", std::to_string(nsnp));

        if (prevalence > 1 || prevalence <= 0)
            throwError(error_loc, "Prevelance should be a value between 0 and 1.", std::to_string(prevalence));

        if (maf > 0.5 || maf <= 0)
            throwError(error_loc, "Minor allele frequency should be a value between 0 and 0.5.", std::to_string(maf));

        if (groups.size() < 2)
            throwError(error_loc,"Simulation requires at least two groups.");

        if (me <= 0)
            throwError(error_loc, "Mean error rate should be a value greater than zero.", std::to_string(me));

        if (sde <= 0)
            throwError(error_loc, "Error rate standard deviation should be a value greater than zero.", std::to_string(sde));

        if (oddsRatio <= 0)
            throwError(error_loc, "Odds ratio should be a value greater than zero.", std::to_string(me));

        //if (!request.useCommonTest && request.nsnp < 5)
            //throw std::domain_error("Rare test requires number of variants to be at least 5 (value given: " +
                //std::to_string(request.nsnp) + ").");

        if (useBootstrap && nboot < 1)
            throwError(error_loc, "Number of bootstrap iterations should be at least 1.", std::to_string(nboot));

        int ncase = 0;
        int ncontrol = 0;

        int populationCaseIndividuals = floor(npop * prevalence);
        int populationControlIndividuals = npop - populationCaseIndividuals;

        int caseIndividuals = 0;
        int controlIndividuals = 0;

        for (SimulationRequestGroup g : groups) {

            if (g.isCase){
                ncase++;
                caseIndividuals+= g.n;
            }
            if (!g.isCase){
                ncontrol++;
                controlIndividuals += g.n;
            }

            if(g.n < 1)
                throwError(error_loc, "Each group should have a sample size of at least 1.");
            if (g.meanDepth <= 0)
                throwError(error_loc, "Mean depth should be a value greater than zero for all groups.");
            if (g.sdDepth <= 0)
                throwError(error_loc, "Read depth standard deviation should be a value greater than zero for all groups.");
        }


        if(caseIndividuals > populationCaseIndividuals)
            throwError(error_loc, "There are not enough case individuals in the population (Population × Prevalence = " +
                       std::to_string(populationCaseIndividuals) + ") to sample " + std::to_string(caseIndividuals) + " individuals.");
        if(controlIndividuals > populationControlIndividuals)
            throwError(error_loc, "There are not enough control individuals in the population (Population × [1-Prevalence] = " +
                       std::to_string(populationControlIndividuals) + ") to sample " + std::to_string(controlIndividuals) + " individuals.");

        if (ncase < 1 || ncontrol < 1)
            throwError(error_loc, "Case/control simulation requires at least one case and one control group (" +
                std::to_string(ncase) + " case and " +
                std::to_string(ncontrol) + " control groups given).");          

    }

    //for debugging
    void print() {
        printInfo("npop = " + std::to_string(npop));
        printInfo("prevalence = " + std::to_string(prevalence));
        printInfo("nsnp = " + std::to_string(nsnp));
        printInfo("me = " + std::to_string(me));
        printInfo("sde = " + std::to_string(sde));
        printInfo("oddsRatio = " + std::to_string(oddsRatio));
        printInfo("MAF = " + std::to_string(maf));
        printInfo("npop = " + std::to_string(npop));
        printInfo("npop = " + std::to_string(npop));

        for (int i = 0; i < groups.size(); i++) {
            printInfo( "group " + std::to_string(i) + ":");
            groups[i].print();
        }
    }
};
/*
SimulationRequest testSimulationRequest() {
    std::vector<SimulationRequestGroup> groups;
    groups.push_back(newSimulationRequestGroup(0, "300", "control", "low", "10", "5"));
    groups.push_back(newSimulationRequestGroup(0, "200", "case", "high", "40", "7"));
    groups.push_back(newSimulationRequestGroup(0, "100", "control", "high", "50", "6"));

    return newSimulationRequest("3000", "0.1", "100", "0.01", "0.025", "1.4", "0.1", groups, "common", false, 0);
}
*/

std::vector<std::vector<Variant>> startSimulation (std::vector<SimulationRequest> simReqs);
std::vector<TestInput> simulate(std::vector<SimulationRequest> simReqs);

VectorXd simulateMinorAlleleFrequency(int nsnp, double min, double max);
VectorXd simulatePopulationY(int npop, int ncase);
MatrixXd simulatePopulationX(int npop, int ncase, double oddsRatio, VectorXd maf);
VectorXd sampleY(SimulationRequest simReq, VectorXd Y, int nsamp, int ncase, int ncase_pop);
MatrixXd sampleX(SimulationRequest simReq, MatrixXd X, int nsamp, int ncase, int ncase_pop);
Variant generateSeqData(VectorXd x, VectorXd y, VectorXd g, std::map<int, SimulationRequestGroup> group,
                                      double me, double sde, Variant &variant);

std::vector<char> baseCall(std::vector<char> trueGenotype, double error, int readDepth);

double inline generateGenotype(double prob_y, double prob_x0, double prob_x1);

Variant randomVariant();
double pSingle(char base, char true1, char true2, double error);
VectorXd calculateLikelihood(std::vector<char> &bases, double error);
MatrixXd calculateLikelihood2(std::vector<char> &bases, double error);
VectorXd calcEM(MatrixXd M);
VectorXd calculateExpectedGenotypes(std::vector<MatrixXd> M, VectorXd p);

