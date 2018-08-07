#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "Parser/BED/Interval.h"

using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::MatrixXd;

struct GenotypeLikelihood {
    bool missing = false;
    double L00;
    double L01;
    double L11;

    std::string toString(){
        return "(" + std::to_string(L00) + ", " +
                std::to_string(L01) + ", " +
                std::to_string(L11) + ")";
    }
};

enum class Test { COMMON_LIKELIHOOD_RVS, COMMON_LIKELIHOOD_NORVS,
                  COMMON_REGULAR_TRUE, COMMON_REGULAR_GTCALL,
                  RARE_LIKELIHOOD_RVS, RARE_LIKELIHOOD_NORVS,
                  RARE_REGULAR_TRUE, RARE_REGULAR_GTCALL, NONE };


struct Variant {
private:
    VectorXd expectedGenotype;
    VectorXd trueGenotype;
    VectorXd genotypeCalls;

    Vector3d P;
    Vector3d GL;

    std::string chrom;
    int pos;
    std::string ref;
    std::string alt;

    std::vector<double> pvalues;    
    std::vector<Test> psource;

    bool valid;
    std::string errorMessage;


public:

    Variant(std::string chromosome, std::string position, std::string reference, std::string alternative) :
        chrom(chromosome), pos(position), ref(reference), alt(alternative) {
        valid = true;
    }
    ~Variant() { }

    inline bool isValid(){ return valid; }
    inline void setInvalid(std::string message) {
        valid = false;
        errorMessage = message;
        reduceSize();
    }
    std::string getErrorMessage() { return errorMessage; }


    inline std::string toString() {
        std::string t = "\t";
        return chrom + t + std::to_string(pos) + t + ref + t + alt;
    }
    inline void print() {
        std::cout << toString() + "\n";
    }
    inline std::string info() {
        std::string t = "\t";
        return "chr:" + chr + " pos:" + std::to_string(pos) + " ref:" + ref + " alt:" + alt;
    }

}










    double trueMaf;


    std::vector<int> readDepths;
    std::vector<std::vector<int>> baseCalls;
    std::vector<std::string> vcfCalls;
    std::string format;
    std::vector<std::string> columnUsed;









    inline void reduceSize() {
        likelihood.clear();
        VectorXd empty;
        P = empty;
        expectedGenotype = empty;
    }

    inline void addPval(double pval, Test psource) {
        this->pvalues.push_back(pval);
        this->psource.push_back(psource);
    }

    inline double getPval(int i) { return pvalues[i]; }
    inline std::string getPvalSource(int i) {

        switch(psource[i])
        {
            case Test::COMMON_LIKELIHOOD_RVS :
               return "RVS Genotype Likelihoods - common";
            case Test::COMMON_LIKELIHOOD_NORVS :
               return "Genotype Likelihoods - common";
            case Test::COMMON_REGULAR_TRUE :
               return "True Genotype - common";
            case Test::COMMON_REGULAR_GTCALL :
                return "Regular Genotype Calls - common";
            case Test::RARE_LIKELIHOOD_RVS :
               return "RVS Genotype Likelihoods - rare";
            case Test::RARE_LIKELIHOOD_NORVS :
               return "Genotype Likelihoods - rare";
            case Test::RARE_REGULAR_TRUE :
               return "True Genotype - rare";
            case Test::RARE_REGULAR_GTCALL :
                return "Regular Genotype Calls - rare";
            default:
                  return "???";
         }
    }

    inline std::string getPvalSourceShort(int i) {

        switch(psource[i])
        {
            case Test::COMMON_LIKELIHOOD_RVS :
               return "RVS";
            case Test::COMMON_LIKELIHOOD_NORVS :
               return "no RVS";
            case Test::COMMON_REGULAR_TRUE :
               return "True GT";
            case Test::COMMON_REGULAR_GTCALL :
                return "GT calls";
            case Test::RARE_LIKELIHOOD_RVS :
                return "RVS";
            case Test::RARE_LIKELIHOOD_NORVS :
                return "no RVS";
            case Test::RARE_REGULAR_TRUE :
                return "True GT";
            case Test::RARE_REGULAR_GTCALL :
                return "GT calls";
            default:
                  return "???";
         }
    }

    inline int nPvals() { return this->pvalues.size(); }

    inline bool operator<(Variant& line) {
        if (this->chr == line.chr)
            return this->pos < line.pos;
        else
            return this->chr < line.chr;
    }

    inline MatrixXd getLikelihoodMatrix() {

        MatrixXd M(likelihood.size(), 3);
        for(int i = 0; i < likelihood.size(); i++){
            M(i,0)=likelihood[i].L00;
            M(i,1)=likelihood[i].L01;
            M(i,2)=likelihood[i].L11;
        }

        return M;
    }

    //Uses EM algorithm to estimate the genotype frequencies in the sample.
    inline void calculateGenotypeFrequency() {

        MatrixXd M = getLikelihoodMatrix();

        double p = 0.15;
        double q = 0.15;
        double qn = 1;
        double pn = 0;
        double dn = 0;
        double d = 0;

        int k = 0;

        while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

            d = 1 - p - q;
            Vector3d v = { p, q, d };
            VectorXd pD = M * v;
            VectorXd Ep = M.col(0).array() * (p / pD.array()).array();
            VectorXd Eq = M.col(1).array() * (q / pD.array()).array();
            pn = p;
            qn = q;
            dn = 1 - q - p;
            p = Ep.sum() / Ep.rows();
            q = Eq.sum() / Eq.rows();

            k++;
            if (k == 1000)
                break;
        }

        VectorXd freq(3);
        freq[0] = std::max(0.0, p);
        freq[1] = std::max(0.0, q);
        freq[2] = std::max(0.0, 1 - p - q);

        P=freq;
    }

};

struct CollapsedVariants{

    std::vector<Variant> variants;

}


inline bool variantCompare(Variant lhs, Variant rhs) { return lhs < rhs; }
