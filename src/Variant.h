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

struct Variant {
private:
    std::vector<double> pvalues;
    std::vector<std::string> psource;

public:
    std::vector<GenotypeLikelihood> likelihood;
    std::string chr;
    int pos;
    double trueMaf;
    std::string ref;
    std::string alt;
    std::string filter;
    VectorXd P;
    VectorXd expectedGenotype;
    VectorXd trueGenotype;
    VectorXd genotypeCalls;

    Interval interval;

    bool valid = true;
    std::string errorMessage;
    inline bool isValid(){ return valid; }
    inline void setInvalid(std::string message) {
        valid = false;
        errorMessage = message;
        reduceSize();
    }
    std::string getErrorMessage() { return errorMessage; }

    inline std::string toString() {
        std::string t = "\t";
        return chr + t + std::to_string(pos) + t + ref + t + alt;
    }

    inline std::string info() {
        std::string t = "\t";
        return "chr:" + chr + " pos:" + std::to_string(pos) + " ref:" + ref + " alt:" + alt;
    }

    inline void print() {
        std::cout << toString() + "\n";
    }

    inline void reduceSize() {
        likelihood.clear();
        VectorXd empty;
        P = empty;
        expectedGenotype = empty;
    }

    inline void addPval(double pval, std::string psource) {
        this->pvalues.push_back(pval);
        this->psource.push_back(psource);
    }

    inline double getPval(int i) { return pvalues[i]; }
    inline std::string getPvalSource(int i) { return psource[i]; }

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

inline bool variantCompare(Variant lhs, Variant rhs) { return lhs < rhs; }
