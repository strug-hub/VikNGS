#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "Parser/BED/Interval.h"

using Eigen::VectorXd;

struct GenotypeLikelihood {
    bool missing = false;
    double L00;
    double L01;
    double L11;
};

struct Variant {
    double pvalue;
    std::string chr;
    int pos;
    std::string ref;
    std::string alt;
    std::string filter;
    std::vector<GenotypeLikelihood> likelihood;
    VectorXd P;
    VectorXd expectedGenotype;

    Interval interval;

    bool valid = true;
    std::string errorMessage;
    inline bool isValid(){ return valid; }
    inline void setInvalid(std::string message) {
        valid = false;
        errorMessage = message;
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

    inline bool operator<(Variant& line) {
        if (this->chr == line.chr)
            return this->pos < line.pos;
        else
            return this->chr < line.chr;
    }
};

inline bool variantCompare(Variant lhs, Variant rhs) { return lhs < rhs; }
