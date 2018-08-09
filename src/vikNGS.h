#pragma once

#include "Request.h"
#include "Variant.h"
#include "SampleInfo.h"
#include "Output/OutputHandler.h"

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>  

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Vector3d;
using Eigen::DiagonalMatrix;

//========================================================
// Main object that contains all the information
// about tests and results
//========================================================

struct Data {
    std::vector<Variant> variants;
    SampleInfo sampleInfo;

    inline int size(){ return variants.size(); }
};

enum class Test { COMMON_LIKELIHOOD_RVS, COMMON_LIKELIHOOD_NORVS,
                  COMMON_REGULAR_TRUE, COMMON_REGULAR_GTCALL,
                  RARE_LIKELIHOOD_RVS_SKAT, RARE_LIKELIHOOD_NORVS_SKAT,
                  RARE_REGULAR_TRUE_SKAT, RARE_REGULAR_GTCALL_SKAT,
                  RARE_LIKELIHOOD_RVS_CALPHA, RARE_LIKELIHOOD_NORVS_CALPHA,
                  RARE_REGULAR_TRUE_CALPHA, RARE_REGULAR_GTCALL_CALPHA, NONE };

enum class CollapseType { COLLAPSE_K, COLLAPSE_GENE, COLLAPSE_EXON, NONE };
enum class Family { NORMAL, BINOMIAL, NONE };
enum class Depth { HIGH, LOW };



//========================================================
// Main functions that have different implementations
// for command line vs GUI
//========================================================

Data startVikNGS(Request req);
std::vector<Variant> runTest(SampleInfo &input, Request &req);

//========================================================
// Global variable for thread stopping
//========================================================

extern bool STOP_RUNNING_THREAD;

//========================================================
// Logging functions
//========================================================

void printInfo(std::string message);
void printWarning(std::string message);
void printError(std::string message);
void throwError(std::string source, std::string message);
void throwError(std::string source, std::string message, std::string valueGiven);
void printWarning(std::string source, std::string message, std::string valueGiven);
void printWarning(std::string source, std::string message);


class variant_exception: public std::exception{
    std::string message;

public:
    variant_exception(std::string msg){
        this->message= msg;
    }

  virtual const char* what() const throw()
  {
        std::string what= "Error in single variant p-value computation: " + this->message;
        return what.c_str();
  }
};


std::vector<Variant> processVCF(SampleInfo input, Request req);

//CommonTest.cpp
std::vector<Variant> runCommonTest(Request &req, SampleInfo &input);
//RareTest.cpp
std::vector<Variant> runRareTest(Request req, SampleInfo input);
