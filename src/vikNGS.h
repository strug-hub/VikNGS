#pragma once

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "SampleInfo.h"
#include "Interval.h"
#include "Variant.h"
#include "Request.h"

struct Test;
enum class Genotype;
enum class Statistic;
bool isRare(Statistic s);
enum class Family;
enum class Depth;
enum class Filter;
enum class CollapseType;

//========================================================
// Main object that contains all the information
// about tests and results
//========================================================

struct Data {
    SampleInfo sampleInfo;
    IntervalSet intervals;
    std::vector<Test> tests;
    std::vector<VariantSet> variants;
    double processingTime;
    double evaluationTime;

    inline int size(){ return static_cast<int>(variants.size()); }
};

//========================================================
// Global variable for thread stopping
//========================================================
extern bool STOP_RUNNING_THREAD;

//========================================================
// Main functions that have different implementations
// for command line vs GUI
//========================================================

Data startVikNGS(Request req);
std::vector<VariantSet> processVCF(Request &req, SampleInfo &input);

//========================================================
// Output functions
//========================================================
void outputDebug(std::string line, std::string outputDir);
void outputMatrix(MatrixXd M, std::string filename);
void outputVector(VectorXd V, std::string filename);
void outputFiltered(std::vector<Variant> variants, std::string explain, std::string outputDir);
void outputFiltered(std::vector<std::string> variantInfo, std::vector<int> failCode,
                           std::vector<std::string> codeMap, std::string outputDir);
void outputPvals(std::vector<VariantSet>& variants, std::string outputDir, std::vector<Test>& test);
void initializeOutputFiles (std::string outputDir);
