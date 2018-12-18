#pragma once
#include "../Math/EigenStructures.h"

#include <map>
#include <vector>

enum class CollapseType;
enum class Depth;
struct File;
struct Variant;
struct Interval;
struct IntervalSet;
struct SampleInfo;
struct Request;

SampleInfo parseSampleInfo(Request & req);

//SampleParser.cpp
bool validateSampleIDs(std::string sampleDir, std::map<std::string, int> &IDmap);
VectorXd parseSamplePhenotype(std::string sampleDir, std::map<std::string, int> &IDmap);
VectorXi parseSampleGroupID(std::string sampleDir, std::map<std::string, int> &IDmap);
std::map<int, Depth> parseSampleReadDepth(std::string sampleDir, std::map<std::string, int> &IDmap, VectorXi& G, int highLowCutOff);
MatrixXd parseSampleCovariates(std::string sampleDir, std::map<std::string, int> &IDmap);

static const int CHROM = 0;
static const int POS = 1;
static const int ID = 2;
static const int REF = 3;
static const int ALT = 4;
static const int FILTER = 6;
static const int FORMAT = 8;
static const char VCF_SEP = '\t';

//VCFParser.cpp
std::map<std::string, int> getSampleIDMap(std::string vcfDir);
std::vector<std::string> extractHeader(File &vcf);
std::string extractHeaderLine(File &vcf);
Variant constructVariant(std::vector<std::string> &columns, bool calculateExpected, bool calculateCalls, bool getVCFCalls);
Vector3d getGenotypeLikelihood(std::string &column, int indexPL, int indexGL, int indexGT);
double getVCFGenotypeCall(std::string &column, int indexGT);

static const char BED_SEP = '\t';

//BEDParser.cpp
IntervalSet parseBEDLines(std::string bedDir, CollapseType collapse);
Interval lineToGene(std::vector<std::string> lineSplit, int lineNumber);
std::vector<Interval> lineToExons(std::vector<std::string> lineSplit, int lineNumber);

//StringTools.cpp
std::vector<std::vector<int>> collapseEveryK(int k, int n);
std::string trim(std::string str);
std::vector<std::string> splitString(std::string &s, char sep);
std::vector<std::string> splitString(std::string &s, char sep, int stop);

//VariantFilter.cpp
void printFilterResults(Request &req, std::vector<std::string> variantInfo, std::vector<int> failCode, int total);
