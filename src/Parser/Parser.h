#pragma once
#include "../vikNGS.h"
#include "../Request.h"
#include "File.h"
#include "../SampleInfo.h"
#include "../Variant.h"
#include "../Math/MathHelper.h"
#include "BED/Interval.h"
#include "VCF/VCFParserUtils.h"

#include <iostream>  
#include <vector>
#include <string>
#include <algorithm>
#include <map>


enum class Filter { VALID, NOT_SNP, NO_PASS, MISSING_DATA, NO_VARIATION, MAF };


SampleInfo parseSampleInfo(Request req);


static const int ID_COL = 0;
static const int PHENOTYPE_COL = 1;
static const int GROUP_COL = 2;
static const int DEPTH_COL = 3;
static const int COV_COL = 4;
static const char SAMPLE_SEP = '\t';

//SampleInfoParser.cpp
bool validateSampleIDs(std::string sampleDir, std::map<std::string, int> &IDmap);
VectorXd parseSamplePhenotype(std::string sampleDir, std::map<std::string, int> &IDmap);
VectorXd parseSampleGroupID(std::string sampleDir, std::map<std::string, int> &IDmap);
std::map<int, Depth> parseSampleReadDepth(std::string sampleDir, std::map<std::string, int> &IDmap, VectorXd &G, int highLowCutOff);
MatrixXd parseSampleCovariates(std::string sampleDir, std::map<std::string, int> &IDmap);


//VCFParser.cpp


//ParserTools.cpp
std::vector<std::vector<int>> collapseEveryK(int k, int n);
std::string extractString(MemoryMapped &charArray, int start, int end);
std::string trim(std::string str);
std::vector<std::string> split(std::string &s, char sep);
void calculateExpectedGenotypes(Variant &variants);
VectorXd calcEG(std::vector<GenotypeLikelihood> &likelihood, VectorXd &p);
VectorXd calcEM(std::vector<GenotypeLikelihood> &likelihood);
std::string determineFamily(VectorXd Y);

//VCFParser.cpp
std::map<std::string, int> getSampleIDMap(std::string vcfDir);
std::vector<Variant> parseVCFLines(SampleInfo &input, Request &req);
std::map<std::string, int> getSampleIDMap(std::string vcfDir);
std::vector<std::string> extractHeader(File &vcf);


//BEDParser.cpp
std::vector<Interval> parseBEDLines(std::string bedDir, bool collapseExon);

//VariantFilter.cpp
int filterVariant(Request &req, Variant &variant, VectorXd &Y, std::string family);
void printFilterResults(Request &req, std::vector<std::string> variantInfo, std::vector<int> failCode, int total);
bool isIn(std::string &vcfLine, int minPos, int maxPos, std::string &chr);
int isIn(std::string &vcfLine, int minPos, int maxPos, std::string &chr, std::vector<Interval> &intervals);

