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

static const int PASS = 0;
static const int SNP_FAIL = 1;
static const int FILTER_FAIL = 2;
static const int MISSING_FAIL = 3;
static const int HOMOZYGOUS_FAIL = 4;
static const int MAF_FAIL = 5;

SampleInfo parseSampleInfo(Request req);

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
std::vector<Variant> parseVCFLines(SampleInfo &input, Request &req);
std::map<std::string, int> getSampleIDMap(std::string vcfDir);
std::vector<std::string> extractHeader(File &vcf);

//SampleParser.cpp
void parseSampleLines(Request req, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup);

//BEDParser.cpp
std::vector<Interval> parseBEDLines(std::string bedDir, bool collapseExon);

//VariantFilter.cpp
int filterVariant(Request &req, Variant &variant, VectorXd &Y, std::string family);
void printFilterResults(Request &req, std::vector<std::string> variantInfo, std::vector<int> failCode, int total);
bool isIn(std::string &vcfLine, int minPos, int maxPos, std::string &chr);
int isIn(std::string &vcfLine, int minPos, int maxPos, std::string &chr, std::vector<Interval> &intervals);

