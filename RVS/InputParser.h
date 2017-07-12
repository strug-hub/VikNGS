#pragma once
#include "RVS.h"
#include "MemoryMapped/MemoryMapped.h"
#include <vector>
#include <iostream>  
#include <string>
#include <algorithm>
#include <regex>
#include <map>

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct GenotypeLikelihood {
	bool missing = false;
	double L00;
	double L01;
	double L11;
};

struct VCFLine {
	std::string chr;
	int loc;
	std::string ref;
	std::string alt;
	std::string filter;
	std::vector<int> readDepth;
	std::vector<GenotypeLikelihood> likelihood;
	std::vector<double> P;
	VectorXd expectedGenotype;

	inline bool operator<(VCFLine& line) {
		if (this->chr == line.chr)
			return this->loc < line.loc;
		else
			return this->chr < line.chr;
	}

	inline void print() {
		std::cout << chr;
		std::cout << "\t";
		std::cout << loc;
		std::cout << "\t";
		std::cout << "ref: ";
		std::cout << ref;
		std::cout << "\t";
		std::cout << "alt: ";
		std::cout << alt;
		std::cout << "\n";
	}
};

struct InfoLine {
	std::string id;
	double Y;
	double G;
	std::vector <double> Z;
};

inline bool lineCompare(VCFLine lhs, VCFLine rhs) { return lhs < rhs; }

//InfoParser.cpp
void parseInfo(std::string sampleInfoDir, std::map<std::string, int> &IDmap,
	VectorXd &Y, VectorXd &G, VectorXd &H, MatrixXd &Z);

//VCFParser.cpp
std::map<std::string, int> getSampleIDMap(std::string vcfDir);
std::vector<VCFLine> parseVCFLines(std::string vcfDir);

//VariantFilter.cpp
std::vector<VCFLine> filterVariants(std::vector<VCFLine> variants, VectorXd &G, double missingThreshold);
std::vector<VCFLine> removeDuplicates(std::vector<VCFLine> variants);
std::vector<VCFLine> filterHomozygousVariants(std::vector<VCFLine> &variants);
std::vector<VCFLine> filterMinorAlleleFrequency(std::vector<VCFLine> &variants, double mafCutoff, bool common);

//ParserHelper.cpp
inline std::string extractString(MemoryMapped &charArray, int start, int end);
std::vector<std::string> readLine(MemoryMapped &charArray, int &pos);
inline std::string trim(std::string str);
inline std::vector<std::string> split(std::string s, char sep);
