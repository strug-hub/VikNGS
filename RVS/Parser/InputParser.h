#pragma once
#include "MemoryMapped/MemoryMapped.h"
#include "../Log.h"
#include "../RVS.h"

#include <iostream>  
#include <vector>
#include <string>
#include <algorithm>
#include <map>

#include "../Eigen/Dense"
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
	std::vector<GenotypeLikelihood> likelihood;
	VectorXd P;
	VectorXd expectedGenotype;

	bool valid = true;
	inline bool isValid(){ return valid; }

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

//ParserHelper.cpp
inline std::string extractString(MemoryMapped &charArray, int start, int end);
std::vector<std::string> readLine(MemoryMapped &charArray, int &pos);
inline std::string trim(std::string str);
std::vector<std::string> split(std::string s, char sep);

//InfoParser.cpp
bool parseInfo(std::string sampleInfoDir, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup,
	int highLowCutOff=30);

//VCFParser.cpp
std::map<std::string, int> getSampleIDMap(std::string vcfDir);
std::vector<VCFLine> parseVCFLines(std::string vcfDir);

//VariantFilter.cpp
std::vector<VCFLine> filterVariants(std::vector<VCFLine> variants, VectorXd &G, double missingThreshold);
std::vector<VCFLine> removeDuplicates(std::vector<VCFLine> variants);
std::vector<VCFLine> filterHomozygousVariants(std::vector<VCFLine> &variants);
std::vector<VCFLine> filterMinorAlleleFrequency(std::vector<VCFLine> &variants, double mafCutoff, bool common);

//BEDParser.cpp
std::vector<std::vector<int>> parseIntervals(std::string bedDir, std::vector<VCFLine> variants);

struct File {
	MemoryMapped mmap;
	int pos;
	int maxPos;
	int lineNumber;

	inline void open(std::string directory) {
		mmap.open(directory);
		pos = 0;
		lineNumber = 0;
		maxPos = mmap.size();
	}

	inline void close() {
		mmap.close();
	}

	inline std::string nextLine() {
		int start = pos;

		while (true) {
			pos++;
			if (pos >= maxPos || mmap[pos] == '\n')
				break;
		}
		pos++;
		lineNumber++;

		return extractString(mmap, start, pos - 1);
	}

	inline int getLineNumber() {
		return lineNumber;
	}

	inline bool hasNext() {
		return pos < maxPos;
	}
};
