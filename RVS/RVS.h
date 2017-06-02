#pragma once
#include "stdafx.h"
#include "MemoryMapped/MemoryMapped.h"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>  
#include <regex>

//========================================================
// Group struct
//========================================================

struct Group {
	std::string ID = "";
	std::vector<size_t> index;
	//double n = 0;
	//double meanY = 0;
	//double meanX = 0;
	//double varX = 0;
	bool hrg = false;
};

//========================================================
// Sample struct
//========================================================

struct Sample {
	std::string ID = "";
	std::string groupID = "";
	double y = 0;
	bool hrg = false;
	std::vector<double> numeric_cov;
	std::vector<std::string> factor_cov;
};


//========================================================
// SNP struct
//========================================================

struct genotypeLikelihood {
	bool valid = false;
	double L00;
	double L01;
	double L11;
};

struct SNP {
	std::string chr;
	int loc;

	std::vector<genotypeLikelihood> gl;
	std::vector<int> rd;

	std::vector<double> p;
	std::vector<double> EG;

	double maf = NAN;

	inline bool operator<(SNP& snp2) { return this->loc < snp2.loc; }
	inline double L00(int index) { return this->gl[index].L00; }
	inline double L01(int index) { return this->gl[index].L01; }
	inline double L11(int index) { return this->gl[index].L11; }

	SNP clone() {
		SNP newSNP;
		newSNP.chr = this->chr;
		newSNP.loc = this->loc;
		newSNP.gl = this->gl;
		newSNP.p = this->p;
		newSNP.EG = this->EG;
		newSNP.maf = this->maf;
		return newSNP;
	}

};

//========================================================
// Interval struct
//========================================================

struct Interval {
	int start = 0;
	int end = 0;
	std::vector<int> indexes;
	std::string chr;

	std::vector<SNP> getInInterval(std::vector<SNP> &snps) {
		std::vector<SNP> in;
		for (size_t i = 0; i < indexes.size(); i++)
			in.push_back(snps[indexes[i]]);

		return in;
	}

	bool addIfIn(SNP snp, int index) {

		if (snp.loc >= start && snp.loc <= end && chr.compare(snp.chr) == 0) {
			indexes.push_back(index);
			return true;
		}
		
		return false;
	}

};


inline bool locCompare(SNP lhs, SNP rhs) { return lhs < rhs; }

//========================================================
// functions
//========================================================

//HelperVCF.cpp
std::vector<Interval> getIntervals(std::string);
std::vector<std::string> parseHeader(MemoryMapped &, int &);
std::vector<Sample> getSampleInfo(std::string, std::string, int);
SNP initSNP(MemoryMapped &, std::vector<int>, int);
std::vector<SNP> parseAndFilter(std::string, int, double, std::vector<Sample> &);
std::vector<double> calcEM(SNP &);
std::vector<double> calcEG(SNP &);
void collapseVariants(std::vector<SNP> &, std::vector<Interval> &);

//Statistics.cpp
double meanX(SNP &, Group &);
double meanY(std::vector<Sample> &, SNP &);
double varX(SNP &, Group &);
double variance(std::vector<double> &);
double chiSquareOneDOF(double);
std::vector<double> randomSample(std::vector<double> &, int);
std::vector<double> CovariateRegression(SNP &, std::vector<Sample> &, std::string = "norm");

//CommonTest.cpp
std::vector<double> RVSasy(std::vector<SNP> &, std::vector<Sample> &, std::vector<Group> &, bool = true);
std::vector<double> RVSbtrap(std::vector<SNP> &, std::vector<Sample> &, std::vector<Group> &, int, bool, bool = true);

//RareTest.cpp
std::vector<double> RVSrare(std::vector<SNP> &, std::vector<Sample> &, std::vector<Group> &, int, bool = true, int = 5, int = 1, int = 1);

//CompQuadForm.cpp
double qfc(std::vector<double>, double, int);

//========================================================
// inline functions
//========================================================

/**
Removes whitespace from a string

@param str String to remove whitespace from.
@return str without whitespace.
*/
inline std::string trim(std::string str)
{
	str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
	return str;
}

/**
Checks if a string is a number

@param str string to check
@return true if str is numeric
*/
inline bool isNumeric(const std::string& str) {
	return (std::regex_match(str, std::regex("-?[1234567890]+(\.[1234567890]+)?")));
}
//-?\d+(\.\d+)?

/**
Extracts string from a MemoryMapped class

@param charArray MemoryMapped to extract string from.
@param start Index to start extracting string from.
@param end Index to stop extracting string from.
@return Extracted string.
*/
inline std::string getString(MemoryMapped &charArray, int start, int end) {
	std::string ret;
	for (; start < end; start++) { ret += charArray[start]; }
	return ret;
}

/**
Separates a string into a vector, splitting at every postion with
the sep character

@param s String to split.
@param sep Character to split the string at.
@return Split string.
*/
inline std::vector<std::string> split(std::string s, char sep) {
	std::vector<std::string> split;
	int start = 0;
	for (int i = 0; i <= s.length(); i++) {
		if (s[i] == sep || i == s.length()) {
			split.push_back(s.substr(start, i- start));
			start = i+1;
		}
	}
	return split;
}

/**
Finds a string in a vector of strings.

@param query The string to find.
@param v The vector to search in.
@return Index of query in v or -1 if query is not found in v.
*/
inline int findIndex(std::string query, std::vector<std::string> v) {
	for (size_t i = 0; i <= v.size(); i++)
		if (query.compare(v[i]) == 0)
			return i;
	return -1;
}

/*
Calculates the robust variance of E(G | D). var(x) = E(x^2) - E(x)^2

@param p Genotype frequency from EM algorithm
@return Robust variance of E(G | D)
*/
inline double calcRobustVar(std::vector<double> p) {
	return (4 * p[2] + p[1]) - pow(2 * p[2] + p[1], 2);
}


/*
Makes a copy of all SNPs

@param snps	SNPs to copy
@return vector of copied SNPs
*/
inline std::vector<SNP> cloneAll(std::vector<SNP> snps) {
	std::vector<SNP> clones;

	for (size_t i = 0; i < snps.size(); i++)
		clones.push_back(snps[i].clone());

	return clones;
}

//========================================================
// timing functions
//========================================================

#include <time.h>

inline clock_t startTime() {
	clock_t t;
	t = clock();
	return t;
}

inline void endTime(clock_t t, std::string label) {
	t = clock() - t;
	std::cout << "Timer " + label + " took ";
	std::cout << (float)t / CLOCKS_PER_SEC;
	std::cout << " seconds\n";
}

inline double endTime(clock_t t) {
	t = clock() - t;
	return (double)t / CLOCKS_PER_SEC;
}