#pragma once
#include "MemoryMapped/MemoryMapped.h"
#include "../Log.h"
#include "../RVS.h"
#include "../Request.h"
#include "../TestInput.h"
#include "../Variant.h"
#include "../Math/MathHelper.h"

#include <iostream>  
#include <vector>
#include <string>
#include <algorithm>
#include <map>


//ParserTools.cpp
std::vector<std::vector<int>> collapseEveryK(int k, int n);
std::string extractString(MemoryMapped &charArray, int start, int end);
std::string trim(std::string str);
std::vector<std::string> split(std::string s, char sep);
std::vector<Variant> calculateExpectedGenotypes(std::vector<Variant> &variants);
VectorXd calcEG(std::vector<GenotypeLikelihood> &likelihood, VectorXd &p);
VectorXd calcEM(std::vector<GenotypeLikelihood> &likelihood);
std::string determineFamily(VectorXd Y);

//VCFParser.cpp
std::vector<Variant> parseVCFLines(std::string vcfDir, int minPos, int maxPos);
std::map<std::string, int> getSampleIDMap(std::string vcfDir);

//SampleParser.cpp
void parseSampleLines(Request req, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup);

//BEDParser.cpp
std::vector<std::vector<int>> parseBEDLines(std::string bedDir, std::vector<Variant> variants,
	bool collapseCoding, bool collapseExon);

//VariantFilter.cpp
std::vector<Variant> filterVariants(Request req, std::vector<Variant> &variants, VectorXd &Y, std::string family);
std::vector<Variant> removeDuplicates(std::vector<Variant> variants);
std::vector<Variant> filterMinorAlleleFrequency(std::vector<Variant> &variants, double mafCutoff, bool common);

struct File {
	MemoryMapped mmap;
	int pos;
	uint64_t maxPos;
	int lineNumber;


	inline void open(std::string directory) {
		try {
			mmap.open(directory);
			pos = 0;
			lineNumber = 0;
			maxPos = mmap.size();
		}
		catch (...) {
			throwError("file struct", "Cannot open file from provided directory.", directory);
		}
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
