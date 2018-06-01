#pragma once
#include "MemoryMapped/MemoryMapped.h"
#include "../Log.h"
#include "../RVS.h"
#include "../Request.h"
#include "../TestInput.h"
#include "../Variant.h"
#include "../Math/MathHelper.h"
#include "BED/Interval.h"

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
std::vector<Variant> parseVCFLines(TestInput &input, Request &req);
std::map<std::string, int> getSampleIDMap(std::string vcfDir);

//SampleParser.cpp
void parseSampleLines(Request req, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup);

//BEDParser.cpp
std::vector<Interval> parseBEDLines(std::string bedDir, bool collapseExon);

//VariantFilter.cpp
int filterVariant(Request &req, Variant &variant, VectorXd &Y, std::string family);
void printFilterResults(Request &req, std::vector<std::string> variantInfo, std::vector<int> failCode, int total);

std::vector<Variant> removeDuplicates(std::vector<Variant> variants);
std::vector<Variant> filterMinorAlleleFrequency(std::vector<Variant> &variants, double mafCutoff, bool common);

struct File {
	MemoryMapped mmap;
    uint64_t pos;

    uint64_t segmentSize;
    uint64_t pageSize;
    uint64_t pagesPerSegment;
    uint64_t currentPage;

    bool lastSegment;
    uint64_t lineNumber;
    uint64_t memory = 4e9; //4gb

	inline void open(std::string directory) {
		try {

            pageSize = mmap.getpagesize();
            pagesPerSegment = memory/pageSize;
            uint64_t one = 1;
            pagesPerSegment = std::max(one, pagesPerSegment);

            mmap.open(directory, pagesPerSegment * pageSize, MemoryMapped::CacheHint::SequentialScan);
            currentPage = pagesPerSegment;

            pos = 0;
			lineNumber = 0;
            segmentSize = mmap.mappedSize();

            lastSegment = false;
            if(segmentSize >= mmap.size())
                lastSegment = true;

		}
		catch (...) {
			throwError("file struct", "Cannot open file from provided directory.", directory);
		}
	}

	inline void close() {

		mmap.close();
	}

	inline std::string nextLine() {

        std::string line = "";

		while (true) {
            if(pos >= segmentSize){

                if(lastSegment)
                    break;
                else
                    remap();
            }

            if (mmap[pos] == '\n')
				break;

            line+= mmap[pos];
            pos++;
		}

		lineNumber++;
        pos++;

        return line;
	}

	inline int getLineNumber() {
		return lineNumber;
	}

	inline bool hasNext() {
        return !lastSegment || (pos < segmentSize);
	}

    void remap() {

        printInfo("Reading another 4GB");
        size_t offset = currentPage * pageSize;
        size_t mapSize = pagesPerSegment * pageSize;

        mmap.remap(offset, mapSize);
        segmentSize = mmap.mappedSize();

        lastSegment = (offset + mapSize) >= mmap.size();
        currentPage += pagesPerSegment;
        pos=0;

    }
};
