#pragma once
#include "../InputParser.h"

static const char BED_SEPARATOR = '\t';

struct Interval {
	bool valid = false;
	std::string id;

	int txStart = -1;
	int txEnd = -1;
	int cdsStart = -1;
	int cdsEnd = -1;
	std::vector<int> exonStart;
	std::vector<int> exonEnd;

	std::vector<int> indexes;
	std::string chr;

    inline void txAddIfIn(Variant variant, int index) {
		if (variant.pos >= txStart && variant.pos <= txEnd && chr.compare(variant.chr) == 0)
			indexes.push_back(index);
	}

    inline void exonAddIfIn(Variant variant, int index) {
		if (variant.pos >= txStart && variant.pos <= txEnd && chr.compare(variant.chr) == 0) {
			for (int i = 0; i < exonStart.size(); i++) {
				if (variant.pos >= exonStart[i] && variant.pos <= exonEnd[i]) {
					indexes.push_back(index);
					return;
				}
			}
		}
	}

    inline void cdsAddIfIn(Variant variant, int index) {
		if (variant.pos >= cdsStart && variant.pos <= cdsEnd && chr.compare(variant.chr) == 0) {
			for (int i = 0; i < exonStart.size(); i++) {
				if (variant.pos >= exonStart[i] && variant.pos <= exonEnd[i]) {
					indexes.push_back(index);
					return;
				}
			}
		}
	}

	inline std::vector<int> getIndexes() {
		return indexes;
	}

	inline int nIndex() { return indexes.size(); }

	void print() {
		std::cout << "Interval " + id + " (" + std::to_string(indexes.size()) + " variants)";
		std::cout << "\n\tvalid = ";
		std::cout << valid;
		std::cout << "\n\tTranscription: ";
		std::cout << std::to_string(txStart) + " to " + std::to_string(txEnd);

		if (cdsStart > -1) {
			std::cout << "\n\tCoding: ";
			std::cout << std::to_string(cdsStart) + " to " + std::to_string(cdsEnd);
		}
		for (int i = 0; i < exonStart.size(); i++) {
			std::cout << "\n\tExon " + std::to_string(i) + ": ";
			std::cout << std::to_string(exonStart[i]) + " to " + std::to_string(exonEnd[i]);
		}
		std::cout << "\n";
	}

};

/*
Takes a line from a BED file and extracts the exons from column 11 and 12. Number of
exons should match number in column 10.

@param lineSplit Line from BED file split by BED_SEPARATOR'.
@param interval Interval object; requires values for txStart and txEnd to be valid.
@param lineNumber Line number of BED file where lineSplit was derived.
@return Interval object with exons or Interval.valid = false if error was encountered.
*/
Interval getExons(std::vector<std::string> &lineSplit, Interval &interval, int lineNumber);


/*
Parses intervals from a a BED file and collapses variants along intervals

@param collapse Vector of intervals on which to collapse variants.
@param variants Lines parsed from a VCF file.
@param collapseCoding Collapse along coding region (defined as thickStart and thickEnd in BED format).
@param collapseExon Collapse along exons only (defined as blocks in BED format).
@return the indexes of variants within each interval
*/
std::vector<std::vector<int>> collapseVariants(std::vector<Interval> &collapse, std::vector<Variant> variants, bool collapseCoding, bool collapseExon);


/*
Takes a line from a BED file and produces an Interval object from the line.

@param lineSplit Line from BED file split by tab '\t'.
@param colsExpected Minimum number of columns expected in each row of in BED file.
@param lineNumber Line number of BED file where lineSplit was derived.
@param collapseCoding Collapse along coding region (defined as thickStart and thickEnd in BED format).
@return Interval object with Interval.valid = false if error was encountered.
*/
Interval lineToInterval(std::vector<std::string> lineSplit, int colsExpected, int lineNumber, bool collapseCoding);
