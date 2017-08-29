#include "../RVS.h"
#include "InputParser.h"

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

	inline void txAddIfIn(VCFLine variant, int index) {
		if (variant.loc >= txStart && variant.loc <= txEnd && chr.compare(variant.chr) == 0) 
			indexes.push_back(index);
	}

	inline void exonAddIfIn(VCFLine variant, int index) {
		if (variant.loc >= txStart && variant.loc <= txEnd && chr.compare(variant.chr) == 0) {
			for (int i = 0; i < exonStart.size(); i++) {
				if (variant.loc >= exonStart[i] && variant.loc <= exonEnd[i]) {
					indexes.push_back(index);
					return;
				}
			}
		}
	}

	inline void cdsAddIfIn(VCFLine variant, int index) {
		if (variant.loc >= cdsStart && variant.loc <= cdsEnd && chr.compare(variant.chr) == 0) {
			for (int i = 0; i < exonStart.size(); i++) {
				if (variant.loc >= exonStart[i] && variant.loc <= exonEnd[i]) {
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
Parses intervals from a a BED file and collapses variants along intervals

@param collapse Vector of intervals on which to collapse variants.
@param variants Lines parsed from a VCF file.
@param collapseCoding Collapse along coding region (defined as thickStart and thickEnd in BED format).
@param collapseExon Collapse along exons only (defined as blocks in BED format).
@return the indexes of variants within each interval
*/
std::vector<std::vector<int>> collapseVariants(std::vector<Interval> &collapse, std::vector<VCFLine> variants, bool collapseCoding, bool collapseExon) {

	if (!collapseCoding && !collapseExon) {
		for (int i = 0; i < variants.size(); i++)
			for (int j = 0; j < collapse.size(); j++)
				collapse[j].txAddIfIn(variants[i], i);
	}
	else if (collapseExon) {
		for (int i = 0; i < variants.size(); i++)
			for (int j = 0; j < collapse.size(); j++)
				collapse[j].exonAddIfIn(variants[i], i);
	}
	else if (collapseCoding) {
		for (int i = 0; i < variants.size(); i++)
			for (int j = 0; j < collapse.size(); j++)
				collapse[j].cdsAddIfIn(variants[i], i);
	}
	
	std::vector<std::vector<int>> interval;

	for (int i = 0; i < collapse.size(); i++) {
		//for debugging:
		//collapse[i].print();
		if (collapse[i].nIndex() > 0)
			interval.push_back(collapse[i].getIndexes());
	}
	
	return interval;
}

/*
Takes a line from a BED file and extracts the exons from column 11 and 12. Number of
exons should match number in column 10.

@param lineSplit Line from BED file split by tab '\t'.
@param interval Interval object; requires values for txStart and txEnd to be valid.
@param lineNumber Line number of BED file where lineSplit was derived.
@return Interval object with exons or Interval.valid = false if error was encountered.
*/
Interval getExons(std::vector<std::string> &lineSplit, Interval &interval, int lineNumber) {

	int nexons;

	try {
		nexons = stoi(lineSplit[9]);
	}
	catch (...) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file: Unexpected value '" + lineSplit[9];
		message += "' in ninth column. Should be an integer value. Skipping line.";
		printWarning(message);
		return interval;
	}

	if (nexons < 1)
		return interval;

	std::vector<std::string> blockSizes = split(lineSplit[10], ',');
	std::vector<std::string> blockStarts = split(lineSplit[11], ',');

	if(blockSizes.size() < nexons) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file: Expected " + std::to_string(nexons);
		message += " exons in eleventh column but only found ";
		message += std::to_string(blockSizes.size()) + ". Skipping line.";
		printWarning(message);
		return interval;
	}

	if (blockStarts.size() < nexons) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file: Expected " + std::to_string(nexons);
		message += " exons in twelfth column but only found ";
		message += std::to_string(blockStarts.size()) + ". Skipping line.";
		printWarning(message);
		return interval;
	}

	int exonSize;
	int exonStart;

	for (int i = 0; i < nexons; i++) {

		try {
			exonSize = stoi(blockSizes[i]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber);
			message += " in BED file: Unexpected value '" + blockSizes[i];
			message += "' in eleventh column. Should be an integer value. Skipping line.";
			printWarning(message);
			return interval;
		}
		try {
			exonStart = stoi(blockStarts[i]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber);
			message += " in BED file: Unexpected value '" + blockStarts[i];
			message += "' in twelfth column. Should be an integer value. Skipping line.";
			printWarning(message);
			return interval;
		}

		interval.exonStart.push_back(interval.txStart + exonStart);
		interval.exonEnd.push_back(interval.txStart + exonStart + exonSize);
	}

	interval.valid = true;
	return interval;
}


/*
Takes a line from a BED file and produces an Interval object from the line.

@param lineSplit Line from BED file split by tab '\t'.
@param colsExpected Minimum number of columns expected in each row of in BED file.
@param lineNumber Line number of BED file where lineSplit was derived.
@param collapseCoding Collapse along coding region (defined as thickStart and thickEnd in BED format).
@return Interval object with Interval.valid = false if error was encountered.
*/
Interval lineToInterval(std::vector<std::string> lineSplit, int colsExpected, int lineNumber, bool collapseCoding) {

	Interval interval;

	if (lineSplit.size() < colsExpected) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file: Expected at least " + std::to_string(colsExpected);
		message += " but only found " + std::to_string(lineSplit.size()) + ". Skipping line.";
		printWarning(message);
		return interval;
	}

	if (lineSplit.size() < 4)
		interval.id = std::to_string(lineNumber);
	else
		interval.id = lineSplit[3];

	//BED files start at 0, VCF start at 1
	try {
		interval.txStart = stoi(lineSplit[1]) + 1;
	}
	catch (...) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file: Unexpected value '" + lineSplit[1];
		message += "' in second column. Should be an integer value. Skipping line.";
		printWarning(message);
		return interval;
	}
	try {
		interval.txEnd = stoi(lineSplit[2]) + 1;
	}
	catch (...) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file: Unexpected value '" + lineSplit[2];
		message += "' in third column. Should be an integer value. Skipping line.";
		printWarning(message);
		return interval;
	}

	if (colsExpected > 4) {
		try {
			interval.cdsStart = stoi(lineSplit[6]) + 1;
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber);
			message += " in BED file: Unexpected value '" + lineSplit[6];
			message += "' in seventh column. Should be an integer value. Skipping line.";
			printWarning(message);
			return interval;
		}
		try {
			interval.cdsEnd = stoi(lineSplit[7]) + 1;
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber);
			message += " in BED file: Unexpected value '" + lineSplit[7];
			message += "' in eigth column. Should be an integer. Skipping line.";
			printWarning(message);
			return interval;
		}

		//non-protein coding gene, don't consider (ex. non-coding RNA genes)
		if (collapseCoding && interval.cdsEnd - interval.cdsStart < 1)
			return interval;
	}

	if (colsExpected > 8) {

		interval = getExons(lineSplit, interval, lineNumber);
		if (!interval.valid)
			return interval;
	}

	interval.valid = true;
	return interval;
}


/*
Parses intervals from a a BED file and collapses variants along intervals. Variants can be collapsed
along exons, exons+introns or any arbitrary region.
Note: BED files start position on a chromosome is 0, VCF start at 1

@param bedDir Directory to a bed file.
@param variants Lines parsed from a VCF file.
@param collapseCoding Collapse along coding region (defined as thickStart and thickEnd in BED format).
@param collapseExon Collapse along exons only (defined as blocks in BED format).
@return The indexes of variants within each interval.
*/
std::vector<std::vector<int>> parseBEDLines(std::string bedDir, std::vector<VCFLine> variants, bool collapseCoding, bool collapseExon) {
	
	File bed;
	bed.open(bedDir);

	std::vector<Interval> collapse;

	int colsExpected = 3;
	if (collapseCoding || collapseExon)
		colsExpected = 12;

	while (bed.hasNext()) {

		std::vector<std::string> lineSplit = split(bed.nextLine(),  '\t');

		Interval interval = lineToInterval(lineSplit, colsExpected, bed.getLineNumber(), collapseCoding);
		
		if(interval.valid)
			collapse.push_back(interval);
	}
	
	return collapseVariants(collapse, variants, collapseCoding, collapseExon);
}
