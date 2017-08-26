#include "../RVS.h"
#include "InputParser.h"

struct Interval {
	std::string id;

	int start = 0;
	int end = 0;
	std::vector<int> indexes;
	std::string chr;

	inline void addIfIn(VCFLine variant, int index) {
		if (variant.loc >= start && variant.loc <= end && chr.compare(variant.chr) == 0) 
			indexes.push_back(index);
	}

	inline std::vector<int> getIndexes() {
		return indexes;
	}

	inline int nIndex() { return indexes.size(); }
};

/*
Parses intervals from a a BED file and collapses variants along intervals

@param collapse Vector of intervals on which to collapse variants.
@param variants Lines parsed from a VCF file.
@return the indexes of variants within each interval
*/
std::vector<std::vector<int>> collapseVariants(std::vector<Interval> &collapse, std::vector<VCFLine> variants) {

	for (int i = 0; i < variants.size(); i++)
		for (int j = 0; j < collapse.size(); j++)
				collapse[j].addIfIn(variants[i], i);

	std::vector<std::vector<int>> interval;

	for (int i = 0; i < collapse.size(); i++)
		if (collapse[i].nIndex() > 0)
			interval.push_back(collapse[i].getIndexes());

	return interval;
}

/*
Parses intervals from a a BED file and collapses variants along intervals

@param bedDir Directory to a bed file.
@param variants Lines parsed from a VCF file.
@return The indexes of variants within each interval.
*/
std::vector<std::vector<int>> parseBEDLines(std::string bedDir, std::vector<VCFLine> variants) {
	
	File bed;
	bed.open(bedDir);

	std::vector<Interval> collapse;

	while (bed.hasNext()) {

		std::vector<std::string> lineSplit = split(bed.nextLine(),  '\t');

		if (lineSplit.size() < 4)
			continue;

		Interval inv;
		
		inv.chr = lineSplit[0];
		inv.id = lineSplit[3];

		try {
			inv.start = stoi(lineSplit[1]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(bed.lineNumber);
			message += " in BED file: Unexpected value '" + lineSplit[1];
			message += "' in second column. Should be a numeric value. Skipping line.";
			printWarning(message);
			continue;
		}
		try {
			inv.end = stoi(lineSplit[2]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(bed.lineNumber);
			message += " in BED file: Unexpected value '" + lineSplit[2];
			message += "' in third column. Should be a numeric value. Skipping line.";
			printWarning(message);
			continue;
		}

		collapse.push_back(inv);
	}
	
	return collapseVariants(collapse, variants);
}
