#include "../RVS.h"
#include "InputParser.h"

struct Interval {
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

std::vector<std::vector<int>> collapseVariants(std::vector<Interval> &collapse, std::vector<VCFLine> variants) {

	for (int i = 0; i < variants.size(); i++)
		for (int j = 0; j < collapse.size(); j++)
				collapse[j].addIfIn(variants[i], i);

	std::vector<std::vector<int>> interval;

	for (int i = 0; i < collapse.size(); i++)
		if (collapse[i].nIndex() > 0)
			interval.push_back(collapse[i].getIndexes());

	//todo: do we want to keep intervals that are not in the bed file???

	return interval;
}

std::vector<std::vector<int>> parseIntervals(std::string bedDir, std::vector<VCFLine> variants) {
	MemoryMapped bed(bedDir);
	int pos = 0;
	std::vector<Interval> collapse;

	//loops once per line
	while (true) {
		std::string chr = "";
		int start = 0;
		int end = 0;

		int startPos = pos;
		bool breakFlag = false;

		//chr column
		while (true) {
			pos++;
			if (pos >= bed.mappedSize()) {
				breakFlag = true;
				break;
			}
			if (bed[pos] == '\t') {
				chr = extractString(bed, startPos, pos);
				pos++;
				break;
			}
		}

		startPos = pos;
		if (breakFlag)
			break;

		//start column
		while (true) {
			pos++;
			if (pos >= bed.mappedSize()) {
				breakFlag = true;
				break;
			}
			if (bed[pos] == '\t') {
				start = stoi(extractString(bed, startPos, pos));
				pos++;
				break;
			}
		}

		startPos = pos;
		if (breakFlag)
			break;

		//end column
		while (true) {
			pos++;
			if (pos >= bed.mappedSize()) {
				breakFlag = true;
				break;
			}
			if (bed[pos] == '\t') {
				end = stoi(extractString(bed, startPos, pos));
				pos++;
				break;
			}
		}

		startPos = pos;
		if (breakFlag)
			break;

		//loop till next line column
		while (true) {
			pos++;
			if (pos >= bed.mappedSize()) {
				breakFlag = true;
				break;
			}
			if (bed[pos] == '\n') {
				pos++;
				if (pos >= bed.mappedSize())
					breakFlag = true;
				break;
			}
		}

		Interval inv;
		inv.chr = chr;
		inv.start = start;
		inv.end = end;

		collapse.push_back(inv);

		if (breakFlag)
			break;
	}
	return collapseVariants(collapse, variants);
}
