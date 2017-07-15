#include "stdafx.h"
#include "../RVS.h"
#include "InputParser.h"

#include <iostream>  
#include <string>
#include <vector>

//todo
/*
std::vector<Interval> getIntervals(std::string bedDir) {

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
				chr = getString(bed, startPos, pos);
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
				start = stoi(getString(bed, startPos, pos));
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
				end = stoi(getString(bed, startPos, pos));
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
	return collapse;
}

void collapseVariants(std::vector<SNP> &snps, std::vector<Interval> &collapse) {

	for (size_t i = 0; i < snps.size(); i++)
		for (size_t j = 0; j < collapse.size(); j++)
			if (!isnan(snps[i].maf))
				collapse[j].addIfIn(snps[i], i);

	return;
}
*/