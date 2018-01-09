#include "BEDParserUtils.h"

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

	if (blockSizes.size() < nexons) {
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