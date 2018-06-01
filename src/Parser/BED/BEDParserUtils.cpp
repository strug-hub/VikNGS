#include "BEDParserUtils.h"

static const std::string BED_PARSER_UTILS = "BED parser utils";


BEDInterval getExons(std::vector<std::string> &lineSplit, BEDInterval &interval, int lineNumber) {

	int nexons;

	try {
		nexons = stoi(lineSplit[9]);
	}
	catch (...) {
		std::string message = "Line " + std::to_string(lineNumber) +
		" in BED file - Unexpected non-interger value in ninth column, skipping line.";
		printWarning(BED_PARSER_UTILS, message);
		return interval;
	}

	if (nexons < 1)
		return interval;

	std::vector<std::string> blockSizes = split(lineSplit[10], ',');
	std::vector<std::string> blockStarts = split(lineSplit[11], ',');

	if (blockSizes.size() < nexons) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file - Expected " + std::to_string(nexons);
		message += " exons in eleventh column but only found ";
		message += std::to_string(blockSizes.size()) + ". Skipping line.";
		printWarning(BED_PARSER_UTILS, message);
		return interval;
	}

	if (blockStarts.size() < nexons) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file - Expected " + std::to_string(nexons);
		message += " exons in twelfth column but only found ";
		message += std::to_string(blockStarts.size()) + ". Skipping line.";
		printWarning(BED_PARSER_UTILS, message);
		return interval;
	}

	int exonSize;
	int exonStart;

	for (int i = 0; i < nexons; i++) {

		try {
			exonSize = stoi(blockSizes[i]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber) +
			 " in BED file - Unexpected non-integer value in eleventh column, skipping line.";
			printWarning(BED_PARSER_UTILS, message, blockSizes[i]);
			return interval;
		}
		try {
			exonStart = stoi(blockStarts[i]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber);
			message += " in BED file - Unexpected non-integer value in twelfth column, skipping line.";
			printWarning(BED_PARSER_UTILS, message, blockStarts[i]);
			return interval;
		}

        interval.exonStart.emplace_back(interval.txStart + exonStart);
        interval.exonEnd.emplace_back(interval.txStart + exonStart + exonSize);
	}

	interval.valid = true;
	return interval;
}

BEDInterval lineToBedInterval(std::vector<std::string> lineSplit, int colsExpected, int lineNumber) {

    BEDInterval interval;

	if (lineSplit.size() < colsExpected) {

		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file - Expected at least " + std::to_string(colsExpected);
		message += " columns but only found " + std::to_string(lineSplit.size()) + ". Skipping line.";
		printWarning(BED_PARSER_UTILS, message);
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
		message += " in BED file - Unexpected non-integer value in second column, skipping line.";
		printWarning(BED_PARSER_UTILS, message, lineSplit[1]);
		return interval;
	}
	try {
		interval.txEnd = stoi(lineSplit[2]) + 1;
	}
	catch (...) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in BED file - Unexpected non-integer value in third column, skipping line.";
		printWarning(BED_PARSER_UTILS, message, lineSplit[2]);
		return interval;
	}

	if (colsExpected > 4) {
		try {
			interval.cdsStart = stoi(lineSplit[6]) + 1;
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber);
			message += " in BED file - Unexpected non-integer value in seventh column, skipping line.";
			printWarning(BED_PARSER_UTILS, message, lineSplit[6]);
			return interval;
		}
		try {
			interval.cdsEnd = stoi(lineSplit[7]) + 1;
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber);
			message += " in BED file - Unexpected non-integer value in eigth column, skipping line.";
			printWarning(BED_PARSER_UTILS, message, lineSplit[7]);
			return interval;
		}

	}

	if (colsExpected > 8) {

		interval = getExons(lineSplit, interval, lineNumber);
		if (!interval.valid)
			return interval;
	}

	interval.valid = true;
	return interval;
}

std::vector<Interval> toExonInterval(std::vector<BEDInterval> bed){

    int index = 0;
    std::vector<Interval> intervals;

    for(int i = 0; i < bed.size(); i++){
        for(int j = 0; j < bed[i].exonStart.size(); j++){

            Interval interval;
            interval.index = index;
            index++;
            interval.id = bed[i].id;
            interval.exon = j+1;

            interval.chr = bed[i].chr;
            interval.end = bed[i].exonEnd[j];
            interval.start = bed[i].exonStart[j];
            intervals.push_back(interval);
        }
    }

    return intervals;
}


std::vector<Interval> toInterval(std::vector<BEDInterval> bed){

    std::vector<Interval> intervals;

    for(int i = 0; i < bed.size(); i++){

        Interval interval;
        interval.index = i;
        interval.id = bed[i].id;
        interval.exon = -1;

        interval.chr = bed[i].chr;
        interval.end = bed[i].txEnd;
        interval.start = bed[i].txStart;
        intervals.push_back(interval);
    }

    return intervals;
}
