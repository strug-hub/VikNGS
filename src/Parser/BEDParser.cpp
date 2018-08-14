#include "Parser.h"

static const std::string ERROR_SOURCE = "BED_PARSER";

/**
Parses intervals from a a BED file.
Note: BED files start position on a chromosome is 0, VCF start at 1

@param bedDir Directory to a bed file.
@param collapse Collapsing method.
@return An interval set, mapping chromosomes to sorted intervals.
*/
IntervalSet parseBEDLines(std::string bedDir, CollapseType collapse) {

    IntervalSet is;
    if(collapse == CollapseType::COLLAPSE_K || collapse == CollapseType::NONE)
        return is;

	File bed;
	bed.open(bedDir);

	while (bed.hasNext()) {

        std::string line = bed.nextLine();
        if (line[0] == '#')
            continue;

        std::vector<std::string> lineSplit = splitString(line, BED_SEP);

        if(collapse == CollapseType::COLLAPSE_GENE){
            Interval interval = lineToGene(lineSplit, bed.getLineNumber());
            is.addInterval(interval);
        }
        if(collapse == CollapseType::COLLAPSE_EXON){
            std::vector<Interval> invs = lineToExons(lineSplit, bed.getLineNumber());
            for(Interval inv : invs)
                is.addInterval(inv);
        }
	}

    is.sort();
    return is;
}

/**
Takes a line from a BED file and produces an Interval object from the line.

@param lineSplit Line from BED file split into a list.
@param lineNumber Line number of BED file where lineSplit was derived.
@return Interval object corresponding to gene.
*/
Interval lineToGene(std::vector<std::string> lineSplit, int lineNumber){

    Interval interval;
    //3 lines expected
    if (lineSplit.size() < 3) {
        std::string message = "Line " + std::to_string(lineNumber);
        message += " in BED file - Expected at least 3 columns but only found ";
        message += std::to_string(lineSplit.size()) + ". Skipping line.";
        printWarning(ERROR_SOURCE, message);
        return interval;
    }

    interval.chr = lineSplit[0];
    if (lineSplit.size() < 4)
        interval.id = std::to_string(lineNumber);
    else
        interval.id = lineSplit[3];

    //BED files start at 0, VCF start at 1
    interval.start = stoi(lineSplit[1]) + 1;
    interval.end = stoi(lineSplit[2]) + 1;

    return interval;
}

/**
Takes a line from a BED file and produces an Interval object from the line exons.

@param lineSplit Line from BED file split by tab '\t'.
@param lineNumber Line number of BED file where lineSplit was derived.
@return Vector of interval objects corresponding to exons.
*/
std::vector<Interval> lineToExons(std::vector<std::string> lineSplit, int lineNumber){

    std::vector<Interval> intervals;
    //12 lines expected
    if (lineSplit.size() < 12) {
        std::string message = "Line " + std::to_string(lineNumber);
        message += " in BED file - Expected at least 12 columns but only found ";
        message += std::to_string(lineSplit.size()) + ". Skipping line.";
        printWarning(ERROR_SOURCE, message);
        return intervals;
    }

    int nexons = stoi(lineSplit[9]);
    if (nexons < 1)
        return intervals;

    //BED files start at 0, VCF start at 1
    int txStart = stoi(lineSplit[1]) + 1;

    std::vector<std::string> blockSizes = splitString(lineSplit[10], ',');
    std::vector<std::string> blockStarts = splitString(lineSplit[11], ',');

    for (size_t i = 0; i < nexons; i++) {

        Interval inv;

        inv.chr = lineSplit[0];
        inv.id = lineSplit[3] + "_" + std::to_string(i);
        int exonSize = stoi(blockSizes[i]);
        int exonStart = stoi(blockStarts[i]);

        inv.start = txStart + exonStart;
        inv.end = txStart + exonStart + exonSize;
        intervals.push_back(inv);
    }

    return intervals;
}
