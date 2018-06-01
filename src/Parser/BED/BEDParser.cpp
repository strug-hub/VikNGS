#include "BEDParserUtils.h"

/*
Parses intervals from a a BED file.
Note: BED files start position on a chromosome is 0, VCF start at 1

@param bedDir Directory to a bed file.
@param collapseExon Collapse along exons only (defined as blocks in BED format).
@return The indexes of variants within each interval.
*/
std::vector<Interval> parseBEDLines(std::string bedDir, bool collapseExon) {
	
	File bed;
	bed.open(bedDir);

    std::vector<BEDInterval> intervals;

	int colsExpected = 3;
    if (collapseExon)
		colsExpected = 12;

	while (bed.hasNext()) {

        std::string line = bed.nextLine();
        if (line[0] == '#')
            continue;

        std::vector<std::string> lineSplit = split(line, BED_SEPARATOR);

        BEDInterval interval = lineToBedInterval(lineSplit, colsExpected, bed.getLineNumber());
		
		if(interval.valid)
            intervals.push_back(interval);
	}
	
    if(collapseExon)
        return toExonInterval(intervals);
    else
        return toInterval(intervals);
}
