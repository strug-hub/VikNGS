#include "BEDParserUtils.h"

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

		std::vector<std::string> lineSplit = split(bed.nextLine(), BED_SEPARATOR);

		Interval interval = lineToInterval(lineSplit, colsExpected, bed.getLineNumber(), collapseCoding);
		
		if(interval.valid)
			collapse.push_back(interval);
	}
	
	return collapseVariants(collapse, variants, collapseCoding, collapseExon);
}