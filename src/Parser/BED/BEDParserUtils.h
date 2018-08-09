#pragma once
#include "../Parser.h"

static const char BED_SEPARATOR = '\t';

/*
Takes a line from a BED file and extracts the exons from column 11 and 12. Number of
exons should match number in column 10.

@param lineSplit Line from BED file split by BED_SEPARATOR'.
@param interval Interval object; requires values for txStart and txEnd to be valid.
@param lineNumber Line number of BED file where lineSplit was derived.
@return Interval object with exons or Interval.valid = false if error was encountered.
*/
BEDInterval getExons(std::vector<std::string> &lineSplit, BEDInterval &interval, int lineNumber);


/*
Takes a line from a BED file and produces an Interval object from the line.

@param lineSplit Line from BED file split by tab '\t'.
@param colsExpected Minimum number of columns expected in each row of in BED file.
@param lineNumber Line number of BED file where lineSplit was derived.
@param collapseCoding Collapse along coding region (defined as thickStart and thickEnd in BED format).
@return Interval object with Interval.valid = false if error was encountered.
*/
BEDInterval lineToBedInterval(std::vector<std::string> lineSplit, int colsExpected, int lineNumber);


std::vector<Interval> toExonInterval(std::vector<BEDInterval> bed);
std::vector<Interval> toInterval(std::vector<BEDInterval> bed);

