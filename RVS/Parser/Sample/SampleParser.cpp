#include "../InputParser.h"
#include "SampleParserUtils.h"

static const std::string SAMPLE_PARSER = "sample parser";

/*
Parses a data file containing tab-separated info for each sample.

@param sampleInfoDir Directory of tab-separated data file.
@param IDmap Map that gives a unique int for each sample name.
@param Y Vector for response variable (1st column).
@param Z Matrix for covariates (4th column onward).
@param G Vector identifing sample ID (2nd column).
@param readGroup Map of group ID to high or low read group.
@param highLowCutOff Cut-off for read depth. Anything >= value is considered high read depth.

@effect Fills Y, Z, G and readGroup with values from sample data file
*/
//returns Y,G,H,Z sorted by 
void parseSampleLines(Request req, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup) {
	
	File sampleInfo;
	sampleInfo.open(req.sampleDir);
	
	int nID = IDmap.size();

	Y = VectorXd(nID);
	G = VectorXd(nID);

	std::vector<std::vector<std::string>> covariates;

	std::string line;
	std::vector<std::string> lineSplit;
	int index;

	std::map<std::string, int> ID;
	std::map<std::string, int> groupIDMap;
	int groupIndex = 0;

	while (sampleInfo.hasNext()) {

		line = sampleInfo.nextLine();
		lineSplit = split(line, '\t');

		if (lineSplit.size() < 2)
			break;

		std::string sampleID = lineSplit[0];

		ID[sampleID] = 1;
		if (IDmap.count(sampleID) < 1) {
			std::string message = "Line " + std::to_string(sampleInfo.lineNumber) +
			" in sample information file - ID does not correspond to an ID in VCF file.";
			throwError(SAMPLE_PARSER, message, sampleID);
		}

		int index = IDmap[sampleID];
		
		try {
			Y[index] = std::stod(lineSplit[1]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(sampleInfo.lineNumber) +
			" in sample information file - Unexpected value non-numeric value in second column.";
			throwError(SAMPLE_PARSER, message, lineSplit[1]);
		}

		std::string groupID = lineSplit[2];

		if (!groupIDMap.count(groupID)) {
			groupIDMap[groupID] = groupIndex;

			std::string depth = lineSplit[3];
		
			if (depth.find("H") != std::string::npos ||
				depth.find("h") != std::string::npos)
				readGroup[groupIndex] = 1;
			else if (depth.find("L") != std::string::npos ||
				depth.find("l") != std::string::npos)
				readGroup[groupIndex] = 0;
			else {
				try {
					int d = std::stoi(depth);
					if(d >= req.highLowCutOff)
						readGroup[groupIndex] = 1;
					else
						readGroup[groupIndex] = 0;
				}
				catch (...) {
					std::string message = "Line " + std::to_string(sampleInfo.lineNumber);
					message += " in sample information file: Unexpected value in fourth column. Should be a numeric value or one of 'L', 'H'.";
					throwError(SAMPLE_PARSER, message, lineSplit[3]);
				}
			}
			groupIndex++;
		}
			
		G[index] = groupIDMap[groupID];

		std::vector<std::string> cov;
		for (size_t i = 4; i < lineSplit.size(); i++)
			cov.push_back(lineSplit[i]);

		covariates.push_back(cov);

	}

	if(covariates.size() > 0)
		Z = handleCovariates(covariates);

	sampleInfo.close();

	std::vector<std::string> missingSample;
	for (auto it = IDmap.begin(); it != IDmap.end(); ++it) {
		std::string id = it->first;
		
		if (ID[id] != 1)
			missingSample.push_back(id);
	}

	if (missingSample.size() > 0) {

		std::string message = "Samples were found in the provided VCF but are missing in the sample information file:";
		for (int i = 0; i < missingSample.size(); i++)
			message += " '" + missingSample[i] + "' ";

		throwError(SAMPLE_PARSER, message);
	}
}
