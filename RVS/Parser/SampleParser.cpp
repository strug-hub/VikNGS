#include "InputParser.h"

/*
Parses a data file containing tab-separated info for each sample.

@param covariate Column of covariates read from sample data file.
@param isNumeric True if covariate is a numeric value and not categorical.
@param columnIndex Index of the covariate column.
@return Matrix of covariate values
*/
MatrixXd handleCovariate(std::vector<std::string> covariate, bool isNumeric, int columnIndex) {

	if (isNumeric) {
		VectorXd z(covariate.size());
		for (int i = 0; i < covariate.size(); i++) {
		
			if (covariate[i] == "NA" || covariate[i] == "")
				z[i] = NAN;
			else
				z[i] = stod(covariate[i]);
		}
		return z;
	}
	else {

		std::map<std::string, int> covIDMap;
		int id = 0;

		for (int i = 0; i < covariate.size(); i++) {
			if (!covIDMap.count(covariate[i])) {
				covIDMap[covariate[i]] = id;
				id++;
			}
		}

		if (covIDMap.size() > 5) {
			std::string message = "In sample data file: Covariate in column ";
			message += std::to_string(columnIndex + 5) + " is being treated as categorical and has ";
			message += std::to_string(covIDMap.size()) + " unique values.";
			printWarning(message);
		}

		MatrixXd z(covariate.size(), covIDMap.size());
		z = z.setZero();

		int col;
		for (int i = 0; i < covariate.size(); i++) {

			if (covariate[i] == "NA" || covariate[i] == "") {
				for (int j = 0; j < covIDMap.size(); j++) 
					z(i, j) = NAN;
			}
			else {
				col = covIDMap[covariate[i]];
				z(i, col) = 1;
			}
		}
		return z;
	}
}

/*
Parses a data file containing tab-separated info for each sample.

@param covariates Column(s) of covariates read from sample data file.
@return Matrix of covariate values
*/
MatrixXd handleCovariates(std::vector<std::vector<std::string>> covariates) {
	MatrixXd Z;

	for (int i = 0; i < covariates[0].size(); i++) {

		std::vector<std::string> covariate;
		bool isNumeric = true;

		for (int j = 0; j < covariates.size(); j++) {
			covariate.push_back(covariates[j][i]);
			if (isNumeric) {
				try {
					std::stod(covariates[j][i]);
				}
				catch (...) {
					std::string cov = covariates[j][i];
					if(cov != "NA" && cov != "")
						isNumeric = false;
				}
			}
		}
		MatrixXd cov = handleCovariate(covariate, isNumeric, i);

		if (i == 0)
			Z = cov;
		else {
			MatrixXd temp = Z;
			MatrixXd Z(temp.rows(), temp.cols() + cov.cols());
			Z << temp, cov;
		}
	}

	return Z;
}

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
void parseSampleLines(std::string sampleInfoDir, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, int highLowCutOff) {
	
	File sampleInfo;
	sampleInfo.open(sampleInfoDir);
	
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
			std::string message = "Line " + std::to_string(sampleInfo.lineNumber);
			message += " in sample data file: ID '" + sampleID;
			message += "' does not correspond to an ID in VCF file.";
			printError(message);
			throw std::runtime_error("Sample data error");
		}

		int index = IDmap[sampleID];
		
		try {
			Y[index] = std::stod(lineSplit[1]);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(sampleInfo.lineNumber);
			message += " in sample data file: Unexpected value '" + lineSplit[1];
			message += "' in second column. Should be a numeric value.";
			printError(message);
			throw std::runtime_error("Sample data error");
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
					if(d >= highLowCutOff)
						readGroup[groupIndex] = 1;
					else
						readGroup[groupIndex] = 0;
				}
				catch (...) {
					std::string message = "Line " + std::to_string(sampleInfo.lineNumber);
					message += " in sample data file: Unexpected value '" + depth;
					message += "' in fourth column. Should be a numeric value or one of 'L', 'H'.";
					printError(message);
					throw std::runtime_error("Sample data error");
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

	Z = handleCovariates(covariates);

	sampleInfo.close();

	std::vector<std::string> missingSample;
	for (auto it = IDmap.begin(); it != IDmap.end(); ++it) {
		std::string id = it->first;
		
		if (ID[id] != 1)
			missingSample.push_back(id);
	}

	if (missingSample.size() > 0) {

		std::string message = "Samples were found in the provided VCF but are missing in the sample data file:";
		for (int i = 0; i < missingSample.size(); i++)
			message += " '" + missingSample[i] + "' ";

		printError(message);
		throw std::runtime_error("Sample data error");
	}
}
