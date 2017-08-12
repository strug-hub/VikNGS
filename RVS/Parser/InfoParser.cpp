#include "InputParser.h"

MatrixXd handleCovariate(std::vector<std::string> covariate, bool isNumeric, int columnNum) {

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
			std::string message = "Covariate in column ";
			message.append(std::to_string(columnNum + 5));
			message.append(" is being treated as a factor and has ");
			message.append(std::to_string(covIDMap.size()));
			message.append(" unique values.");
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
		std::cout << cov;

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

@return Vector that indicates which sample data was read correctly.
@effect Fills Y, Z, G and readGroup with values from input file
*/
//returns Y,G,H,Z sorted by 
bool parseInfo(std::string sampleInfoDir, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup, int highLowCutOff) {
	
	MemoryMapped sampleInfo(sampleInfoDir);
	int pos = 0;
	int nID = IDmap.size();

	Y = VectorXd(nID);
	G = VectorXd(nID);

	std::vector<std::vector<std::string>> covariates;


	//extract every sample from sampleInfo file, assumes one sample per line
	std::vector<std::string> lineSplit;
	int index;

	std::map<std::string, int> groupIDMap;
	int groupIndex = 0;
	int lineNumber = 0;

	while (true) {

		lineSplit = readLine(sampleInfo, pos);

		if (lineSplit.size() < 2)
			break;

		int index = IDmap[lineSplit[0]];
		
		try {
			Y[index] = std::stod(lineSplit[1]);
		}
		catch (...) {
			std::string message = "Unable to parse data file, sample ";
			message.append(lineSplit[0]);
			message.append(" (row ");
			message.append(std::to_string(lineNumber));
			message.append(") has unexpected value in second column: ");
			message.append(lineSplit[1]);
			printError(message);
			return false;
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
					std::string message = "Unable to parse data file, sample ";
					message.append(lineSplit[0]);
					message.append(" (row ");
					message.append(std::to_string(lineNumber));
					message.append(") has unexpected value in fourth column: ");
					message.append(depth);
					printError(message);
					return false;
				}
			}
			groupIndex++;
		}
			
		G[index] = groupIDMap[groupID];

		std::vector<std::string> cov;
		for (size_t i = 4; i < lineSplit.size(); i++)
			cov.push_back(lineSplit[i]);

		covariates.push_back(cov);

		lineNumber++;
	}

	Z = handleCovariates(covariates);


	//TODO: dummyproofing parsing here before finalized product
	/*
	if (lineCounter != countCase)
	std::cout << "Warning: " + std::to_string(lineCounter) + " lines were counted in case ID file but only " +
	std::to_string(countCase) + " were found to correspond to columns in .vcf file.\n";
	*/

	
	sampleInfo.close();
	return true;
}
