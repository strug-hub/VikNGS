#include "stdafx.h"
#include "InputParser.h"

std::vector<double> handleCovariates(std::vector<std::string> cov) {
	std::vector<double> covariates;

	//todo: handle non-numeric covariates

	if (cov.size() > 0)
		for (size_t i = 0; i < cov.size(); i++)
			covariates.push_back(stod(cov[i]));


	return covariates;
}

//returns Y,G,H,Z sorted by 
void parseInfo(std::string sampleInfoDir, std::map<std::string, int> &IDmap,
	VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup) {
	
	MemoryMapped sampleInfo(sampleInfoDir);
	int pos = 0;
	int nID = IDmap.size();

	Y = VectorXd(nID);
	G = VectorXd(nID);

	bool initializeZ = false;

	//extract every sample from sampleInfo file, assumes one sample per line
	std::vector<std::string> lineSplit;
	int index;

	std::map<std::string, int> groupIDMap;
	int groupIndex = 0;

	while (true) {

		lineSplit = readLine(sampleInfo, pos);

		if (lineSplit.size() < 2)
			break;

		int index = IDmap[lineSplit[0]];
		Y[index] = std::stod(lineSplit[1]);

		std::string groupID = lineSplit[2];

		if (!groupIDMap.count(groupID)) {
			groupIDMap[groupID] = groupIndex;

			std::string depth = lineSplit[3];
		
			//TODO:check to see if readGroup is the consistent for every sample in a group
			if (depth.find("H") != std::string::npos ||
				depth.find("h") != std::string::npos)
				readGroup[groupIndex] = 1;
			else
				readGroup[groupIndex] = 0;

			groupIndex++;
		}
			
		G[index] = groupIDMap[groupID];

		std::vector<std::string> cov;
		for (size_t i = 4; i < lineSplit.size(); i++)
			cov.push_back(lineSplit[i]);

		std::vector<double> covariates = handleCovariates(cov);

		if (!initializeZ) {
			initializeZ = true;
			Z = MatrixXd(nID, 1 + covariates.size());
		}

		Z(index, 0) = 1;
		for (int i = 0; i < covariates.size(); i++)
			Z(index, i+1) = covariates[i];
	}

	//TODO: dummyproofing parsing here before finalized product
	/*
	if (lineCounter != countCase)
	std::cout << "Warning: " + std::to_string(lineCounter) + " lines were counted in case ID file but only " +
	std::to_string(countCase) + " were found to correspond to columns in .vcf file.\n";
	*/

	sampleInfo.close();
}
