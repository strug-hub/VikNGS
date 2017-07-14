#include "stdafx.h"
#include "InputParser.h"

/*
Extracts the last header row from a VCF file

@param vcf MemoryMapped object containing a VCF file.
@param pos Index of a character in MemoryMapped file.
@return The names of each column in the last header row from vcf.
@effect pos Is the index of the first character after the header upon return.
*/
std::vector<std::string> parseHeader(MemoryMapped &vcf, int &pos) {
	int lastPos = 0;

	//look for last header line
	//assumes header and only header begins with '#'
	while (true) {

		if (vcf[pos] == '\n') {
			if (vcf[pos + 1] != '#')
				break;
			lastPos = pos + 2;
		}
		pos++;
	}

	pos = lastPos;
	std::vector<std::string> colNames;

	//parse the header
	//assumes column names are tab-separated
	//assumes last line of header and first variant are separated by return character
	while (true) {
		if (vcf[pos] == '\t') {
			colNames.push_back(extractString(vcf, lastPos, pos));
			lastPos = pos + 1;
		}
		else if (vcf[pos] == '\n') {
			colNames.push_back(extractString(vcf, lastPos, pos));
			pos++;
			break;
		}
		pos++;
	}

	return colNames;
}

GenotypeLikelihood parsePL(int pos, MemoryMapped &vcf) {
	int startPos = pos;
	double l00;
	double l01;
	double l11;

	while (vcf[pos] != ',')
		pos++;
	pos++;
	l00 = stod(extractString(vcf, startPos, pos));
	startPos = pos;

	while (vcf[pos] != ',')
		pos++;
	pos++;
	l01 = stod(extractString(vcf, startPos, pos));
	startPos = pos;

	while (pos < vcf.mappedSize() && vcf[pos + 1] != '\t' && vcf[pos + 1] != ':' && vcf[pos + 1] != '\n')
		pos++;

	l11 = stod(extractString(vcf, startPos, pos + 1));

	//convert from Phred-scaled likelihoods
	GenotypeLikelihood gl;
	gl.L00 = pow(10, -l00*0.1);
	gl.L01 = pow(10, -l01*0.1);
	gl.L11 = pow(10, -l11*0.1);
	return gl;

}

int parseDP(int pos, MemoryMapped &vcf) {
	int startPos = pos;
	while (pos < vcf.mappedSize() && vcf[pos + 1] != '\t' && vcf[pos + 1] != ':' && vcf[pos + 1] != '\n')
		pos++;
	return stoi(extractString(vcf, startPos, pos + 1));
}

VCFLine extractLine(MemoryMapped &vcf, std::vector<int> &colPos) {

	//TODO: get this data from reading header!
	//assumes order of columns in VCF file
	int chr = 0;
	int loc = 1;
	int ref = 3;
	int alt = 4;
	int filter = 6;
	int format = 8;

	VCFLine variant;

	variant.chr = extractString(vcf, colPos[chr], colPos[chr + 1] - 1);
	variant.loc = std::stoi(extractString(vcf, colPos[loc], colPos[loc + 1] - 1));
	variant.ref = extractString(vcf, colPos[ref], colPos[ref + 1] - 1);
	variant.alt = extractString(vcf, colPos[alt], colPos[alt + 1] - 1);
	variant.filter = extractString(vcf, colPos[filter], colPos[filter + 1] - 1);

	//finds the index of "PL" and "DP" from the FORMAT column
	//assumes FORMAT column is the 9th column

	int index = 0;
	int indexPL = -1;
	int indexDP = -1;

	int pos = colPos[format];
	while (true) {
		if (vcf[pos] == ':')
			index++;
		else if (vcf[pos] == 'P' && vcf[pos + 1] == 'L')
			indexPL = index;
		else if (vcf[pos] == 'D' && vcf[pos + 1] == 'P')
			indexDP = index;
		else if (vcf[pos] == '\t') {
			if (indexPL == -1) {
				std::cout << "ERROR: Phred-scaled likelihoods (PL) not found in FORMAT column ";
				std::cout << "or FORMAT is not the 9th column in the .vcf file.";
			}
			if (indexDP == -1) {
				std::cout << "ERROR: Read depth (DP) not found in FORMAT column ";
				std::cout << "or FORMAT is not the 9th column in the .vcf file.";
			}

			break;
		}
		pos++;
	}

	std::vector<GenotypeLikelihood> likelihood;
	std::vector<int> readDepth;
	int counter;


	//get genotype likelihood for every sample
	for (int i = format + 1; i < colPos.size(); i++) {
		pos = colPos[i];
		counter = 0;

		if (vcf[pos] == '.') {
			GenotypeLikelihood gl;
			gl.L00 = NAN;
			gl.L01 = NAN;
			gl.L11 = NAN;
			gl.missing = true;
			likelihood.push_back(gl);
			readDepth.push_back(0);

			continue;
		}

		while (true) {

			if (pos >= vcf.mappedSize() || vcf[pos] == '\t' || vcf[pos] == '\n')
				break;

			if (vcf[pos] == ':') {
				counter++;
				if (counter == indexPL) {
					pos++;
					likelihood.push_back(parsePL(pos, vcf));
				}

				if (counter == indexDP) {
					pos++;
					readDepth.push_back(parseDP(pos, vcf));
				}

			}
			pos++;
		}

	}

	variant.likelihood = likelihood;
	variant.readDepth = readDepth;

	return variant;
}

std::vector<VCFLine> parseVCFLines(std::string vcfDir) {

	MemoryMapped vcf(vcfDir);
	std::vector<VCFLine> variants;

	int pos = 0;

	//skips header
	std::vector<std::string> colNames = parseHeader(vcf, pos);

	int end = (int)vcf.mappedSize();

	while (pos < end) {
		std::vector<int> colPos;
		colPos.push_back(pos);

		while (true) {

			//parsing the columns before ID columns
			if (vcf[pos] == '\t')
				colPos.push_back(pos + 1);

			pos++;
			if (pos >= vcf.mappedSize() || vcf[pos] == '\n')
				break;

		}

		if (colPos.size() > 1)
			variants.push_back(extractLine(vcf, colPos));

		pos++;
	}

	vcf.close();
	return variants;
}

std::map<std::string, int> getSampleIDMap(std::string vcfDir) {

	//open VCF file and find the line with column names
	MemoryMapped vcf(vcfDir);
	int pos = 0;
	std::vector<std::string> ID = parseHeader(vcf, pos);
	vcf.close();

	bool flag = false;
	int count = 0;
	std::map<std::string, int> IDmap;

	for (int i = 0; i < ID.size(); i++) {
		if (flag) {
			IDmap[ID[i]] = count;
			count++;
		}
		else if (ID[i] == "FORMAT")
			flag = true;
	}

	return IDmap;
}