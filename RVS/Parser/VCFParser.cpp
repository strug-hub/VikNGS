#include "InputParser.h"

/*
Extracts the last header row from a VCF file

@param vcf A VCF file.
@return The names of each column in the last header row from vcf.
*/
std::vector<std::string> extractHeader(File &vcf) {
	
	std::string line;

	while (true) {
		if (vcf.hasNext()) {
			line = vcf.nextLine();

			if (line.substr(0, 2) == "##")
				continue;
			else if (line.substr(0, 1) == "#")
				break;
			else {
				printError("In VCF file: Problem identifying header. Ensure header begins with '#'.");
				throw std::runtime_error("No VCF header");
			}
		}
	}
	
	return split(line, '\t');
}

/*
Calculates genotype likelihood from PL or GL for a single sample.

@param column VCF field that contains PL or GL value for one sample.
@param indexPL Index of PL values after column is split by ':'.
@param indexPL Index of GL values after column is split by ':'.
@param indexGT Index of GT values after column is split by ':'.
@return A GenotypeLikelihood object.
*/
GenotypeLikelihood getGenotypeLikelihood(std::string column, int indexPL, int indexGL, int indexGT) {

	GenotypeLikelihood gl;
	gl.L00 = NAN;
	gl.L01 = NAN;
	gl.L11 = NAN;
	gl.missing = true;

	std::vector<std::string> parts = split(column, ':');

	if (indexGT > -1) {
		std::string gt = parts[indexGT];

		if (gt.size() > 1) {

			if (gt[0] == '0' && gt[2] == '0') {
				gl.L00 = 1;
				gl.L01 = 0;
				gl.L11 = 0;
				gl.missing = false;
			}
			else if (gt[0] == '0' && gt[2] == '1') {
				gl.L00 = 0;
				gl.L01 = 1;
				gl.L11 = 0;
				gl.missing = false;
			}
			else if (gt[0] == '1' && gt[2] == '0') {
				gl.L00 = 0;
				gl.L01 = 1;
				gl.L11 = 0;
				gl.missing = false;
			}
			else if (gt[0] == '1' && gt[2] == '1') {
				gl.L00 = 0;
				gl.L01 = 0;
				gl.L11 = 1;
				gl.missing = false;
			}
		}
	}

	if (indexGL > -1) {
		std::vector<std::string> l = split(parts[indexGL], ',');

		if (l[0][0] == '.') {
			return gl;
		}

		if (l.size() < 3)
			throw std::runtime_error("Failed parsing GL");

		gl.L00 = pow(10, std::stod(l[0]));
		gl.L01 = pow(10, std::stod(l[1]));
		gl.L11 = pow(10, std::stod(l[2]));
		gl.missing = false;
	}

	else if (indexPL > -1) {
		std::vector<std::string> l = split(parts[indexPL], ',');

		if (l[0][0] == '.') {
			return gl;
		}

		if (l.size() < 3)
			throw std::runtime_error("Failed parsing PL");

		gl.L00 = pow(10, -std::stod(l[0])*0.1);
		gl.L01 = pow(10, -std::stod(l[1])*0.1);
		gl.L11 = pow(10, -std::stod(l[2])*0.1);
		gl.missing = false;
	}

	return gl;
}


/*
Extracts a single line (variant) from a VCF file.

@param columns The values in each column in a VCF file line.
@param lineNumber The number of the line in the VCF file.
@return A single VCFLine.
*/
VCFLine extractLine(std::vector<std::string> columns, int lineNumber) {

	//TODO:? get this data from reading header
	//assumes order of columns in VCF file
	int chr = 0;
	int loc = 1;
	int ref = 3;
	int alt = 4;
	int filter = 6;
	int format = 8;

	if (columns.size() < 8) {
		std::string message = "Line " + std::to_string(lineNumber);
		message += " in VCF file: Expected at least 8 columns but found " + std::to_string(columns.size()) + ".";
		printError(message);

		throw std::runtime_error("Columns missing");
	}

	VCFLine variant;

	variant.chr = columns[chr];
	variant.loc = std::stoi(columns[loc]);
	variant.ref = columns[ref];
	variant.alt = columns[alt];
	variant.filter = columns[filter];

	//finds the index of "PL" and "GL" from the FORMAT column
	std::string fmt = columns[format];

	int index = 0;
	int indexPL = -1;
	int indexGL = -1;
	int indexGT = -1;

	int pos = 0;
	while (pos < fmt.size()) {
		if (fmt[pos] == ':')
			index++;
		else if (fmt[pos] == 'P' && fmt[pos + 1] == 'L')
			indexPL = index;
		else if (fmt[pos] == 'G' && fmt[pos + 1] == 'L')
			indexGL = index;
		else if (fmt[pos] == 'G' && fmt[pos + 1] == 'T')
			indexGT = index;

		else if (fmt[pos] == '\t') {
			break;
		}
		pos++;
	}

	if (indexPL == -1 && indexGL == -1 && indexGT) {
		std::string message = "Line " + std::to_string(lineNumber) + " in VCF file: ";
		message += "Phred-scaled likelihoods (PL or GL) and genotype calls (GT) not found in FORMAT column. Skipping variant.";
		printWarning(message);

		variant.valid = false;
		return variant;
	}

	//get genotype likelihood for every sample
	for (int i = format + 1; i < columns.size(); i++) {
		try {
			GenotypeLikelihood gl = getGenotypeLikelihood(columns[i], indexPL, indexGL, indexGT);
			variant.likelihood.push_back(gl);
		}
		catch (...) {
			std::string message = "Line " + std::to_string(lineNumber) + " in VCF file: ";
			message += "Problem while trying to parse PL or GL on column " + std::to_string(i);
			message += ". Column value read from file: \"" + columns[i] + "\". Skipping variant.";
			printWarning(message);

			variant.valid = false;
			return variant;
		}
	}

	return variant;
}

/*
Reads every line (variant) from a VCF file.

@param vcfDir Directory of multisample VCF file.
@return A list of VCFLines parsed from the file at vcfDir.
*/
std::vector<VCFLine> parseVCFLines(std::string vcfDir) {

	File vcf;
	vcf.open(vcfDir);

	std::vector<VCFLine> variants;

	int pos = 0;

	try {
		//skips header
		extractHeader(vcf);
		
		while(vcf.hasNext()) {

			std::string line = vcf.nextLine();
			std::vector<std::string> columns = split(line, '\t');

			VCFLine variant = extractLine(columns, vcf.getLineNumber());
			
			if (vcf.getLineNumber() % 1000 == 0) {
				std::cout << vcf.getLineNumber();
				std::cout << '\n';
			}

			if(variant.isValid())
				variants.push_back(variant);

			if (vcf.getLineNumber() > 20000)
				break;
		}

		vcf.close();

	}
	catch (...) {
		std::string errorString = "VCF line ";
		errorString += std::to_string(vcf.getLineNumber());
		errorString += " caused an error while parsing.";
		printError(errorString);
		throw;
	}
	return variants;
}

/*
Reads every sample ID (columns after FORMAT) from a multisample VCF and stores it in a map.

@param vcfDir Directory of multisample VCF file.
@return A map from sample name to a unique integer.
*/
std::map<std::string, int> getSampleIDMap(std::string vcfDir) {

	//open VCF file and extract header
	std::vector<std::string> ID;
	
	try {
		File vcf;
		vcf.open(vcfDir);
		ID = extractHeader(vcf);
		vcf.close();
	}
	catch (...) {
		throw;
	}

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

	if (!flag) {
		printError("In VCF file: Cannot find FORMAT column in header!");
		throw std::runtime_error("Cannot find FORMAT column");
	}

	return IDmap;
}