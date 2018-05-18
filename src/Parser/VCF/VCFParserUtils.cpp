#include "VCFParserUtils.h"

static const std::string VCF_PARSER_UTILS = "VCF parser utils";

std::vector<std::string> extractHeader(File &vcf) {

	std::string line;

	while (vcf.hasNext()) {
		line = vcf.nextLine();

		if (line.substr(0, 2) == "##")
			continue;
		else if (line.substr(0, 1) == "#")
			break;
		else
			throwError(VCF_PARSER_UTILS, "Problem identifying header. Ensure header begins with a single '#'.");
	}

	return split(line, VCF_SEPARATOR);
}

GenotypeLikelihood getGT(std::string gt, GenotypeLikelihood gl) {

	double error1 = randomDouble(0.9995, 1.0);
	double error2 = randomDouble(0, (1 - error1));
	double p1 = error1;
	double p0_1 = error2;
	double p0_2 = 1 - (p1 + p0_1);

	if (gt[0] == '0' && gt[2] == '0') {
		gl.L00 = p1;
		gl.L01 = p0_1;
		gl.L11 = p0_2;
		gl.missing = false;
	}
	else if (gt[0] == '0' && gt[2] == '1') {
		gl.L00 = p0_1;
		gl.L01 = p1;
		gl.L11 = p0_2;
		gl.missing = false;
	}
	else if (gt[0] == '1' && gt[2] == '0') {
		gl.L00 = p0_1;
		gl.L01 = p1;
		gl.L11 = p0_2;
		gl.missing = false;
	}
	else if (gt[0] == '1' && gt[2] == '1') {
		gl.L00 = p0_1;
		gl.L01 = p0_2;
		gl.L11 = p1;
		gl.missing = false;
	}

	//todo : make optional?, missing = reference
	//else if (gt[0] == '.') {
	//	gl.L00 = p1;
	//	gl.L01 = p0_1;
	//	gl.L11 = p0_2;
	//	gl.missing = false;
	//}

	return gl;
}

GenotypeLikelihood getGL(std::vector<std::string> l, GenotypeLikelihood gl) {

	try {
		gl.L00 = pow(10, std::stod(l[0]));
		gl.L01 = pow(10, std::stod(l[1]));
		gl.L11 = pow(10, std::stod(l[2]));

		if(std::isnan(gl.L00) || std::isnan(gl.L01) || std::isnan(gl.L11))
			return gl;
		if(gl.L00 + gl.L01 + gl.L11 > 2.9)
		  return gl;
	}
	catch (...) { return gl; }

	gl.missing = false;
	return gl;
}

GenotypeLikelihood getPL(std::vector<std::string> l, GenotypeLikelihood gl) {

	try {
		gl.L00 = pow(10, -std::stod(l[0])*0.1);
		gl.L01 = pow(10, -std::stod(l[1])*0.1);
		gl.L11 = pow(10, -std::stod(l[2])*0.1);

		if(std::isnan(gl.L00) || std::isnan(gl.L01) || std::isnan(gl.L11))
			return gl;
		if(gl.L00 + gl.L01 + gl.L11 > 2.9)
		  return gl;
	}
	catch (...) { return gl; }

	gl.missing = false;
	return gl;
}

GenotypeLikelihood getGenotypeLikelihood(std::string column, int indexPL, int indexGL, int indexGT) {

	GenotypeLikelihood gl;
	gl.L00 = NAN;
	gl.L01 = NAN;
	gl.L11 = NAN;
	gl.missing = true;

	std::vector<std::string> parts = split(column, ':');

	//if GT is missing, return NAN
	if (indexGT > -1 && parts.size() > indexGT) {
		std::string gt = parts[indexGT];
		if (gt[0] == '.') {
			return gl;
			//todo print warning?
		}
	}

	//try to get GL
	if (indexGL > -1 && parts.size() > indexGL) {
		std::vector<std::string> l = split(parts[indexGL], ',');

		if (l[0][0] != '.' && l.size() == 3) {
			gl = getGL(l, gl);

			if (!gl.missing){
			  return gl;
			}
			//else
				//printWarning(VCF_PARSER_UTILS, "Expected GL at index " + std::to_string(indexGL) + " but parsing failed.", column);
		}
	}

	//try to get PL
	if (indexPL > -1 && parts.size() > indexPL) {
		std::vector<std::string> l = split(parts[indexPL], ',');

		if (l[0][0] != '.' && l.size() == 3) {
			gl = getPL(l, gl);

			if (!gl.missing){
			  return gl;
			}
			//else
			//	printWarning(VCF_PARSER_UTILS, "Expected PL at index " + std::to_string(indexPL) + " but parsing failed.", column);
		}
	}

	//try to get GT
	if (indexGT > -1 && parts.size() > indexGT) {
		std::string gt = parts[indexGT];
		gl = getGT(gt, gl);
	}

	
	return gl;
}

std::map<std::string, int> getSampleIDMap(std::string vcfDir) {

	//open VCF file and extract header
	std::vector<std::string> ID;

	File vcf;
	vcf.open(vcfDir);
	ID = extractHeader(vcf);
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

	if (!flag) 
		throwError(VCF_PARSER_UTILS, "Cannot find FORMAT column in header.");
	if (IDmap.size() <= 0)
		throwError(VCF_PARSER_UTILS, "No sample IDs were found in the VCF file header.");

	return IDmap;
}

Variant constructVariant(std::vector<std::string> columns) {

    Variant variant;

		if (columns.size() < 8) {
			variant.setInvalid("Expected at least 8 columns in VCF line but found " + std::to_string(columns.size()) + ".");
			return variant;
	}

	variant.chr = columns[CHR];
	variant.pos = std::stoi(columns[LOC]);
	variant.ref = columns[REF];
	variant.alt = columns[ALT];
	variant.filter = columns[FILTER];

	//finds the index of "PL" and "GL" from the FORMAT column
	std::string fmt = columns[FORMAT];

	int indexPL = -1;
	int indexGL = -1;
	int indexGT = -1;

	std::vector<std::string> format = split(fmt, ':');
	
	for (int i = 0; i < format.size(); i++) {

		if (format[i].size() == 2) {
			if (format[i][0] == 'P' && format[i][1] == 'L')
				indexPL = i;
			else if (format[i][0] == 'G' && format[i][1] == 'L')
				indexGL = i;
			else if (format[i][0] == 'G' && format[i][1] == 'T')
				indexGT = i;
		}
	}

	if (indexPL < 0 && indexGL < 0 && indexGT < 0) {
		variant.setInvalid("Phred-scaled likelihoods (PL or GL) and genotype calls (GT) not found in FORMAT column.");
		return variant;
	}

	//get genotype likelihood for every sample
	for (int i = FORMAT + 1; i < columns.size(); i++) {		
		GenotypeLikelihood gl = getGenotypeLikelihood(columns[i], indexPL, indexGL, indexGT);
		variant.likelihood.push_back(gl);
	}

	return variant;
}
