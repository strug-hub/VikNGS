#include "Parser.h"
#include "../Math/Math.h"
#include "File.h"
#include "../Variant.h"
#include "../Log.h"
static const std::string ERROR_SOURCE = "VCF_PARSER";


/**
Extracts the last header row from a VCF file

@param vcf A VCF file.
@throws Error if cannot find a header in file.
@return Vector of names of each column in the last header row from VCF.
*/
std::vector<std::string> extractHeader(File &vcf) {
    std::string header = extractHeaderLine(vcf);
    return splitString(header, VCF_SEP);

}

std::string extractHeaderLine(File &vcf) {

    std::string line;

    while (vcf.hasNext()) {
        line = vcf.nextLine();

        if (line.substr(0, 2) == "##")
            continue;
        else if (line.substr(0, 1) == "#")
            break;
        else
            throwError(ERROR_SOURCE, "Problem identifying header. Ensure column header begins with a single '#'.");

    }

    return line;
}

/**
Reads every sample ID (columns after FORMAT) from a multisample VCF and stores it in a map.

@param vcfDir Directory of multisample VCF file.
@return A map from sample name to a unique integer.
*/
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

    for (size_t i = 0; i < ID.size(); i++) {
        if (flag) {
            IDmap[ID[i]] = count;
            count++;
        }
        else if (ID[i] == "FORMAT")
            flag = true;
    }

    if (!flag)
        throwError(ERROR_SOURCE, "Cannot find FORMAT column in VCF header.");
    if (IDmap.size() <= 0)
        throwError(ERROR_SOURCE, "No sample IDs were found in the VCF file header.");

    return IDmap;
}

/**
Builds variant from a VCF line.

@param columns VCF file line split into columns.
@param getLikelihoods Extract genotype likelihoods from PL/GL.
@param calculateCalls Produce genotype calls from genotype likelihood.
@param getVCFCall Extract genotype calls from GT.

@return A Variant object corresponding to VCF line.
*/
Variant constructVariant(std::vector<std::string> &columns, bool calculateExpected, bool calculateCalls, bool getVCFCalls){

    if (columns.size() < 8) {
        printWarning(ERROR_SOURCE, "Found a variant with " + std::to_string(columns.size()) +
                     " columns (at least 8 expected). Skipping variant.");
        return Variant();
    }

    //finds the index of "PL" and "GL" from the FORMAT column
    std::string fmt = columns[FORMAT];

    int indexPL = -1;
    int indexGL = -1;
    int indexGT = -1;

    std::vector<std::string> format = splitString(fmt, ':');

    for (size_t i = 0; i < format.size(); i++) {

        if (format[i].size() == 2) {
            if (format[i][0] == 'P' && format[i][1] == 'L')
                indexPL = static_cast<int>(i);
            else if (format[i][0] == 'G' && format[i][1] == 'L')
                indexGL = static_cast<int>(i);
            else if (format[i][0] == 'G' && format[i][1] == 'T')
                indexGT = static_cast<int>(i);
        }
    }

    if (indexPL < 0 && indexGL < 0 && indexGT < 0) {
        printWarning(ERROR_SOURCE, "Genotype likelihoods (PL or GL) and genotype calls (GT) not found in FORMAT column. Skipping variant.");
        return Variant();
    }

    try{
        Variant variant(columns[CHROM], std::stoi(columns[POS]), columns[ID], columns[REF], columns[ALT]);

        size_t nsamp = columns.size() - (FORMAT + 1);

        if(calculateExpected || calculateCalls){
            std::vector<Vector3d> likelihoods;
            likelihoods.reserve(nsamp);
            for (size_t i = FORMAT + 1; i < columns.size(); i++)
                likelihoods.emplace_back(getGenotypeLikelihood(columns[i], indexPL, indexGL, indexGT));

            if(calculateExpected)
                variant.setExpectedGenotypes(likelihoods);
            if(calculateCalls)
                variant.setCallGenotypes(likelihoods);
        }
        if(getVCFCalls){
            VectorXd calls(nsamp);
            int index = 0;
            for (size_t i = FORMAT + 1; i < columns.size(); i++){
                calls[index] = getVCFGenotypeCall(columns[i], indexGT);
                index++;
            }

            variant.setVCFCallGenotypes(calls);
        }
        return variant;

    }catch(...){
        printWarning(ERROR_SOURCE, "Issue when trying to parse variant " +
                     columns[CHROM] + " " + columns[POS] + ". Skipping variant.");
        return Variant();
    }
}

Vector3d getGT(std::string gt) {

	double error1 = randomDouble(0.9995, 1.0);
	double error2 = randomDouble(0, (1 - error1));
	double p1 = error1;
	double p0_1 = error2;
	double p0_2 = 1 - (p1 + p0_1);
    Vector3d gl;

    if (gt[0] == '0'){
        if(gt[2] == '0') {
            gl[0] = p1;
            gl[1] = p0_1;
            gl[2] = p0_2;
            return gl;
        }
        else if (gt[2] == '1'){
            gl[0] = p0_1;
            gl[1] = p1;
            gl[2] = p0_2;
            return gl;
        }
    }
    else if (gt[0] == '1'){
        if(gt[2] == '0') {
            gl[0] = p0_2;
            gl[1] = p1;
            gl[2] = p0_1;
            return gl;
        }
        else if (gt[2] == '1'){
            gl[0] = p0_1;
            gl[1] = p0_2;
            gl[2] = p1;
            return gl;
        }
    }

	//todo : make optional?, missing = reference
	//else if (gt[0] == '.') {
	//	gl.L00 = p1;
	//	gl.L01 = p0_1;
	//	gl.L11 = p0_2;
	//	gl.missing = false;
	//}

    gl[0] = NAN;
    gl[1] = NAN;
    gl[2] = NAN;
	return gl;
}

/**
Calculates genotype likelihood from PL or GL (or GT) for a single sample.

@param column VCF field that contains PL/GL/GT values for one sample.
@param indexPL Index of PL values after column is split by ':'.
@param indexPL Index of GL values after column is split by ':'.
@param indexGT Index of GT values after column is split by ':'.
@return A vector with the 3 genotype likelihoods. Vector of NAN if issue in parsing.
*/
Vector3d getGenotypeLikelihood(std::string &column, int indexPL, int indexGL, int indexGT) {

    Vector3d gl;
    gl[0] = NAN;
    gl[1] = NAN;
    gl[2] = NAN;

    std::vector<std::string> split = splitString(column, ':');
    int size = static_cast<int>(split.size());

	//if GT is missing, return NAN
    if (indexGT > -1 && size > indexGT) {
        std::string gt = split[static_cast<size_t>(indexGT)];

        if (gt[0] == '.')
			return gl;
	}

	//try to get GL
    if (indexGL > -1 && size > indexGL) {
        std::vector<std::string> l = splitString(split[static_cast<size_t>(indexGL)], ',');

        if (l.size() == 3 && l[0][0] != '.') {
            try {

                if(l[0] == "-nan" || l[0] == "-1.4013e-45") l[0] = "-100";
                if(l[1] == "-nan" || l[1] == "-1.4013e-45") l[1] = "-100";
                if(l[2] == "-nan" || l[2] == "-1.4013e-45") l[2] = "-100";

                gl[0] = pow(10, std::stod(l[0]));
                gl[1] = pow(10, std::stod(l[1]));
                gl[2] = pow(10, std::stod(l[2]));

                if(gl.sum() > 0)
                    return gl;
            }
            catch (...) {
                gl[0] = NAN;
                gl[1] = NAN;
                gl[2] = NAN;
            }
        }
	}

	//try to get PL
    if (indexPL > -1 && size > indexPL) {
        std::vector<std::string> l = splitString(split[static_cast<size_t>(indexPL)], ',');

        if (l.size() == 3 && l[0][0] != '.') {
            try {
                if(l[0] == "-nan" || l[0] == "-1.4013e-45") l[0] = "10";
                if(l[1] == "-nan" || l[1] == "-1.4013e-45") l[1] = "10";
                if(l[2] == "-nan" || l[2] == "-1.4013e-45") l[2] = "10";

                gl[0] = pow(10, -std::stod(l[0])*0.1);
                gl[1] = pow(10, -std::stod(l[1])*0.1);
                gl[2] = pow(10, -std::stod(l[2])*0.1);

                if(gl.sum() > 0)
                    return gl;
            }
            catch (...) {
                gl[0] = NAN;
                gl[1] = NAN;
                gl[2] = NAN;
            }
        }
    }

	//try to get GT
    gl = getGT(split[static_cast<size_t>(indexGT)]);
	
	return gl;
}

/**
Extracts genotype call from GT for a single sample.

@param column VCF field that contains PL or GL value for one sample.
@param indexPL Index of PL values after column is split by ':'.
@param indexPL Index of GL values after column is split by ':'.
@param indexGT Index of GT values after column is split by ':'.
@return Genotype call (0, 1 or 2). NAN if issue in parsing.
*/
double getVCFGenotypeCall(std::string &column, int indexGT) {

    std::vector<std::string> split = splitString(column, ':');
    int size = static_cast<int>(split.size());

    if (indexGT > -1 && size > indexGT) {
        std::string gt = split[static_cast<size_t>(indexGT)];

        if (gt[0] == '0'){
            if(gt[2] == '0')
                return 0;
            else if (gt[2] == '1')
                return 1;
        }
        else if (gt[0] == '1'){
            if(gt[2] == '0')
                return 1;
            else if (gt[2] == '1')
                return 2;
        }
    }

    return NAN;
}
