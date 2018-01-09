#include "VCFParserUtils.h"

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
	}
	catch (...) { 
		throw;
	}

	while (vcf.hasNext()) {

		std::string line = vcf.nextLine();
		std::vector<std::string> columns = split(line, VCF_SEPARATOR);

		VCFLine variant = constructVariant(columns);

		if (variant.isValid())
			variants.push_back(variant);
		else {
			std::string warningMessage = "Error while parsing variant on line " + std::to_string(vcf.getLineNumber()) + ": ";
			warningMessage += variant.getErrorMessage();
			printWarning(warningMessage);
		}
	}

	vcf.close();
	return variants;
}