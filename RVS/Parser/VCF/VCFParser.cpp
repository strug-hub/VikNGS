#include "VCFParserUtils.h"

static const std::string VCF_PARSER = "VCF parser";

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

	//skips header
	extractHeader(vcf);

	while (vcf.hasNext()) {


		std::string line = vcf.nextLine();
		std::vector<std::string> columns = split(line, VCF_SEPARATOR);

		VCFLine variant = constructVariant(columns);

		if (variant.isValid())
			variants.push_back(variant);
		else {
			std::string lineNumber = std::to_string(vcf.getLineNumber());
			printWarning(VCF_PARSER, "Skipping line " + lineNumber + " of VCF file - " + variant.getErrorMessage());
		}
	}

	vcf.close();
	return variants;
}