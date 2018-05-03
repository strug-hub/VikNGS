#include "VCFParserUtils.h"

static const std::string VCF_PARSER = "VCF parser";

/*
Reads every line (variant) from a VCF file.

@param vcfDir Directory of multisample VCF file.
@param minPos Variants less than this POS will not get parsed.
@param maxPos Variants greater than this POS will not get parsed.
@return A list of Variant parsed from the file at vcfDir.
*/
std::vector<Variant> parseVCFLines(std::string vcfDir, int minPos, int maxPos) {

	File vcf;
	vcf.open(vcfDir);

    std::vector<Variant> variants;

	//skips header
	extractHeader(vcf);

	int lineCount = 0;

	if(minPos > -1){
		while (vcf.hasNext()) {
			std::string line = vcf.nextLine();
			std::vector<std::string> columns = split(line, VCF_SEPARATOR);
		
		int pos = std::stoi(columns[LOC]);	
		if(pos >= minPos)
			break;
		}
	}
	
	while (vcf.hasNext()) {


        if(lineCount % 1000 == 0){
            if(lineCount == 0)
                printInfo("Parsing VCF file...");
            else
                printInfo( std::to_string(lineCount) + " variant lines have been parsed");
        }

		std::string line = vcf.nextLine();
		std::vector<std::string> columns = split(line, VCF_SEPARATOR);
			
        Variant variant = constructVariant(columns);
		lineCount++;

		if (variant.isValid()){
			if(maxPos > -1 && variant.pos > maxPos)
				break;

			variants.push_back(variant);
		}
		else {
			std::string lineNumber = std::to_string(vcf.getLineNumber());
			printWarning(VCF_PARSER, "Skipping line " + lineNumber + " of VCF file - " + variant.getErrorMessage());
		}
	}

    printInfo( std::to_string(lineCount) + " variants read from VCF");

	vcf.close();
	return variants;
}
