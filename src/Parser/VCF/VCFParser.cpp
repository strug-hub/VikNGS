#include "VCFParserUtils.h"

static const std::string VCF_PARSER = "VCF parser";

std::vector<Variant> parseVCFLines(Request &req, VectorXd &Y, std::string family) {

    std::string vcfDir = req.vcfDir;
    int minPos = req.minPos;
    int maxPos = req.maxPos;

	File vcf;
	vcf.open(vcfDir);
    std::vector<Variant> variants;
    std::vector<std::string> filterInfo;
    std::vector<int> filterCode;

	//skips header
	extractHeader(vcf);

	int lineCount = 0;

    //skip until POS in VCF >= minPos
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

        if(lineCount % 5000 == 0){
            if(lineCount == 0)
                printInfo("Parsing VCF file...");
            else
                printInfo( std::to_string(lineCount) + " variant lines have been parsed");
        }

        std::string line = vcf.nextLine();
        std::vector<std::string> columns = split(line, VCF_SEPARATOR);

        lineCount +=1;

        Variant variant = constructVariant(columns);
		lineCount++;

		if (variant.isValid()){
			if(maxPos > -1 && variant.pos > maxPos)
				break;

            variant = calculateExpectedGenotypes(variant);
            int code = filterVariant(req, variant, Y, family);
            if(code > 0){
                filterInfo.push_back(variant.toString());
                filterCode.push_back(code);
            }
            else
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
