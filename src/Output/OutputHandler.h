#pragma once

#include "../vikNGS.h"
#include <iostream>
#include <fstream>

static const std::string dfile = "/debug.txt";
static const std::string pfile = "/pvalues_.txt";
static const std::string ffile = "/filtered.txt";

void initializeOutputFiles (std::string outputDir){

    std::ofstream debug(outputDir + dfile);
    debug.close();

    std::ofstream pvals(outputDir + pfile);
	pvals.close();

    std::ofstream filtered(outputDir + ffile);
	filtered.close();
}

void outputPvals(std::vector<VariantSet> &variants, std::string outputDir) {
	
    //std::string pfile = outputDir + "/pvalues_.txt";
    std::ofstream pvals(outputDir + pfile, std::ios_base::app);

	if (pvals.is_open())
	{
        for (size_t i = 0; i < variants.size(); i++) {

            //todo

            //pvals << variants[i].toString() << '\t';
            //pvals << variants[i].getPval(0) << '\t';
            //pvals << variants[i].getPvalSourceShort(0) << '\n';

		}
		pvals.close();
	}
}


void outputFiltered(std::vector<std::string> variantInfo, std::vector<int> failCode,
                           std::vector<std::string> codeMap, std::string outputDir) {

    //todo
    //std::string ffile = outputDir + "/filtered.txt";
    std::ofstream filtered(outputDir + ffile, std::ios_base::app);

	if (filtered.is_open())
	{
        for (size_t i = 0; i < variantInfo.size(); i++) {
            filtered << variantInfo[i] << '\t';
            filtered << codeMap[failCode[i]] << std::endl;
		}
		filtered.close();
	}
}

void outputFiltered(std::vector<Variant> variants, std::string explain, std::string outputDir) {

    //todo
    //std::string ffile = outputDir + "/filtered.txt";
    std::ofstream filtered(outputDir + ffile, std::ios_base::app);

    if (filtered.is_open())
    {
        for (size_t i = 0; i < variants.size(); i++) {

            //filtered << variants[i].toString() << '\t';
            filtered << explain << '\n';
        }
        filtered.close();
    }
}

void outputDebug(std::string line, std::string outputDir) {

    //std::string dfile = outputDir + "/debug.txt";
    std::ofstream debug(outputDir + dfile, std::ios_base::app);

    if (debug.is_open())
    {
        debug << line << std::endl;
        debug.close();
    }
}

