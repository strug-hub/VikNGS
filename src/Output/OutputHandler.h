#pragma once

#include "../Variant.h"
#include <iostream>  
#include <fstream>


inline void createFile(std::string outputDir){

    std::string dfile = outputDir + "/debug.txt";
    std::string pfile = outputDir + "/pvalues.txt";
    std::string ffile = outputDir + "/filtered.txt";

    std::ofstream debug(dfile);
    debug.close();

	std::ofstream pvals(pfile);	
	pvals.close();

	std::ofstream filtered(ffile);	
	filtered.close();
}

inline void outputPvals(std::vector<Variant> &variants, std::string outputDir) {
	
    std::string pfile = outputDir + "/pvalues.txt";
	std::ofstream pvals(pfile, std::ios_base::app);

	if (pvals.is_open())
	{
        for (size_t i = 0; i < variants.size(); i++) {

            pvals << variants[i].toString() << '\t';
            pvals << variants[i].getPval(0) << '\n';
		}
		pvals.close();
	}
}


inline void outputFiltered(std::vector<std::string> variantInfo, std::vector<int> failCode,
                           std::vector<std::string> codeMap, std::string outputDir) {

    std::string ffile = outputDir + "/filtered.txt";
	std::ofstream filtered(ffile, std::ios_base::app);

	if (filtered.is_open())
	{
        for (size_t i = 0; i < variantInfo.size(); i++) {
            filtered << variantInfo[i] << '\t';
            filtered << codeMap[failCode[i]] << std::endl;
		}
		filtered.close();
	}
}

inline void outputFiltered(std::vector<Variant> variants, std::string explain, std::string outputDir) {

    std::string ffile = outputDir + "/filtered.txt";
    std::ofstream filtered(ffile, std::ios_base::app);

    if (filtered.is_open())
    {
        for (size_t i = 0; i < variants.size(); i++) {

            filtered << variants[i].toString() << '\t';
            filtered << explain << '\n';
        }
        filtered.close();
    }
}

inline void outputDebug(std::string line, std::string outputDir) {

    std::string dfile = outputDir + "/debug.txt";
    std::ofstream debug(dfile, std::ios_base::app);

    if (debug.is_open())
    {
        debug << line << std::endl;
        debug.close();
    }
}

