#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "../Variant.h"

static const std::string dfile = "/debug.txt";
static const std::string pfile = "/pvalues.txt";
static const std::string ffile = "/filtered.txt";


inline void initializeOutputFiles (std::string outputDir){

    //std::ofstream debug(fileName(outputDir, dfile));
    //debug.close();

    std::ofstream pvals(outputDir + pfile);
	pvals.close();

    std::ofstream filtered(outputDir + ffile);
	filtered.close();
}

inline void outputPvals(std::vector<VariantSet>& variants, std::string outputDir, std::vector<Test>& test) {
	
    std::ofstream pvals((outputDir + pfile), std::ios_base::app);

	if (pvals.is_open())
	{
        for (size_t i = 0; i < variants.size(); i++)
            pvals << variants[i].toString(test[0].toShortString(), 0, i) << '\n';

		pvals.close();
	}
}


inline void outputFiltered(std::vector<Variant> variants, std::string outputDir) {

    //std::string ffile = outputDir + "/filtered.txt";
    std::ofstream filtered((outputDir+ ffile), std::ios_base::app);

	if (filtered.is_open())	{
        for (size_t i = 0; i < variants.size(); i++)
            filtered << variants[i].toString() + "\t" + filterToString(variants[i].getFilter()) << std::endl;

		filtered.close();
	}
}


inline void outputDebug(std::string line, std::string outputDir) {

    //std::string dfile = outputDir + "/debug.txt";
    std::ofstream debug((outputDir+ dfile), std::ios_base::app);

    if (debug.is_open())
    {
        debug << line << std::endl;
        debug.close();
    }
}


inline void outputMatrix(MatrixXd M, std::string filename) {
    std::ofstream m(filename + ".txt", std::ios_base::app);

    for(int i=0; i < M.rows(); i++){
        std::string line = "";

        for(int j=0; j < M.cols(); j++){

            m << std::setprecision(16) << M(i,j);

            if(j < M.cols() -1)
                m << "\t";
        }
        m << std::endl;
    }
    if (m.is_open())
        m.close();

}

inline void outputVector(VectorXd V, std::string filename) {
    std::ofstream v(filename + ".txt", std::ios_base::app);

    for(int i=0; i < V.rows(); i++){
        v << std::setprecision(24) << V[i] << std::endl;
    }

    if (v.is_open())
        v.close();

}



