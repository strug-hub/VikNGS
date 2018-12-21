#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "../Variant.h"
#include "../Enum/TestSettings.h"
#include <time.h>

static const std::string dfile = "/debug";
static const std::string pfile = "/pvalues";
static const std::string ffile = "/filtered";
static const std::string txt = ".txt";


inline std::string fileName(std::string outputDir, std::string base, std::string name){
    if(name.size() > 0)
        return outputDir + base + "_" + name + txt;
    else
        return outputDir + base + txt;
}

inline void initializeOutputFiles (std::string outputDir, std::string name = ""){

    std::ofstream pvals(fileName(outputDir, pfile, name));
    pvals.close();

    std::ofstream filtered(fileName(outputDir, ffile, name));
    filtered.close();
}

inline void outputPvals(std::vector<VariantSet>& variants, std::string outputDir, std::vector<TestSettings>& test, std::string name="") {

    std::ofstream pvals(fileName(outputDir, pfile, name), std::ios_base::app);

    if (pvals.is_open())
    {
        for (size_t i = 0; i < variants.size(); i++)
            pvals << variants[i].toString(test[0].toShortString(), 0, i) << '\n';

        pvals.close();
    }
}


inline void outputFiltered(std::vector<Variant> variants, std::string outputDir, std::string name = "") {

    std::ofstream filtered(fileName(outputDir, ffile, name), std::ios_base::app);

    if (filtered.is_open())	{
        for (size_t i = 0; i < variants.size(); i++)
            filtered << variants[i].toString() + "\t" + filterToString(variants[i].getFilter()) << std::endl;

        filtered.close();
    }
}


inline void outputDebug(std::string line, std::string outputDir, std::string name = "") {

    //std::string dfile = outputDir + "/debug.txt";
    std::ofstream debug(fileName(outputDir, dfile, name), std::ios_base::app);

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

