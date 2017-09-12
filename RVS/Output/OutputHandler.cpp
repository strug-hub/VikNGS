#include "../RVS.h"
#include <iostream>  
#include <fstream>

void outputPvals(std::vector<std::vector<double>> pvalues, std::string outputDir) {
	std::ofstream out(outputDir);

	int precise = 12;

	if (out.is_open())
	{
		if(pvalues[0].size() > 1)
			out << "Linear p-value\tQuadratic p-value\n";

		for (size_t i = 0; i < pvalues.size(); i++) {

			out << pvalues[i][0];
			if (pvalues[i].size() > 1) {
				out << '\t';
				out << pvalues[i][1];
			}

			out << '\n';
		}
		out.close();
	}


}
