#include "../RVS.h"

#include <iostream>  
#include <fstream>

void outputPvals(std::vector<std::string> info, std::vector<double> pvalues, std::string outputDir) {
	std::ofstream out(outputDir);

	int precise = 12;

	if (out.is_open())
	{
		for (size_t i = 0; i < pvalues.size(); i++) {

			out << info[i];
			out << '\t';
			out << pvalues[i];
			out << '\n';
		}
		out.close();
	}


}
