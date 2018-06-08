#include "SampleParserUtils.h"

static const std::string SAMPLE_PARSER_UTILS = "sample parser utils";


MatrixXd handleCovariate(std::vector<std::string> covariate, bool isNumeric, int columnIndex) {

	if (isNumeric) {
		VectorXd z(covariate.size());
		for (int i = 0; i < covariate.size(); i++) {

			if (covariate[i] == "NA" || covariate[i] == "")
				z[i] = NAN;
			else
				z[i] = stod(covariate[i]);
		}
		return z;
	}
	else {
		std::map<std::string, int> covIDMap;
		int id = 0;

		for (int i = 0; i < covariate.size(); i++) {
			if (!covIDMap.count(covariate[i])) {
				covIDMap[covariate[i]] = id;
				id++;
			}
		}

		if (covIDMap.size() > 5) {
			std::string message = "In sample data file: Covariate in column ";
			message += std::to_string(columnIndex + 5) + " is being treated as categorical and has ";
			message += std::to_string(covIDMap.size()) + " unique values. This could have a negative impact on computation time.";
			printWarning(message);
		}

        MatrixXd z(covariate.size(), covIDMap.size() - 1);
		z = z.setZero();

		int col;
		for (int i = 0; i < covariate.size(); i++) {

			if (covariate[i] == "NA" || covariate[i] == "") {
                for (int j = 0; j < covIDMap.size() - 1; j++)
					z(i, j) = NAN;
			}
			else {
				col = covIDMap[covariate[i]];

                //don't include last categorical covariate
                if(col < id-1)
                    z(i, col) = 1;
			}
		}
		return z;
	}
}


MatrixXd handleCovariates(std::vector<std::vector<std::string>> covariates) {
	MatrixXd Z;

	for (int i = 0; i < covariates[0].size(); i++) {

		std::vector<std::string> covariate;
		bool isNumeric = true;

		for (int j = 0; j < covariates.size(); j++) {
			covariate.push_back(covariates[j][i]);
			if (isNumeric) {
				try {
					std::stod(covariates[j][i]);
				}
				catch (...) {
					std::string cov = covariates[j][i];
					if (cov != "NA" && cov != "")
						isNumeric = false;
				}
			}
		}
		MatrixXd cov = handleCovariate(covariate, isNumeric, i);
		
		if (i == 0)
			Z = cov;
		else {
			MatrixXd temp = Z;
			MatrixXd T(temp.rows(), temp.cols() + cov.cols());
			Z=T;
			Z << temp, cov;
		}
	}

    MatrixXd interceptZ(Z.rows(), Z.cols()+1);
    interceptZ << Eigen::VectorXd::Constant(Z.rows(), 1), Z;
    return interceptZ;
}
