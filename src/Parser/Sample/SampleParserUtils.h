#pragma once
#include "../InputParser.h"

/*
Parses a data file containing tab-separated info for each sample.

@param covariate Column of covariates read from sample data file.
@param isNumeric True if covariate is a numeric value and not categorical.
@param columnIndex Index of the covariate column.
@return Matrix of covariate values
*/
MatrixXd handleCovariate(std::vector<std::string> covariate, bool isNumeric, int columnIndex);

/*
Parses a data file containing tab-separated info for each sample.

@param covariates Column(s) of covariates read from sample data file.
@return Matrix of covariate values
*/
MatrixXd handleCovariates(std::vector<std::vector<std::string>> covariates);