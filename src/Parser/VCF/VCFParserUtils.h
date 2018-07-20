#pragma once
#include "../InputParser.h"

static const char VCF_SEPARATOR = '\t';

static const int CHR = 0;
static const int LOC = 1;
static const int REF = 3;
static const int ALT = 4;
static const int FILTER = 6;
static const int FORMAT = 8;

/*
Extracts the last header row from a VCF file

@param vcf A VCF file.
@throws runtime_error if cannot find a header in file.
@return The names of each column in the last header row from vcf.
*/
std::vector<std::string> extractHeader(File &vcf);
std::string extractHeaderLine(File &vcf);


/*
Calculates genotype likelihood from PL or GL for a single sample.

@param column VCF field that contains PL or GL value for one sample.
@param indexPL Index of PL values after column is split by ':'.
@param indexPL Index of GL values after column is split by ':'.
@param indexGT Index of GT values after column is split by ':'.
@throws runtime_error if column cannot be parsed
@return A GenotypeLikelihood object.
*/
GenotypeLikelihood getGenotypeLikelihood(std::string &column, int indexPL, int indexGL, int indexGT, Variant& variant);
double getGenotypeCall(std::string &column, int indexGT);


/*
Reads every sample ID (columns after FORMAT) from a multisample VCF and stores it in a map.

@param vcfDir Directory of multisample VCF file.
@return A map from sample name to a unique integer.
*/
std::map<std::string, int> getSampleIDMap(std::string vcfDir);

/*
Takes the columns of a VCF line and builds a Variant object.

@param columns Column values corresponding to a VCF line.
@param onlyGT Extract genotype calls only, no likelihood.

@return Variant object.
*/
Variant constructVariant(std::vector<std::string> &columns);
