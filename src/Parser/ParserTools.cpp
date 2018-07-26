#include "InputParser.h"
#include <regex>
#include <fstream>

/**
Create a vector of indexes used for collapsing by k

@param k number of variants to collapse.
@param n total number of variants.
@return Collapsed indexes.
*/
std::vector<std::vector<int>> collapseEveryK(int k, int n){
  std::vector<std::vector<int>> collapsedVariants;

  int index = 0;

  while(index < n){
    std::vector<int> group;

    for (int i = 0; i < k; i++) {
	group.push_back(index);
	index++;
	if(index >= n)
	  break;
    }

    collapsedVariants.push_back(group);
  }

  return collapsedVariants;
}

/**
Extracts string from a MemoryMapped class

@param charArray MemoryMapped to extract string from.
@param start Index to start extracting string from.
@param end Index to stop extracting string from.
@return Extracted string.
*/
std::string extractString(MemoryMapped &charArray, int start, int end) {
	std::string ret;
	for (; start < end; start++) { ret += charArray[start]; }
	return ret;
}

/**
Determines if response variable is case/control(1/0) or quantitative.
Returns family type (case/control = binomial, quantitative = normal)

@param VectorXd Y vector of response variable.
@return Distribution family for response variable.
*/
std::string determineFamily(VectorXd Y) {

    //if a value not 0 or 1 is found, assume quantitative data
    for(int i = 0; i < Y.rows(); i++){
        if(Y[i] != 0 && Y[i] != 1){
            printInfo("Quantitative data detected");
            return "normal";
        }
    }

    printInfo("Case/control data detected");
    return "binomial";
}


/**
Removes whitespace from a string

@param str String to remove whitespace from.
@return str without whitespace.
*/
std::string trim(std::string str){
	str.erase(std::remove(str.begin(), str.end(), '\t'), str.end());
	str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
	return str;
}

/**
Separates a string into a vector, splitting at every postion with
the sep character

@param s String to split.
@param sep Character to split the string at.
@return Split string.
*/
std::vector<std::string> split(std::string &s, char sep) {
	std::vector<std::string> split;
	int start = 0;
	for (int i = 0; i <= s.length(); i++) {
		if (s[i] == sep || i == s.length()) {
            split.emplace_back(s.substr(start, i - start));
			start = i + 1;
		}
	}
	return split;
}

/**
Checks if a string is a number

@param str string to check
@return true if str is numeric
*/
inline bool isNumeric(const std::string& str) {
	return (std::regex_match(str, std::regex("-?[1234567890]+(\\.[1234567890]+)?")));
}

/*
Uses EM algorithm to estimate the genotype frequencies in the sample

@param likelihood A vector with genotype likelihoods.
@return A vector with probability of 0, 1 or 2 minor alleles.
*/
VectorXd calcEM(std::vector<GenotypeLikelihood> &likelihood) {
	double p = 0.15;
	double q = 0.15;
	double qn = 1;
	double pn = 0;
	double dn = 0;
	double d = 0;

	int glCounter = 0;

	double Ep;
	double Eq;
	double pD;

	int k = 0;

	while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

		d = 1 - p - q;
		glCounter = 0;
		Ep = 0;
		Eq = 0;

		for (int i = 0; i < likelihood.size(); i++) {
			if (likelihood[i].missing)
				continue;

			glCounter++;

			pD = 1 / (p * likelihood[i].L00 + q * likelihood[i].L01 + d * likelihood[i].L11);
			Ep += p * likelihood[i].L00 * pD;
			Eq += q * likelihood[i].L01 * pD;
		}

		pn = p;
		qn = q;
		dn = 1 - q - p;
		p = Ep / glCounter;
		q = Eq / glCounter;

		k++;
		if (k == 1000)
			break;
	}

	VectorXd freq(3);
    freq[0] = std::max(1e-14, p);
    freq[1] = std::max(1e-14, q);
    freq[2] = std::max(1e-14, 1 - p - q);

	return freq;
}

/*
Calculates the conditional expected genotype probability E( P(G_ij | D_ij) )
given the genotype likelihoods P(D_ij | G_ij = g) and frequencies.

E( P(G_ij | D_ij) ) = sum from g=0 to 2 P(G_ij = g | D_ij),
where P(G_ij = g | D_ij) = P(D_ij | G_ij = g) * P(G_ij = g)/P(D_ij).

@param likelihood A vector with genotype likelihoods.
@param p Vector with probability of 0, 1 or 2 minor alleles.
@return Vector containing conditional expectation probability E( P(G_ij | D_ij) ).
*/
VectorXd calcEG(std::vector<GenotypeLikelihood> &likelihood, VectorXd &p) {

	double m0;
	double m1;
	double m2;
	double m;

	VectorXd EG(likelihood.size());

	for (int i = 0; i < likelihood.size(); i++) {

		if (!likelihood[i].missing) {

            m0 = likelihood[i].L00 * p[0];
            m1 = likelihood[i].L01 * p[1];
            m2 = likelihood[i].L11 * p[2];
			m = 1 / (m0 + m1 + m2);

            EG[i] = m1*m + (2 * m2*m);
		}
		else
			EG[i] = NAN;
	}

	return EG;
}

/*
Generates the expected probabilities of the genotypes E(G_ij | D_ij).
Using the genotype likelihood for case and controls from VCF file to generates the population frequency by calling function calcEM
and then use it to calculate the expected genotype probabilities E(G_ij | D_ij) by calling function calcEG.
Variants with homozygous call in the whole sample (standard deviation of their E(G_ij | D_ij) < 10^4) are removed.

@param variants Variant from VCF file
@return none
@effect Expected genotypes are added to variant object
*/
void calculateExpectedGenotypes(Variant &variant) {

    variant.P = calcEM(variant.likelihood);
    variant.expectedGenotype = calcEG(variant.likelihood, variant.P);
}
