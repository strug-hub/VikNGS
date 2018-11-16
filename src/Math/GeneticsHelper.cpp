#include "Math.h"
#include <iostream>

/**
Calculates the robust variance of E(G | D). var(x) = E(x^2) - E(x)^2

@param p Genotype frequency from EM algorithm
@return Robust variance of E(G | D)
*/
double calcRobustVar(Vector3d &P) {
    return (4 * P[2] + P[1]) - pow(2 * P[2] + P[1], 2);
}

double calcRobustVar(double P1, double P2) {
    return (4 * P2 + P1) - pow(2 * P2 + P1, 2);
}

/**
Generates the genotype calls using simple bayesian genotyper (maximum likelihood).

@param  gl A vector of genotype likelihoods, one per individual.
@param  p Population genotype frequency.

@return Vector of genotype calls.
*/
VectorXd calculateGenotypeCalls(std::vector<Vector3d>& gl, Vector3d &P){

    VectorXd GC(gl.size());

    int call;
    int index;
    for (size_t i = 0; i < gl.size(); i++) {
        index = static_cast<int>(i);
        if(std::isnan(gl[i][0])){
            GC[index] = NAN;
            continue;
        }
        call = maxValue(gl[i][0] * P[0], gl[i][1] * P[1], gl[i][2] * P[2]);

        if(call < 0)
            GC[index]=NAN;
        else
            GC[index] = call;
    }

    return GC;
}

/**
Estimates the genotype frequencies from provided genotypes


@param likelihood A vector with genotypes.
@return A vector with probability of 0, 1 or 2 minor alleles.
*/
Vector3d calculateGenotypeFrequencies(VectorXd & gt) {

    Vector3d P;
    P[0] = 0; P[1] = 0; P[2] = 0;
    double n = 0;
    //expecting genotypes to be 0, 1, or 2
    //using 0.1, 1.1 and 2.1 for double comparison
    for(int i = 0; i < gt.rows(); i++){
        if(std::isnan(gt[i]))
            continue;
        if(gt[i] < 0.1)
            P[0]++;
        else if (gt[i] < 1.1)
            P[1]++;
        else if (gt[i] < 2.1)
            P[2]++;
        n++;
    }

    P[0]=P[0]/n;
    P[1]=P[1]/n;
    P[2]=P[2]/n;

    return P;
}


/**
Uses EM algorithm to estimate the genotype frequencies in the sample

@param likelihood A vector with genotype likelihoods.
@return A vector with probability of 0, 1 or 2 minor alleles.
*/
Vector3d calculateGenotypeFrequencies(std::vector<Vector3d>& likelihood) {
    double p = 0.15;
    double q = 0.15;
    double qn = 1;
    double pn = 0;
    double dn = 0;
    double d = 0;

    double n = 0;

    double Ep;
    double Eq;
    double pD;

    int k = 0;

    while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

        d = 1 - p - q;
        n = 0;
        Ep = 0;
        Eq = 0;

        for (size_t i = 0; i < likelihood.size(); i++) {
            if (std::isnan(likelihood[i][0]))
                continue;

            n++;

            pD = 1.0 / ((p * likelihood[i][0]) + (q * likelihood[i][1]) + (d * likelihood[i][2]));
            Ep += p * likelihood[i][0] * pD;
            Eq += q * likelihood[i][1] * pD;
        }

        pn = p;
        qn = q;
        dn = 1 - q - p;
        p = Ep / n;
        q = Eq / n;

        k++;
        if (k == 1000)
            break;
    }

    Vector3d freq;
    freq[0] = std::max(1e-14, p);
    freq[1] = std::max(1e-14, q);
    freq[2] = std::max(1e-14, 1 - p - q);

    return freq;
}

/**
Calculates the conditional expected genotype probability E( P(G_ij | D_ij) )
given the genotype likelihoods P(D_ij | G_ij = g) and frequencies.

E( P(G_ij | D_ij) ) = sum from g=0 to 2 P(G_ij = g | D_ij),
where P(G_ij = g | D_ij) = P(D_ij | G_ij = g) * P(G_ij = g)/P(D_ij).

@param likelihood A vector with genotype likelihoods.
@param p Vector with probability of 0, 1 or 2 minor alleles.
@return Vector containing conditional expectation probability E( P(G_ij | D_ij) ).
*/
VectorXd calculateExpectedGenotypes(std::vector<Vector3d> &likelihood, Vector3d &P) {

    double m0;
    double m1;
    double m2;
    double m;

    VectorXd EG(likelihood.size());

    for (size_t i = 0; i < likelihood.size(); i++) {

        if (!std::isnan(likelihood[i][0])){

            m0 = likelihood[i][0] * P[0];
            m1 = likelihood[i][1] * P[1];
            m2 = likelihood[i][2] * P[2];
            m = 1 / (m0 + m1 + m2);

            EG[i] = m1*m + (2 * m2*m);
        }
        else
            EG[i] = NAN;
    }


    return EG;
}
