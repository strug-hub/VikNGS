#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <math.h> 

/*
Uses EM algorithm to estimate the genotype frequencies in the sample

@param snp A SNP with genotype likelihoods.
@return A vector with three doubles to stands for probability of 0, 1 or 2 minor alleles.
*/
std::vector<double> calcEM(SNP &snp) {
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

		for (size_t i = 0; i < snp.gl.size(); i++) {
			if (snp.L00(i) == NULL)
				continue;

			glCounter++;

			pD = 1 / (p * snp.L00(i) + q * snp.L01(i) + d * snp.L11(i));
			Ep += p * snp.L00(i) * pD;
			Eq += q * snp.L01(i) * pD;
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

	std::vector<double> freq;
	freq.push_back(p);
	freq.push_back(q);
	freq.push_back(1 - p - q);

	return freq;
}


/*
Calculates the conditional expected genotype probability E( P(G_ij | D_ij) ) 
given the genotype likelihoods P(D_ij | G_ij = g) and frequencies.

E( P(G_ij | D_ij) ) = sum from g=0 to 2 P(G_ij = g | D_ij),
where P(G_ij = g | D_ij) = P(D_ij | G_ij = g) * P(G_ij = g)/P(D_ij).

@param snp SNP with genotype likelihoods P(D_ij | G=AA, Aa or aa)} for one locus.
@return A vector containing conditional expectation probability E( P(G_ij | D_ij) ).
*/
std::vector<double> calcEG(SNP &snp) {

	double m0;
	double m1;
	double m2;
	double m;

	std::vector<double> EG;

	for (size_t i = 0; i < snp.gl.size(); i++) {

		if (snp.L00(i) != NULL) {
			m0 = snp.L00(i) * snp.p[0];
			m1 = snp.L01(i) * snp.p[1];
			m2 = snp.L11(i) * snp.p[2];
			m = 1 / (m0 + m1 + m2);

			EG.push_back(m1*m + 2 * m2*m);
		}
		else
			EG.push_back(NULL);

	}

	return EG;
}


/*
Applies sigmoid function s(z) = 1/(1+exp^z)

@param z Input to sigmoid function
@return s(z).
*/
inline double sigmoid(double z){
	return (1.0 / (1.0 + exp(z))); 
}

/*
Performs logistic regression using Newton's method. Estimates intercept (B0) and
coefficient for explanatory variable (B1).

@param y Response variable (vector of 0 and 1).
@param x Explanatory variable.
@return A vector containing the estimates of B0 and B1 as the first and second elements respectively.
*/
std::vector<double> logisticRegression(std::vector<bool> &y, std::vector<double> &x) {
	double B0 = 0;
	double B = 0;
	double lastB0 = 1;
	double lastB = 1;
	double temp;
	double s;
	double det;
	double gradB0;
	double gradB;
	double a;
	double bc;
	double d;
	int iteration = 0;

	while (true) {
		iteration++;

		gradB0 = 0;
		gradB = 0;
		a = 0;
		bc = 0;
		d = 0;

		for (size_t i = 0; i < y.size(); i++) {
			if (x[i] != NULL) {
				s = sigmoid(-(B * x[i] + B0));

				//calculate 1st derivative of likelihood function
				temp = s - y[i];
				gradB0 += temp;
				gradB += temp * x[i];

				//calculate Hessian matrix
				temp = s * (1 - s);
				a += temp;
				bc += temp * x[i];
				d += temp * x[i] * x[i];
			}
		}

		temp = a*d - bc*bc;
		if (temp == 0) {
			std::cout << "Warning: Hessian not invertible (logistic regression).\n";
			std::vector<double> out;
			out.push_back(B0);
			out.push_back(B);
		}

		det = 1 / temp;

		B0 -= det * ((d * gradB0) - (bc * gradB));
		B -= det * ((a * gradB) - (bc * gradB0));

		if (iteration > 100 || (abs(B - lastB) < 1e-7 && abs(B0 - lastB0) < 1e-7 && iteration > 1)) {

			if(iteration > 100)
				std::cout << "Warning: Logistic regression failed to converge after 100 iterations.\n";

			std::vector<double> out;
			out.push_back(B0);
			out.push_back(B);
			return out;
		}
		lastB = B;
		lastB0 = B0;
	}
}

/*
Performs logistic regression using Newton's method. Estimates intercept (B0) only.

@param y Response variable (vector of 0 and 1).
@param x Explanatory variable.
@return Estimate of B0.
*/
double logisticRegressionInterceptOnly(std::vector<bool> &y, std::vector<double> &x) {
	double B0 = 0;
	double B = 0;
	double lastB0 = 1;
	double s;
	double gradB0;
	double hessB0;
	int iteration = 0;

	while (true) {
		iteration++;
		gradB0 = 0;
		hessB0 = 0;
		s = sigmoid(-B0);

		for (size_t i = 0; i < y.size(); i++) {
			if (x[i] != NULL) {

				//calculate 1st derivative of likelihood function
				gradB0 += s - y[i];

				//calculate 2nd derivative of likelihood function
				hessB0 += s * (1 - s);
			}
		}

		B0 -= gradB0/hessB0;

		if (iteration > 100 || (abs(B0 - lastB0) < 1e-7 && iteration > 1)) {
			if (iteration > 100)
				std::cout << "Warning: Logistic regression failed to converge after 100 iterations.\n";
			return B0;
		}
		lastB0 = B0;
	}
}

/*
Calculates score test statistic.
score = sum((Y - Yhat) * X)
testStastic = score^2 / [sum((Y - Yhat)^2 * var(X))]

@param y Vector with phenotypes (case/control).
@param x Vector with E(G | D).
@return Test statistic following a chi-squared distribution with 1 degree of freedom.
*/
double scoreTest(std::vector<bool> &y, std::vector<double> &x) {
	
	double ybar = 0;
	double xbar = 0;

	double n = 0;

	for (size_t i = 0; i < y.size(); i++)
		if (x[i] != NULL) {
			ybar += y[i];
			xbar += x[i];
			n++;
		}

	ybar /= n;
	xbar /= n;
	double xvar = 0;

	for (size_t i = 0; i < x.size(); i++)
		if (x[i] != NULL)
			xvar += pow(x[i] - xbar, 2);

	xvar /= (n);
	double score = 0;
	double dnom = 0;
	double temp;

	for (size_t i = 0; i < y.size(); i++)
		if (x[i] != NULL) {
			temp = y[i] - ybar;

			score += temp * x[i];
			dnom += pow(temp, 2) * xvar;
		}
	return pow(score, 2) / dnom;
}

/*
Finds p-value for test statistic using a chi-squared distribution with one degree of freedom
using chi-squared probability density function.

@param statistic Test statistic.
@return p-value.
*/
double chiSquareOneDOF(double statistic){
	double z = statistic * 0.5;
	double sc = 2 * sqrt(z) * exp(-z);

	double sum = 1;
	double prevSum = sum;
	double nom = 1;
	double dnom = 1;
	double s = 0.5;

	for (int i = 0; i < 200; i++)
	{
		nom *= z;
		s++;
		dnom *= s;
		sum += nom / dnom;
		if (prevSum == sum) break;
		prevSum = sum;
	}

	double p = sum * sc;
	if (isnan(p) || isinf(p) || p <= 1e-8)
		return 1e-14;

	p /= tgamma(0.5);

	return 1 - p;
}


//Logistic regression with gradient descent (slower).
/*
std::vector<double> LogisticRegression(std::vector<bool> &y, std::vector<double> &x) {

	double w = 0;
	double w0 = 0;
	double learningRate = 0.02;

	double wlast = 100;
	double w0last = 100;

	double wgrad;
	double w0grad;

	double tmp = 0;
	int iteration = 0;

	while (true) {
		iteration++;
		wgrad = 0;
		w0grad = 0;
		for (size_t i = 0; i < y.size(); i++) {
			if (x[i] != NULL) {
				tmp = y[i] - (1.0 / (1.0 + exp(-(w * x[i] + w0))));
				wgrad += x[i] * tmp;
				w0grad += tmp;
			}
		}

		w += learningRate * wgrad;
		w0 += learningRate * w0grad;

		if (iteration > 1000 || (abs(w - wlast) < 0.0000001 && abs(w0 - w0last) < 0.0000001)) {
			std::vector<double> out;
			out.push_back(w0);
			out.push_back(w);
			return out;
		}

		wlast = w;
		w0last = w0;
	}
}


double LogisticRegressionInterceptOnly(std::vector<bool> &y, std::vector<double> &x) {

	double w0 = 0;
	double learningRate = 0.02;

	double w0last = 100;
	double w0grad;

	int iteration = 0;

	while (true) {
		iteration++;
		w0grad = 0;
		for (size_t i = 0; i < y.size(); i++)
			if (x[i] != NULL)
				w0grad += y[i] - (1.0 / (1.0 + exp(-w0)));

		w0 += learningRate * w0grad;

		if (iteration > 1000 || abs(w0 - w0last) < 0.0000001)
			return w0;

		w0last = w0;
	}
}
*/