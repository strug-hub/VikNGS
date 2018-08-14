//CompQuadForm: Distribution Function of Quadratic Forms in Normal Variables
//
//	Computes the distribution function of quadratic forms in normal variables using Imhof's method, Davies's algorithm, Farebrother's algorithm or Liu et al.'s algorithm.
//	Version : 1.4.2
//	Published : 2016 - 09 - 20
//	Author: 	P.Lafaye de Micheaux
//  License: 	GPL-2 | GPL-3 [expanded from: GPL (≥ 2)]
//  https://CRAN.R-project.org/package=CompQuadForm 
//  Motified by Scott Mastromatteo
//  assumes 1 degree of freedom, non-centrality parameter = 0

#include "Math.h"

#include <setjmp.h>
#include <iostream>  
#define UseDouble 0             /* all floating point double */

#define TRUE  1
#define FALSE 0
typedef int BOOL;

#define pi 3.14159265358979
#define log28 0.0866  /*  log(2.0) / 8.0  */

static double sigsq, lmax, lmin, mean, c;
static double intl, ersm;
static int count, r, lim;  static BOOL ndtsrt, fail;
static int *th; static double *lb, *nc;
static jmp_buf env;

static inline double exp1(double x)               /* to avoid underflows  */
{
	return x < -50.0 ? 0.0 : std::exp(x);
}

static void counter(void)
	/*  count number of calls to errbd, truncation, cfe */
{
	extern int count, lim;
	count = count + 1;
	if (count > lim) longjmp(env, 1);
}

static inline double square(double x) { return x * x; }

static inline double cube(double x) { return x * x * x; }

static double log1(double x, BOOL first)
	/* if (first) std::log(1 + x) ; else  std::log(1 + x) - x */
{
	if (fabs(x) > 0.1)
	{
		return (first ? std::log(1.0 + x) : (std::log(1.0 + x) - x));
	}
	else
	{
		double s, s1, term, y, k;
		y = x / (2.0 + x);  term = 2.0 * cube(y);  k = 3.0;
		s = (first ? 2.0 : -x) * y;
		y = square(y);
		for (s1 = s + term / k; s1 != s; s1 = s + term / k)
		{
			k = k + 2.0; term = term * y; s = s1;
		}
		return s;
	}
}

static void order(void)
	/* find order of absolute values of lb */
{
	int j, k; double lj;
	extern double *lb; extern int *th; extern int r; extern BOOL ndtsrt;
	for (j = 0; j < r; j++)
	{
		lj = std::fabs(lb[j]);
		for (k = j - 1; k >= 0; k--)
		{
			if (lj > std::fabs(lb[th[k]]))  th[k + 1] = th[k];
			else goto l1;
		}
		k = -1;
	l1:
		th[k + 1] = j;
	}
	ndtsrt = FALSE;
}


static double errbd(double u, double* cx)
	/*  find bound on tail probability using mgf, cutoff
	point returned to *cx */
{
	double sum1, lj, x, y, xconst; int j;
	extern double sigsq, *lb; extern int r;
	counter();
	xconst = u * sigsq;  sum1 = u * xconst;  u = 2.0 * u;
	for (j = r - 1; j >= 0; j--)
	{
		lj = lb[j];
		x = u * lj; y = 1/(1.0 - x);
		xconst = xconst + lj * y;
		sum1 += square(x) * y + log1(-x, FALSE);
	}
	*cx = xconst; return exp1(-0.5 * sum1);
}

static double ctff(double accx, double* upn)
	/*  find ctff so that p(qf > ctff) < accx  if (upn > 0,
	p(qf < ctff) < accx otherwise */
{
	double u1, u2, u, rb, xconst, c1, c2;
	extern double lmin, lmax, mean;
	u2 = *upn;   u1 = 0.0;  c1 = mean;
	rb = 2.0 * ((u2 > 0.0) ? lmax : lmin);
	for (u = u2 / (1.0 + u2 * rb); errbd(u, &c2) > accx;
		u = u2 / (1.0 + u2 * rb))
	{
		u1 = u2;  c1 = c2;  u2 = 2.0 * u2;
	}
	for (u = (c1 - mean) / (c2 - mean); u < 0.9;
		u = (c1 - mean) / (c2 - mean))
	{
		u = (u1 + u2) / 2.0;
		if (errbd(u / (1.0 + u * rb), &xconst) > accx)
		{
			u1 = u; c1 = xconst;
		}
		else
		{
			u2 = u;  c2 = xconst;
		}
	}
	*upn = u2; return c2;
}

static double truncation(double u, double tausq)
	/* bound integration error due to truncation at u */
{
	double sum1, sum2, prod1, prod2, prod3, lj,
		x, y, err1, err2;
	int j, s;
	extern double sigsq, *lb; extern int r;

	counter();
	sum1 = 0.0; prod2 = 0.0;  prod3 = 0.0;  s = 0;
	sum2 = (sigsq + tausq) * square(u); prod1 = 2.0 * sum2;
	u = 2.0 * u;
	for (j = 0; j < r; j++)
	{
		lj = lb[j];
		x = square(u * lj);
		if (x > 1.0)
		{
			prod2 = prod2 + log(x);
			prod3 = prod3 + log1(x, TRUE);
			s++;
		}
		else  prod1 = prod1 + log1(x, TRUE);
	}
	sum1 = 0.5 * sum1;
	prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
	x = exp1(-sum1 - 0.25 * prod2) / pi;
	y = exp1(-sum1 - 0.25 * prod3) / pi;
	err1 = (s == 0) ? 1.0 : x * 2.0 / s;
	err2 = (prod3 > 1.0) ? 2.5 * y : 1.0;
	if (err2 < err1) err1 = err2;
	x = 0.5 * sum2;
	err2 = (x <= y) ? 1.0 : y / x;
	return  (err1 < err2) ? err1 : err2;
}

static void findu(double* utx, double accx)
	/*  find u such that truncation(u) < accx and truncation(u / 1.2) > accx */
{
	double u, ut; int i;
	static double divis[] = { 2.0, 1.4, 1.2, 1.1 };
	ut = *utx; u = ut / 4.0;
	if (truncation(u, 0.0) > accx)
	{
		for (u = ut; truncation(u, 0.0) > accx; u = ut) ut = ut * 4.0;
	}
	else
	{
		ut = u;
		for (u = u / 4.0; truncation(u, 0.0) <= accx; u = u / 4.0)
			ut = u;
	}
	for (i = 0; i < 4; i++)
	{
		u = ut / divis[i]; if (truncation(u, 0.0) <= accx)  ut = u;
	}
	*utx = ut;
}


static void integrate(int nterm, double interv, double tausq, BOOL mainx)
	/*  carry out integration with nterm terms, at stepsize
	interv.  if (! mainx) multiply integrand by
	1.0 - exp(-0.5 * tausq * u ^ 2) */
{
	double inpi, u, sum1, sum2, sum3, x, y, z;
	int k, j;
	extern double intl, ersm; extern double sigsq, c;
	extern double *lb; extern int r;
	inpi = interv / pi;
	for (k = nterm; k >= 0; k--)
	{
		u = (k + 0.5) * interv;
		sum1 = -2.0 * u * c;  sum2 = std::fabs(sum1);
		sum3 = -0.5 * sigsq * square(u);
		for (j = r - 1; j >= 0; j--)
		{
			x = 2.0 * lb[j] * u;  y = square(x);
			sum3 = sum3 - 0.25 * log1(y, TRUE);
			z = std::atan(x);
			sum1 = sum1 + z;   sum2 = sum2 + fabs(z);
		}
		x = inpi * exp1(sum3) / u;
		if (!mainx)
			x = x * (1.0 - exp1(-0.5 * tausq * square(u)));
		sum1 = std::sin(0.5 * sum1) * x;  sum2 = 0.5 * sum2 * x;
		intl = intl + sum1; ersm = ersm + sum2;
	}
}

static double cfe(double x)
	/*  coef of tausq in error when convergence factor of
	exp1(-0.5 * tausq * u ^ 2) is used when df is evaluated at x */
{
	double axl, axl1, axl2, sxl, sum1, lj; int j, k, t;
	extern BOOL ndtsrt, fail; extern int *th; extern double *lb;
	extern int r;
	counter();
	if (ndtsrt) order();
	axl = fabs(x);  sxl = (x > 0.0) ? 1.0 : -1.0;  sum1 = 0.0;
	for (j = r - 1; j >= 0; j--)
	{
		t = th[j];
		if (lb[t] * sxl > 0.0)
		{
			lj = std::fabs(lb[t]);
			axl1 = axl - lj;  axl2 = lj / log28;
			if (axl1 > axl2)  axl = axl1; else
			{
				if (axl > axl2)  axl = axl2;
				sum1 = (axl - axl1) / lj;
				for (k = j - 1; k >= 0; k--)
					sum1++;
				goto  l;
			}
		}
	}
l:
	if (sum1 > 100.0)
	{
		fail = TRUE; return 1.0;
	}
	else
		return pow(2.0, (sum1 / 4.0)) / (pi * square(axl));
}

/*
Computes P[Q > q] where Q = sum lambda_j X_j where X_j are independent random variables 
having a non-central chi^2 distribution with 1 degree of freedom and non-centrality 
parameter of 0; X_0 having a standard Gaussian distribution.

@param coef A vector with coefficients of j-th chi-squared variable.
@param evalpoint point at which df is to be evaluated.
@return P[Q > q].
*/
double qfc(std::vector<double> lambda, double evalpoint, int ncoef) {
	double acc = 1e-04;

    extern double sigsq, lmax, lmin, mean;
    extern double intl, ersm;
    extern int r, lim; extern double c;
    extern int *th; extern double *lb;

	int j, nt, ntm; double almx, xlim, xnt, xntm;
	double utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lj;

	double qfval = -1.0;

	r = ncoef; lim = 10000; c = evalpoint;
	lb = &lambda[0];
	count = 0;
	intl = 0.0; ersm = 0.0;
	qfval = -1.0; ndtsrt = TRUE;  fail = FALSE;
    xlim =  static_cast<double>(lim);
    th =  static_cast<int *>(malloc( static_cast<size_t>(r)*(sizeof(int))));

	/* find mean, sd, max and min of lb,
	check that parameter values are valid */
	sigsq = 0; sd = sigsq;
	lmax = 0.0; lmin = 0.0; mean = 0.0;
	for (j = 0; j < r; j++)
	{
		lj = lb[j];
		sd = sd + square(lj) * 2;
		mean += lj;
		if (lmax < lj) lmax = lj; else if (lmin > lj) lmin = lj;
	}
	if (sd == 0.0)
	{
		qfval = (c > 0.0) ? 1.0 : 0.0; goto  endofproc;
	}

	sd = std::sqrt(sd);
	almx = (lmax < -lmin) ? -lmin : lmax;

	/* starting values for findu, ctff */
	utx = 16.0 / sd;  up = 4.5 / sd;  un = -up;
	/* truncation point with no convergence factor */
	findu(&utx, .5 * acc);
	/* does convergence factor help */
	if (c != 0.0 && (almx > 0.07 * sd))
	{
		tausq = .25 * acc / cfe(c);
		if (fail) fail = FALSE;
		else if (truncation(utx, tausq) < .2 * acc)
		{
			sigsq = sigsq + tausq;
			findu(&utx, .25 * acc);
		}
	}
	acc *= 0.5;

	/* find RANGE of distribution, quit if outside this */
l1:
	d1 = ctff(acc, &up) - c;
	if (d1 < 0.0) { qfval = 1.0; goto endofproc; }
	d2 = c - ctff(acc, &un);
	if (d2 < 0.0) { qfval = 0.0; goto endofproc; }
	/* find integration interval */
	intv = 2.0 * pi / ((d1 > d2) ? d1 : d2);
	/* calculate number of terms required for main and
	auxillary integrations */
	xnt = utx / intv;  xntm = 3.0 / std::sqrt(acc);
	if (xnt > xntm * 1.5)
	{
		/* parameters for auxillary integration */
        ntm = (int)std::floor(xntm + 0.5);
		intv1 = utx / ntm;  x = 2.0 * pi / intv1;
		if (x <= fabs(c)) goto l2;
		/* calculate convergence factor */
		tausq = .33 * acc / (1.1 * (cfe(c - x) + cfe(c + x)));
		if (fail) goto l2;
		acc *= .67;
		/* auxillary integration */
		integrate(ntm, intv1, tausq, FALSE);
		xlim = xlim - xntm;  sigsq = sigsq + tausq;
		/* find truncation point with new convergence factor */
		findu(&utx, .25 * acc);  acc *= 0.75;
		goto l1;
	}

	/* main integration */
l2:
    nt =  static_cast<int>(std::floor(xnt + 0.5));
	integrate(nt, intv, 0.0, TRUE);
	qfval = 0.5 - intl;

	/* test whether round-off error could be significant
	allow for radix 8 or 16 machines */
	up = ersm; x = up + acc / 10.0;

endofproc:
    free((char*)th);
	return std::abs(1 - qfval);
}
