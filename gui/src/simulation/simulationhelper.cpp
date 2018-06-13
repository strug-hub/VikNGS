#include "simulation.h"


Variant generateSeqData(VectorXd x, VectorXd g, std::map<int, SimulationRequestGroup> group, Variant &variant) {
    double mean, sd, error;
    int rd;

    MatrixXd M(x.rows(), 3);
    std::vector<MatrixXd> mM;
    std::vector<GenotypeLikelihood> likelihoods;
    VectorXd genotypeCalls(x.rows());

    for (int i = 0; i < x.rows(); i++) {

        mean = group[g[i]].meanDepth;
        sd = group[g[i]].sdDepth;
        error = group[g[i]].errorRate;

        rd = std::round(randomNormal(mean, sd));
        rd = std::max(rd, 1);

        int k = x[i];

        std::vector<char> trueGenotype;

        //just assume that T is reference and C is alt for simplicity
        if (k == 0) {
            trueGenotype.push_back('T');
            trueGenotype.push_back('T');
        }
        if (k == 1) {
            trueGenotype.push_back('C');
            trueGenotype.push_back('T');
        }
        if (k == 2) {
            trueGenotype.push_back('C');
            trueGenotype.push_back('C');
        }

        std::vector<char> bases = baseCall(trueGenotype, error, rd);
        Vector3d lh = calculateLikelihood(bases, error);
        M.row(i) =  lh;

        GenotypeLikelihood l;
        l.L00 = lh[0];
        l.L01 = lh[1];
        l.L11 = lh[2];
        likelihoods.push_back(l);

        if(lh[0] > lh[1])
            if (lh[0] > lh[2])
                genotypeCalls[i] = 0;
            else
                genotypeCalls[i] = 2;
        else if (lh[1] > lh[2])
            genotypeCalls[i] = 1;
        else
            genotypeCalls[i] = 2;

        mM.push_back(calculateLikelihood2(bases, error));
    }

    VectorXd p = calcEM(M);
    VectorXd EG = calculateExpectedGenotypes(mM, p);

    variant.likelihood = likelihoods;
    variant.P = p;

    variant.expectedGenotype = EG;
    variant.trueGenotype = x;
    variant.genotypeCalls = genotypeCalls;

    return variant;
}


/*
Produces a vector of simulated minor allele frequencies (uniformly random between min and max).

@param nsnp Number of MAF to produce.
@param min Minimum value for MAF.
@param max Maximum value for MAF.

@return Vector of minor allele frequencies.
*/
VectorXd simulateMinorAlleleFrequency(int nsnp, double min, double max) {
	VectorXd maf(nsnp);

	for (int i = 0; i < nsnp; i++) 
		maf[i] = randomDouble(min, max);

	return maf;
}

double inline generateGenotype(double prob_y, double prob_x0, double prob_x1) {

    //prob(x=0, y=1) + prob(x=1, y=1) + prob(x=2, y=1) = prob(y=1)
	double rand = randomDouble(0, prob_y);

	if (rand > prob_x0)
		if (rand < prob_x0 + prob_x1)
			return 1;
		else
			return 2;
	else
		return 0;
}

/*
Produces a matrix of genotypes using minor allele frequency and odds ratio.

@param simReq Provides information about which groups are case/control.
@param oddsRatio Odds ratio of being affected given harmful variant.
@param maf Vector of minor allele frequencies for each variant.

@return Matrix of expected genotypes.
*/
MatrixXd simulateX(SimulationRequest simReq, double oddsRatio, VectorXd maf) {
	int nsnp = maf.size();
    int ncont = simReq.ncontrol();
    int ncase = simReq.ncase();
    int nsamp = ncont + ncase;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> X(nsamp, nsnp);

    double p_y1 = ncase / (double)nsamp;
    double p_y0 = 1 - p_y1 ;

	for (int i = 0; i < nsnp; i++) {

		double beta0 = ncase / (ncont * pow(maf[i] * (1 - oddsRatio) - 1, 2));

		double odds_y1_x0 = beta0;
        //double odds_y1_x1 = beta0 * oddsRatio;
        //double odds_y1_x2 = beta0 * oddsRatio * oddsRatio;


        double p_x0_y0 = ncont / (double) nsamp * pow(1 - maf[i], 2);
        double p_x1_y0 = 2 * maf[i] * (1 - maf[i]) * ncont / (double) nsamp;
		double p_x0_y1 = odds_y1_x0 * p_x0_y0;
		double p_x1_y1 = odds_y1_x0 * oddsRatio * p_x1_y0;

        int index = 0;
        for (int j = 0; j < simReq.groups.size(); j++){

            if(simReq.groups[j].isCase){
                for (int k = 0; k < simReq.groups[j].n; k++){
                    X(index, i) = generateGenotype(p_y1, p_x0_y1, p_x1_y1);
                    index++;
                }
            }
            else{
                for (int k = 0; k < simReq.groups[j].n; k++){
                    X(index, i) = generateGenotype(p_y0, p_x0_y0, p_x1_y0);
                    index++;
                }
            }
        }
    }

	return X;
}

/*
Produces a matrix of genotypes using minor allele frequency and odds ratio. Adjusts odds ratio with respect to collapsed region.

@param simReq Provides information about which groups are case/control.
@param oddsRatio Odds ratio of being affected given harmful variant.
@param maf Vector of minor allele frequencies for each variant.

@param collapse Number of variants per collapsed region.

@return Matrix of expected genotypes.
*/
MatrixXd simulateX(SimulationRequest simReq, double oddsRatio, VectorXd maf, int collapse) {
    int nsnp = maf.size();
    int ncont = simReq.ncontrol();
    int ncase = simReq.ncase();
    int nsamp = ncont + ncase;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> X(nsamp, nsnp);

    double p_y1 = ncase / (double)nsamp;
    double p_y0 = 1 - p_y1 ;

    std::vector<double> oddsRatios(collapse);

    for (int i = 0; i < nsnp; i++) {

        if(i % collapse == 0){
            int effectIndex = randomInt(0, collapse-1);
            oddsRatios[effectIndex] = oddsRatio;

            double lastOR = oddsRatio;
            for(int j = effectIndex-1; j >= 0; j--){
                lastOR = randomDouble(std::min(1.0, lastOR), std::max(1.0, lastOR));
                oddsRatios[j] = lastOR;
            }
            lastOR = oddsRatio;
            for(int j = effectIndex+1; j < collapse; j++){
                lastOR = randomDouble(std::min(1.0, lastOR), std::max(1.0, lastOR));
                oddsRatios[j] = lastOR;
            }
        }

        double odds = oddsRatios[i%collapse];

        double beta0 = ncase / (ncont * pow(maf[i] * (1 - odds) - 1, 2));

        double odds_y1_x0 = beta0;
        //double odds_y1_x1 = beta0 * oddsRatio;
        //double odds_y1_x2 = beta0 * oddsRatio * oddsRatio;


        double p_x0_y0 = ncont / (double) nsamp * pow(1 - maf[i], 2);
        double p_x1_y0 = 2 * maf[i] * (1 - maf[i]) * ncont / (double) nsamp;
        double p_x0_y1 = odds_y1_x0 * p_x0_y0;
        double p_x1_y1 = odds_y1_x0 * odds * p_x1_y0;

        int index = 0;
        for (int j = 0; j < simReq.groups.size(); j++){

            if(simReq.groups[j].isCase){
                for (int k = 0; k < simReq.groups[j].n; k++){
                    X(index, i) = generateGenotype(p_y1, p_x0_y1, p_x1_y1);
                    index++;
                }
            }
            else{
                for (int k = 0; k < simReq.groups[j].n; k++){
                    X(index, i) = generateGenotype(p_y0, p_x0_y0, p_x1_y0);
                    index++;
                }
            }
        }
    }

    return X;
}

/*
Produces a vector of case/control status.

@param simReq Provides information about which groups are case/control.

@return Vector of case/control status.
*/
VectorXd simulateY(SimulationRequest simReq) {
    int ncont = simReq.ncontrol();
    int ncase = simReq.ncase();
    int nsamp = ncont + ncase;

    VectorXd Y(nsamp);

    int index = 0;
    for (int j = 0; j < simReq.groups.size(); j++){

        if(simReq.groups[j].isCase){
            for (int k = 0; k < simReq.groups[j].n; k++){
                Y[index] = 1;
                index++;
            }
        }
        else{
            for (int k = 0; k < simReq.groups[j].n; k++){
                Y[index] = 0;
                index++;
            }
        }
    }
	return Y;
}

/*
Produces a vector of characters representing the called bases given a
true genotype, read depth and sequencing platform error rate.

@param trueGenotype The base corresponding to the true genotype.
@param error The probability of error for each read.
@param readDepth Number of reads.

@return Vector of called bases.
*/
std::vector<char> baseCall(std::vector<char> trueGenotype, double error, int readDepth) {
	
	std::vector<char> bases;
	double r;
	char trueBase;
	std::vector<char> errorBase;

    double e3 = error/3.0;
    double e2 = 2.0 * e3;

	for (int i = 0; i < readDepth; i++) {
		trueBase = trueGenotype[randomInt(0, 1)];
		if (trueBase == 'C')
			errorBase = {'A', 'G', 'T'};
		else if (trueBase == 'T')
			errorBase = { 'A', 'G', 'C' };
		else if (trueBase == 'G')
			errorBase = { 'A', 'T', 'C' };
		else if (trueBase == 'A')
			errorBase = { 'G', 'T', 'C' };

		r = randomDouble(0, 1);

		if (r <= e3)
			bases.push_back(errorBase[0]);
        else if (r <= e2)
			bases.push_back(errorBase[1]);
        else if (r <= error)
			bases.push_back(errorBase[2]);
		else
			bases.push_back(trueBase);
	}

	return bases;
}

/*
Produces a vector of 3 probabilities as follows:
0: probability of calling every base in bases in a homozygous individual (2 major alleles).
1: probability of calling every base in bases in a heterzygous individual (1 minor allele).
2: probability of calling every base in bases in a homozygous individual (2 minor alleles).

@param bases Vector of called bases.
@param error Probability of error for each read.

@return Vector of 3 probabilities.
*/
VectorXd calculateLikelihood(std::vector<char> &bases, double error) {
	
	VectorXd ll(3);
	ll[0] = 1;
	ll[1] = 1;
	ll[2] = 1;

	for (int i = 0; i < bases.size(); i++) {
        ll[0] = ll[0] * pSingle(bases[i], 'T', 'T', error);
        ll[1] = ll[1] * pSingle(bases[i], 'C', 'T', error);
        ll[2] = ll[2] * pSingle(bases[i], 'C', 'C', error);
	}
	return ll;
}

/*
Produces a matrix where each column has 3 probabilities as follows:
0: probability of calling base in a homozygous individual (2 major alleles).
1: probability of calling base in a heterzygous individual (1 minor allele).
2: probability of calling base in a homozygous individual (2 minor alleles).

Each row of the vector is a base call corresponding to a different read.

@param bases Vector of called bases.
@param error Probability of error for each read.

@return nx3 matrix of probabilities (n = read depth).
*/
MatrixXd calculateLikelihood2(std::vector<char> &bases, double error) {
	MatrixXd ll(bases.size(), 3);

	for (int i = 0; i < bases.size(); i++) {
        ll(i, 0) = pSingle(bases[i], 'T', 'T', error);
        ll(i, 1) = pSingle(bases[i], 'C', 'T', error);
        ll(i, 2) = pSingle(bases[i], 'C', 'C', error);
	}
	return ll;
}

/*
Calculates the probability of calling a base given the true genotype and base call error rate.

@param base Called base.
@param true1 True base 1.
@param true2 True base 2 (same as true1 = homozygote, different = heterozygote).
@param error Probability of a base call error.

@return Probability of calling base.
*/
double pSingle(char base, char true1, char true2, double error) {
	
	double p1, p2;

	if (base == true1)
		p1 = 1 - error;
	else
        p1 = error / 3.0;
	
	if (base == true2)
		p2 = 1 - error;
	else
        p2 = error / 3.0;

	return 0.5 * p1 + 0.5 * p2;
}

/*
Uses EM algorithm to estimate the genotype frequencies in the sample.

@param  M A vector with genotype likelihoods (calculated by function calculateLikelihood).
@return A vector with probability of 0, 1 or 2 minor alleles.
*/
VectorXd calcEM(MatrixXd M) {

	double p = 0.15;
	double q = 0.15;
	double qn = 1;
	double pn = 0;
	double dn = 0;
	double d = 0;

	int k = 0;

	while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

		d = 1 - p - q;
		Vector3d v = { p, q, d };
		VectorXd pD = M * v;
		VectorXd Ep = M.col(0).array() * (p / pD.array()).array();
		VectorXd Eq = M.col(1).array() * (q / pD.array()).array();
		pn = p;
		qn = q;
		dn = 1 - q - p;
		p = Ep.sum() / Ep.rows();
		q = Eq.sum() / Eq.rows();

		k++;
		if (k == 1000)
			break;
	}

	VectorXd freq(3);
	freq[0] = std::max(0.0, p);
	freq[1] = std::max(0.0, q);
	freq[2] = std::max(0.0, 1 - p - q);

	return freq;
}

/*
Generates the expected probabilities of the genotypes E(G_ij | D_ij).

@param  M A vector with genotype likelihoods (calculated by function calculateLikelihood2).
@param  p Population frequency (calculated by function calcEM).

@return Vector of expected genotypes.
*/
VectorXd calculateExpectedGenotypes(std::vector<MatrixXd> M, VectorXd p){
		
	VectorXd m(3);
	Vector3d g = { 0, 1, 2 };
	VectorXd EG( M.size());
	double pm;

	for (int i = 0; i < M.size(); i++) {
		for (int j = 0; j < 3; j++) {
			double L = 1;

			for (int k = 0; k < M[i].rows(); k++)
				L = L * M[i](k, j);

			m[j] = L * p[j];			
		}

		pm = (m.array() / m.sum() * g.array()).sum();
		pm = std::max(0.0, pm);
		pm = std::min(2.0, pm);
		EG[i] = pm;
	}

	return EG;
}


Variant randomVariant(){

    Variant v;

    int rref = randomInt(0,3);
    int ralt = randomInt(0,3);
    int rchr = randomInt(0,23);

    if(rref==0) v.ref = "T";
    if(rref==1) v.ref = "C";
    if(rref==2) v.ref = "A";
    if(rref==3) v.ref = "G";

    if(ralt==0) v.alt = "T";
    if(ralt==1) v.alt = "C";
    if(ralt==2) v.alt = "A";
    if(ralt==3) v.alt = "G";

    if(rchr==0){ v.chr = "1"; v.pos = randomInt(0, 248956422);}
    if(rchr==1){ v.chr = "2"; v.pos = randomInt(0, 242193529);}
    if(rchr==2){ v.chr = "3"; v.pos = randomInt(0, 198295559);}
    if(rchr==3){ v.chr = "4"; v.pos = randomInt(0, 190214555);}
    if(rchr==4){ v.chr = "5"; v.pos = randomInt(0, 181538259);}
    if(rchr==5){ v.chr = "6"; v.pos = randomInt(0, 170805979);}
    if(rchr==6){ v.chr = "7"; v.pos = randomInt(0, 159345973);}
    if(rchr==7){ v.chr = "8"; v.pos = randomInt(0, 145138636);}
    if(rchr==8){ v.chr = "9"; v.pos = randomInt(0, 138394717);}
    if(rchr==9){ v.chr = "10"; v.pos = randomInt(0, 133797422);}
    if(rchr==10){ v.chr = "11"; v.pos = randomInt(0, 135086622);}
    if(rchr==11){ v.chr = "12"; v.pos = randomInt(0, 133275309);}
    if(rchr==12){ v.chr = "13"; v.pos = randomInt(0, 114364328);}
    if(rchr==13){ v.chr = "14"; v.pos = randomInt(0, 107043718);}
    if(rchr==14){ v.chr = "15"; v.pos = randomInt(0, 101991189);}
    if(rchr==15){ v.chr = "16"; v.pos = randomInt(0, 90338345);}
    if(rchr==16){ v.chr = "17"; v.pos = randomInt(0, 83257441);}
    if(rchr==17){ v.chr = "18"; v.pos = randomInt(0, 80373285);}
    if(rchr==18){ v.chr = "19"; v.pos = randomInt(0, 58617616);}
    if(rchr==19){ v.chr = "20"; v.pos = randomInt(0, 64444167);}
    if(rchr==20){ v.chr = "21"; v.pos = randomInt(0, 46709983);}
    if(rchr==21){ v.chr = "22"; v.pos = randomInt(0, 50818468);}
    if(rchr==22){ v.chr = "X"; v.pos = randomInt(0, 156040895);}
    if(rchr==23){ v.chr = "Y"; v.pos = randomInt(0, 57227415);}

    return v;
}

