#include "simulation.h"

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

/*
Produces a set of matrices (one per step) of true genotypes using minor allele frequency and odds ratio. Adjusts odds ratio with respect to collapsed region.

@param simReq Provides information about which groups are case/control and sample size.
@param oddsRatio Odds ratio of being affected given harmful variant.
@param maf Vector of minor allele frequencies for each variant.
@param collapse Decay odds ratio outward from causative, one causative variant every collapse.

@return Vector of MatrixXd of true genotypes (nsnp x nsamp).
*/
std::vector<MatrixXd> simulateX(SimulationRequest& simReq, double oddsRatio, VectorXd& maf, int collapse) {
    int nsnp = maf.size();

    std::vector<MatrixXd> X;
    //for collapsing only
    std::vector<double> oddsRatios(collapse);
    double OR = oddsRatio;

    for(int i = 0; i < simReq.steps; i++){

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> x(nsnp, simReq.stepIncrementSize(i));

        for (int h = 0; h < nsnp; h++) {

            if(collapse > 1){

                if(h % collapse == 0){
                    //define causative variant in collapsed region
                    int effectIndex = randomInt(0, collapse-1);
                    oddsRatios[effectIndex] = oddsRatio;

                    //decay odds ratio outward from causitive variant
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

                 OR = oddsRatios[i%collapse];
            }

            double p_x0_y0 = (1 - maf[h]) * (1 - maf[h]);
            double p_x1_y0 = 2 * maf[h] * (1 - maf[h]);

            double beta0 = 1 / (pow(maf[h] * (1 - OR) - 1, 2));

            double p_x0_y1 = beta0 * p_x0_y0;
            double p_x1_y1 = (beta0 * OR) * p_x1_y0;

            int index = 0;
            for (SimulationRequestGroup srg : simReq.groups) {

                if(srg.isCase){
                    for (int j = 0; j < srg.getStepSize(i); j++){
                        x(h, index) = generateGenotype(p_x0_y1, p_x1_y1);
                        index++;
                    }
                }
                else{
                    for (int j = 0; j < srg.getStepSize(i); j++){
                        x(h, index) = generateGenotype(p_x0_y0, p_x1_y0);
                        index++;
                    }
                }
            }
        }
        X.emplace_back(x);
    }

    return X;
}

/*
Produces a set of VectorXd holding the group ID for every step of the simulation.

@param simReq Provides information about sample size and group ID

@return Vector of VectorXd for group ID
*/
std::vector<VectorXd> simulateG(SimulationRequest& simReq) {

    std::vector<VectorXd> G;

    for(int i = 0; i < simReq.steps; i++){

        VectorXd g(simReq.stepIncrementSize(i));
        int counter = 0;
        for (SimulationRequestGroup srg : simReq.groups) {
            for(int j = counter; j < counter + srg.getStepSize(i); j++)
                g[j] = srg.index;
            counter += srg.getStepSize(i);

        }

        G.emplace_back(g);
    }

    return G;
}

/*
Produces a set of VectorXd for phenotype of each sampe for every step.

@param simReq Provides information about group phenotype and step sample size.
@return Vector of VectorXd of case/control status.
*/
std::vector<VectorXd> simulateY(SimulationRequest& simReq) {

    std::vector<VectorXd> Y;

    for(int i = 0; i < simReq.steps; i++){

        VectorXd y(simReq.stepIncrementSize(i));
        int index = 0;

        for (SimulationRequestGroup srg : simReq.groups) {

            for (int j = 0; j < srg.getStepSize(i); j++){
                y[index] = srg.generatePhenotype();
                index++;
            }
        }
        Y.push_back(y);
    }

    return Y;
}

double inline generateGenotype(double prob_x0, double prob_x1) {

    double rand = randomDouble(0, 1);

    if (rand > prob_x0)
        if (rand < prob_x0 + prob_x1)
            return 1;
        else
            return 2;
    else
        return 0;
}

/*
Produces a vector of characters representing the called bases given a
true genotype, read depth and sequencing platform error rate.

@param trueGenotype The true genotype (0, 1 or 2).
@param error The probability of error for each read.
@param readDepth Number of reads.

@return Count vector of called bases [ref, alt, error1, error2]
*/
std::vector<int> baseCalls(int trueGenotype, double error, int readDepth) {

    std::vector<int> bases;
    int correct = randomBinomial(readDepth, 1-error);
    int mistake = readDepth - correct;

    int errorBase1=0;
    int errorBase2=0;
    int errorBase3=0;

    //split errors uniformly in 3
    if(mistake > 0){
        errorBase1 = randomInt(0, mistake);
        errorBase2 = randomInt(0, mistake);
        errorBase3 = mistake - (errorBase1+errorBase2);

        if(errorBase3 < 0){
            errorBase1 = mistake - errorBase1;
            errorBase2 = mistake - errorBase2;
            errorBase3 = - errorBase3;
        }
    }

    if(trueGenotype == 0){
        bases.push_back(correct);
        bases.push_back(errorBase3);
    }
    else if (trueGenotype == 1){

        int ref = randomBinomial(correct, 0.5);
        int alt = correct - ref;

        int luckyRef=0;
        int luckyAlt=0;

        if(errorBase3 > 0){
            luckyRef = randomInt(0, errorBase3);
            luckyAlt = errorBase3 - luckyRef;
        }

        bases.push_back(ref + luckyRef);
        bases.push_back(alt + luckyAlt);
    }
    else {
        bases.push_back(errorBase3);
        bases.push_back(correct);
    }


    bases.push_back(errorBase1);
    bases.push_back(errorBase2);

    return bases;
}

/*
Produces a vector of 3 probabilities as follows:
0: probability of calling bases in a homozygous individual (2 major alleles).
1: probability of calling bases in a heterozygous individual (1 minor allele).
2: probability of calling every base in bases in a homozygous individual (2 minor alleles).

@param bases Count vector of called bases [ref, alt, error1, error2].
@param error Probability of error for each read.

@return Vector of 3 probabilities.
*/
GenotypeLikelihood calculateLikelihood(std::vector<int> &bases, double error) {

    GenotypeLikelihood gl;

    int refCalls = bases[0];
    int altCalls = bases[1];
    int errorCalls = bases[2] + bases[3];

    gl.L00 = std::pow(1-error, refCalls) * std::pow(error/3, altCalls + errorCalls);
    gl.L01 = std::pow(0.5 * (1-error + (error/3)), refCalls + altCalls) * std::pow(error/3, errorCalls);
    gl.L11 = std::pow(1-error, altCalls) * std::pow(error/3, refCalls + errorCalls);

    return gl;
}

std::vector<GenotypeLikelihood> generateSeqData(VectorXd& x, SimulationRequestGroup& group) {
    double mean = group.meanDepth;
    double sd = group.sdDepth;
    double error = group.errorRate;
    int rd;

    std::vector<GenotypeLikelihood> likelihoods;

    for (int i = 0; i < x.rows(); i++) {

        rd = std::round(randomNormal(mean, sd));
        rd = std::max(rd, 1);
        std::vector<int> bases = baseCalls(x[i], error, rd);

        likelihoods.emplace_back(calculateLikelihood(bases, error));
    }

    return likelihoods;
}

/*
Generates the expected probabilities of the genotypes E(G_ij | D_ij).

@param  gl A vector of genotype likelihoods, one per individual.
@param  p Population genotype frequency.

@return Vector of expected genotypes.
*/
VectorXd calculateExpectedGenotypes(std::vector<GenotypeLikelihood>& gl, VectorXd& p){

    VectorXd m(3);
    Vector3d g = { 0, 1, 2 };
    VectorXd EG( gl.size());

    for (int i = 0; i < gl.size(); i++) {
        m[0] = gl[i].L00 * p[0];
        m[1] = gl[i].L01 * p[1];
        m[2] = gl[i].L11 * p[2];

        double msum = m[0] + m[1] + m[2];
        EG[i] = (m[1]/msum) + 2 * (m[2]/msum);
        EG[i] = std::max(0.0, EG[i]);
        EG[i] = std::min(2.0, EG[i]);

    }

    return EG;
}

/*
Generates the genotype calls using simple bayesian genotyper (maximum likelihood).

@param  gl A vector of genotype likelihoods, one per individual.
@param  p Population genotype frequency.

@return Vector of genotype calls.
*/
VectorXd calculateGenotypeCalls(std::vector<GenotypeLikelihood>& gl){

    VectorXd GC(gl.size());

    for (int i = 0; i < gl.size(); i++) {

        if(gl[i].L00 > gl[i].L01)
             if (gl[i].L00 > gl[i].L11)
                 GC[i] = 0;
             else
                 GC[i] = 2;
         else if (gl[i].L01 > gl[i].L11)
             GC[i] = 1;
         else
        GC[i] = 2;
    }

    return GC;
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

