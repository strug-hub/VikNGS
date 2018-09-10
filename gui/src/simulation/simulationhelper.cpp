#include "Simulation.h"

/*
Produces a VectorXd holding the group ID for all samples.

@param simReq Provides information about sample size and group ID
@return VectorXd for group ID
*/
VectorXi simulateG(SimulationRequest& simReq) {

    VectorXi G(simReq.maxSize());
    int index = 0;

    for(int i = 0; i < simReq.steps; i++)
        for (SimulationRequestGroup srg : simReq.groups){
            int n = srg.getIncreaseSize(i);
            for (int j = 0; j < n; j++){
                G[index] = srg.index;
                index++;
            }
        }

    return G;
}

/*
Produces a VectorXd of phenotypes for each sample.

@param simReq Provides information about group phenotype and step sample size.
@return VectorXd of phenotypes.
*/
VectorXd simulateY(SimulationRequest& simReq) {

    VectorXd Y(simReq.maxSize());
    int index = 0;

    for(int i = 0; i < simReq.steps; i++)
        for (SimulationRequestGroup srg : simReq.groups){
            int n = srg.getIncreaseSize(i);
            for (int j = 0; j < n; j++){
                Y[index] = srg.generatePhenotype();
                index++;
            }
        }

    return Y;
}

/*
Produces a set of matrices (one per step) of true genotypes using minor allele frequency and odds ratio. Adjusts odds ratio with respect to collapsed region.

@param simReq Provides information about which groups are case/control and sample size.
@param oddsRatio Odds ratio of being affected given harmful variant.
@param maf Vector of minor allele frequencies for each variant.

@return MatrixXd of true genotypes (nsnp x nsamp).
*/
MatrixXd simulateX(SimulationRequest& simReq, double oddsRatio, VectorXd& maf) {
    int nsnp = maf.size();

    MatrixXd X(simReq.maxSize(), nsnp);
    std::vector<double> oddsRatios(static_cast<size_t>(nsnp));
    oddsRatios[0] = oddsRatio;

    //decay odds ratio outward from causitive variant
    if(nsnp > 1){
        int effectIndex = randomInt(0, nsnp-1);
        oddsRatios[static_cast<size_t>(effectIndex)] = oddsRatio;

        double lastOR = oddsRatio;
        for(int j = effectIndex-1; j >= 0; j--){
            lastOR = randomDouble(std::min(1.0, lastOR), std::max(1.0, lastOR));
            oddsRatios[static_cast<size_t>(j)] = lastOR;
        }
        lastOR = oddsRatio;
        for(int j = effectIndex+1; j < nsnp; j++){
            lastOR = randomDouble(std::min(1.0, lastOR), std::max(1.0, lastOR));
            oddsRatios[static_cast<size_t>(j)] = lastOR;
        }
    }

    for (int h = 0; h < nsnp; h++) {

            double OR = oddsRatios[static_cast<size_t>(h)];

            double p_x0_y0 = (1 - maf[h]) * (1 - maf[h]);
            double p_x1_y0 = 2 * maf[h] * (1 - maf[h]);

            double beta0 = 1 / (pow(maf[h] * (1 - OR) - 1, 2));

            double p_x0_y1 = beta0 * p_x0_y0;
            double p_x1_y1 = (beta0 * OR) * p_x1_y0;

        generate:
            int index = 0;
            bool ok = false;

            for(int i = 0; i < simReq.steps; i++){
                for (SimulationRequestGroup srg : simReq.groups){
                    int n = srg.getIncreaseSize(i);

                    if(srg.isCase)
                        for (int j = 0; j < n; j++){
                            X(index, h) = generateGenotype(p_x0_y1, p_x1_y1);
                            ok = ok || abs(X(index, h) - X(0, h)) > 1e-4;
                            index++;
                        }
                    else
                        for (int j = 0; j < n; j++){
                            X(index, h) = generateGenotype(p_x0_y0, p_x1_y0);
                            ok = ok || abs(X(index, h) - X(0, h)) > 1e-4;
                            index++;
                        }
                }

                if(!ok)
                    goto generate;

                double sum = 0;
                for(int p = 0; p< index; p++)
                    sum += X(p, h);

                if(sum < 1e-4)
                    std::cout << "why";
            }
        }

    return X;
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

VectorXi generateReadDepths(SimulationRequest& simReq){

    VectorXi readDepths(simReq.maxSize());

    int index = 0;
    for(int i = 0; i < simReq.steps; i++)
        for (SimulationRequestGroup srg : simReq.groups){
            int n = srg.getIncreaseSize(i);
            for (int j = 0; j < n; j++){
                readDepths[index] = static_cast<int>(std::round(randomNormal(srg.meanDepth, srg.sdDepth)));
                readDepths[index] = std::max(1, readDepths[index]);
                index++;
            }
        }

    return readDepths;
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

std::vector<std::vector<int>> generateBaseCalls(SimulationRequest& simReq, VectorXd& X, VectorXi& readDepths){

    std::vector<std::vector<int>> calls(static_cast<size_t>(simReq.maxSize()));

    int index = 0;
    for(int i = 0; i < simReq.steps; i++)
        for (SimulationRequestGroup srg : simReq.groups){
            int n = srg.getIncreaseSize(i);
            for (int j = 0; j < n; j++){
                calls[static_cast<size_t>(index)] =
                        baseCalls(static_cast<int>(X[index]), srg.errorRate, readDepths[index]);
                index++;
            }
        }

    return calls;
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
Vector3d calculateLikelihood(std::vector<int> &bases, double error) {

    Vector3d gl;

    int refCalls = bases[0];
    int altCalls = bases[1];
    int errorCalls = bases[2] + bases[3];

    gl[0] = std::pow(1-error, refCalls) * std::pow(error/3, altCalls + errorCalls);
    gl[1] = std::pow(0.5 * (1-error + (error/3)), refCalls + altCalls) * std::pow(error/3, errorCalls);
    gl[2] = std::pow(1-error, altCalls) * std::pow(error/3, refCalls + errorCalls);

    return gl;
}

std::vector<Vector3d> generateLikelihoods(SimulationRequest& simReq, std::vector<std::vector<int>>& baseCalls){

    std::vector<Vector3d> likelihoods(static_cast<size_t>(simReq.maxSize()));

    size_t index = 0;
    for(int i = 0; i < simReq.steps; i++)
        for (SimulationRequestGroup srg : simReq.groups){
            int n = srg.getIncreaseSize(i);
            for (int j = 0; j < n; j++){
                likelihoods[index] = calculateLikelihood(baseCalls[index], srg.errorRate);
                index++;
            }
        }

    return likelihoods;
}

Variant randomVariant(){

    int refalt = randomInt(0,11);
    int rchr = randomInt(0,23);
    std::string ref;
    std::string alt;
    std::string chr;
    int pos = 0;

    if(refalt==0) ref = "A";
    else if(refalt==1) ref = "A";
    else if(refalt==2) ref = "A";
    else if(refalt==3) ref = "T";
    else if(refalt==4) ref = "T";
    else if(refalt==5) ref = "T";
    else if(refalt==6) ref = "C";
    else if(refalt==7) ref = "C";
    else if(refalt==8) ref = "C";
    else if(refalt==9) ref = "G";
    else if(refalt==10) ref = "G";
    else if(refalt==11) ref = "G";

    if(refalt==0) alt = "T";
    else if(refalt==1) alt = "G";
    else if(refalt==2) alt = "C";
    else if(refalt==3) alt = "A";
    else if(refalt==4) alt = "G";
    else if(refalt==5) alt = "C";
    else if(refalt==6) alt = "T";
    else if(refalt==7) alt = "A";
    else if(refalt==8) alt = "G";
    else if(refalt==9) alt = "C";
    else if(refalt==10) alt = "A";
    else if(refalt==11) alt = "T";

    if(rchr==0){ chr = "1"; pos = randomInt(0, 248956422);}
    if(rchr==1){ chr = "2"; pos = randomInt(0, 242193529);}
    if(rchr==2){ chr = "3"; pos = randomInt(0, 198295559);}
    if(rchr==3){ chr = "4"; pos = randomInt(0, 190214555);}
    if(rchr==4){ chr = "5"; pos = randomInt(0, 181538259);}
    if(rchr==5){ chr = "6"; pos = randomInt(0, 170805979);}
    if(rchr==6){ chr = "7"; pos = randomInt(0, 159345973);}
    if(rchr==7){ chr = "8"; pos = randomInt(0, 145138636);}
    if(rchr==8){ chr = "9"; pos = randomInt(0, 138394717);}
    if(rchr==9){ chr = "10"; pos = randomInt(0, 133797422);}
    if(rchr==10){ chr = "11"; pos = randomInt(0, 135086622);}
    if(rchr==11){ chr = "12"; pos = randomInt(0, 133275309);}
    if(rchr==12){ chr = "13"; pos = randomInt(0, 114364328);}
    if(rchr==13){ chr = "14"; pos = randomInt(0, 107043718);}
    if(rchr==14){ chr = "15"; pos = randomInt(0, 101991189);}
    if(rchr==15){ chr = "16"; pos = randomInt(0, 90338345);}
    if(rchr==16){ chr = "17"; pos = randomInt(0, 83257441);}
    if(rchr==17){ chr = "18"; pos = randomInt(0, 80373285);}
    if(rchr==18){ chr = "19"; pos = randomInt(0, 58617616);}
    if(rchr==19){ chr = "20"; pos = randomInt(0, 64444167);}
    if(rchr==20){ chr = "21"; pos = randomInt(0, 46709983);}
    if(rchr==21){ chr = "22"; pos = randomInt(0, 50818468);}
    if(rchr==22){ chr = "X"; pos = randomInt(0, 156040895);}
    if(rchr==23){ chr = "Y"; pos = randomInt(0, 57227415);}

    return Variant(chr, pos, "simulated", ref, alt);
}

