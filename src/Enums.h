#pragma once
#include <string>

enum class Statistic { NONE, COMMON, CAST, SKAT, CALPHA };
enum class Variance { NONE, REGULAR, RVS, RVSFALSE };
enum class Genotype { NONE, EXPECTED, TRUE, CALL, VCF_CALL };
enum class Family { NONE, NORMAL, BINOMIAL };
enum class Depth { HIGH, LOW };
enum class CollapseType { COLLAPSE_K, COLLAPSE_GENE, COLLAPSE_EXON, NONE };
enum class Filter { VALID, INVALID, IGNORE, NOT_SNP, NO_PASS, MISSING_DATA, NO_VARIATION, MAF };

inline bool isRare(Statistic s) {return s == Statistic::CAST || s == Statistic::SKAT || s == Statistic::CALPHA;}

inline std::string genotypeToString(Genotype g){
    switch(g) {
        case Genotype::EXPECTED: return "Expected";
        case Genotype::TRUE: return "True";
        case Genotype::CALL: return "Call";
        case Genotype::VCF_CALL: return "VCF Call";
        default: return "_";
    }
}

struct Test {
private:
    Genotype genotype;
    Statistic statistic;
    Variance variance;
    int nboot = -1;
    int nsamples = -1;

public:
    Test(Genotype g, Statistic s, Variance v, int nbootstrap=-1, int nsamples=-1) : genotype(g), statistic(s), variance(v){
        this->nsamples = nsamples;
        this->nboot = nbootstrap;
    }

    inline bool isRareTest(){ return isRare(statistic); }
    inline bool isCommonTest(){ return statistic == Statistic::COMMON; }
    inline bool isRVS(){ return variance == Variance::RVS || variance == Variance::RVSFALSE; }
    inline bool isRVSFalse(){ return variance == Variance::RVSFALSE; }
    inline void setRVSFalse(){ if(variance == Variance::RVS) variance = Variance::RVSFALSE; }
    inline void setSampleSize(int size){ nsamples=size; }
    inline void setGenotype(Genotype gt){ genotype=gt; }

    inline bool needVCFCalls(){ return genotype == Genotype::VCF_CALL; }
    inline bool needGenotypeCalls(){ return genotype == Genotype::CALL; }
    inline bool needExpectedGenotypes(){ return genotype == Genotype::EXPECTED; }
    inline bool isExpectedGenotypes(){ return genotype == Genotype::EXPECTED; }
    inline bool isBootstrap(){ return nboot>1; }

    inline Genotype getGenotype(){ return genotype; }
    inline Statistic getStatistic(){ return statistic; }
    inline Variance getVariance(){ return variance; }

    inline int getSampleSize(){ return nsamples; }
    inline int getbootstrapSize(){ return nboot; }

    inline std::string toString(){
        std::string part1, part2;
        switch(genotype)
        {
            case Genotype::EXPECTED :
                part1 = "Genotype Likelihoods";
                if(isRVS()) part1 = "RVS " + part1;
                break;
            case Genotype::CALL : part1 = "Genotype Calls"; break;
            case Genotype::VCF_CALL : part1 = "VCF Calls"; break;
            case Genotype::TRUE : part1 = "True Genotypes"; break;
            default: part1 = "???"; break;
        }
        switch(statistic)
        {
            case Statistic::COMMON : part2 = "common"; break;
            case Statistic::CAST : part2 = "CAST"; break;
            case Statistic::SKAT : part2 = "SKAT"; break;
            case Statistic::CALPHA : part2 = "Calpha"; break;
            default: part2 = "???"; break;
        }

        return part1 + " - " + part2;
    }

    inline std::string toShortString(){
        std::string part1, part2;
        switch(genotype)
        {
            case Genotype::EXPECTED :
                if(isRVS()) part1 = "RVS"; else part1 = "No RVS";
                break;
            case Genotype::CALL : part1 = "Call"; break;
            case Genotype::VCF_CALL : part1 = "VCF"; break;
            case Genotype::TRUE : part1 = "True"; break;
            default: part1 = "???"; break;
        }
        switch(statistic)
        {
            case Statistic::COMMON : part2 = "common"; break;
            case Statistic::CAST : part2 = "CAST"; break;
            case Statistic::SKAT : part2 = "SKAT"; break;
            case Statistic::CALPHA : part2 = "Calpha"; break;

            default: return "???";
        }

        return part1 + " " + part2;
    }
};
