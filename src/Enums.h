#pragma once
#include <string>

enum class Statistic { NONE, COMMON, CAST, SKAT, CALPHA };
enum class Variance { NONE, REGULAR, RVS, RVSFALSE };
enum class GenotypeSource { NONE, EXPECTED, TRUE, CALL, VCF_CALL };
enum class Family { NONE, NORMAL, BINOMIAL };
enum class Depth { HIGH, LOW };
enum class CollapseType { COLLAPSE_K, COLLAPSE_GENE, COLLAPSE_EXON, NONE };
enum class Filter { VALID, INVALID, IGNORE, NOT_SNP, NO_PASS, MISSING_DATA, NO_VARIATION, MAF };
//enum class Bootstrap { NONE, PERMUTE,  };

inline bool isRare(Statistic s) {return s == Statistic::CAST || s == Statistic::SKAT || s == Statistic::CALPHA;}

inline std::string genotypeToString(GenotypeSource g){
    switch(g) {
        case GenotypeSource::EXPECTED: return "Expected";
        case GenotypeSource::TRUE: return "True";
        case GenotypeSource::CALL: return "Call";
        case GenotypeSource::VCF_CALL: return "VCF Call";
        default: return "_";
    }
}

struct Test {
private:
    GenotypeSource genotype;
    Statistic statistic;
    Variance variance;
    int nboot = -1;
    int nsamples = -1;
    bool earlyStopping = false;

public:
    Test(GenotypeSource g, Statistic s, Variance v, int nbootstrap=-1, int nsamples=-1) : genotype(g), statistic(s), variance(v){
        this->nsamples = nsamples;
        this->nboot = nbootstrap;
    }

    inline bool isRareTest(){ return isRare(statistic); }
    inline bool isCommonTest(){ return statistic == Statistic::COMMON; }
    inline bool isRVS(){ return variance == Variance::RVS || variance == Variance::RVSFALSE; }
    inline bool isRVSFalse(){ return variance == Variance::RVSFALSE; }
    inline void setRVSFalse(){ if(variance == Variance::RVS) variance = Variance::RVSFALSE; }
    inline void setSampleSize(int size){ nsamples=size; }
    inline void setRegularVariance(){ variance = Variance::REGULAR; }
    inline void setGenotype(GenotypeSource gt){ genotype=gt; }
    inline void setEarlyStopping(bool value) { earlyStopping = value; }

    inline bool needVCFCalls(){ return genotype == GenotypeSource::VCF_CALL; }
    inline bool needGenotypeCalls(){ return genotype == GenotypeSource::CALL; }
    inline bool needExpectedGenotypes(){ return genotype == GenotypeSource::EXPECTED; }
    inline bool isExpectedGenotypes(){ return genotype == GenotypeSource::EXPECTED; }
    inline bool isBootstrap(){ return nboot>1; }

    inline GenotypeSource getGenotype(){ return genotype; }
    inline Statistic getStatistic(){ return statistic; }
    inline Variance getVariance(){ return variance; }

    inline int getSampleSize(){ return nsamples; }
    inline int getBootstrapSize(){ return nboot; }
    inline bool useEarlyStopping(){ return earlyStopping; }

    inline std::string toString(){
        std::string part1, part2;
        switch(genotype)
        {
            case GenotypeSource::EXPECTED :
                part1 = "Genotype Likelihoods";
                if(isRVS()) part1 = "RVS " + part1;
                break;
            case GenotypeSource::CALL : part1 = "Genotype Calls"; break;
            case GenotypeSource::VCF_CALL : part1 = "VCF Calls"; break;
            case GenotypeSource::TRUE : part1 = "True Genotypes"; break;
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
            case GenotypeSource::EXPECTED :
                if(isRVS()) part1 = "RVS"; else part1 = "No RVS";
                break;
            case GenotypeSource::CALL : part1 = "Call"; break;
            case GenotypeSource::VCF_CALL : part1 = "VCF"; break;
            case GenotypeSource::TRUE : part1 = "True"; break;
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
