#pragma once
#include <string>
#include "Statistic.h"
#include "GenotypeSource.h"
#include "Variance.h"

struct TestSettings {
private:
    GenotypeSource genotype;
    Statistic statistic;
    Variance variance;
    int nboot = -1;
    int nsamples = -1;
    bool earlyStopping = false;

public:
    TestSettings(GenotypeSource g, Statistic s, Variance v, int nbootstrap=-1, int nsamples=-1) : genotype(g), statistic(s), variance(v){
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
                if(isRVS()) part1 = "vRVS " + part1;
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
                if(isRVS()) part1 = "vRVS"; else part1 = "No vRVS";
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
