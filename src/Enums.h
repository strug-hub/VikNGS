#pragma once
#include <string>

enum class Statistic { NONE, COMMON, SKAT, CALPHA };
enum class Genotype { NONE, EXPECTED, TRUE, CALL, VCF_CALL };
enum class Family { NORMAL, BINOMIAL, NONE };
enum class Depth { HIGH, LOW };

enum class Filter { VALID, INVALID, IGNORE, NOT_SNP, NO_PASS, MISSING_DATA, NO_VARIATION, MAF };
enum class CollapseType { COLLAPSE_K, COLLAPSE_GENE, COLLAPSE_EXON, NONE };

struct Test {
private:
    Genotype genotype;
    Statistic statistic;
    int rvs = 1;

public:
    Test(Genotype g, Statistic s) : genotype(g), statistic(s){ rvs=1; }
    Test(Genotype g, Statistic s, bool rvs) : genotype(g), statistic(s){ if(rvs) this->rvs=1; else this->rvs = 0; }

    inline bool isRareTest(){ return statistic == Statistic::SKAT || statistic == Statistic::CALPHA; }
    inline bool isCommonTest(){ return statistic == Statistic::COMMON; }
    inline bool useRVS(){ return rvs == 1; }

    inline bool needVCFCalls(){ return genotype == Genotype::VCF_CALL; }
    inline bool needGenotypeCalls(){ return genotype == Genotype::CALL; }
    inline bool needExpectedGenotypes(){ return genotype == Genotype::EXPECTED; }

    inline std::string toString(){
        std::string part1, part2;
        switch(genotype)
        {
            case Genotype::EXPECTED :
                part1 = "Genotype Likelihoods";
                if(useRVS()) part1 = "RVS " + part1;
                break;
            case Genotype::CALL : part1 = "Genotype Calls"; break;
            case Genotype::VCF_CALL : part1 = "VCF Calls"; break;
            case Genotype::TRUE : part1 = "True Genotypes"; break;
            default: part1 = "???"; break;
        }
        switch(statistic)
        {
            case Statistic::COMMON : part2 = "common"; break;
            case Statistic::CALPHA : part2 = "Calpha"; break;
            case Statistic::SKAT : part2 = "SKAT"; break;
            default: return "???";
        }

        return part1 + " - " + part2;
    }

    inline std::string toShortString(){
        std::string part1, part2;
        switch(genotype)
        {
            case Genotype::EXPECTED :
                if(useRVS()) part1 = "RVS"; else part1 = "No RVS";
                break;
            case Genotype::CALL : part1 = "Call"; break;
            case Genotype::VCF_CALL : part1 = "VCF"; break;
            case Genotype::TRUE : part1 = "True"; break;
            default: part1 = "???"; break;
        }
        switch(statistic)
        {
            case Statistic::COMMON : part2 = "common"; break;
            case Statistic::CALPHA : part2 = "Calp"; break;
            case Statistic::SKAT : part2 = "SKAT"; break;
            default: return "???";
        }

        return part1 + " " + part2;
    }
};
