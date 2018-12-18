#pragma once
#include <string>

enum class GenotypeSource { NONE, EXPECTED, TRUE, CALL, VCF_CALL };

inline std::string genotypeToString(GenotypeSource g){
    switch(g) {
        case GenotypeSource::EXPECTED: return "Expected";
        case GenotypeSource::TRUE: return "True";
        case GenotypeSource::CALL: return "Call";
        case GenotypeSource::VCF_CALL: return "VCF Call";
        default: return "_";
    }
}
