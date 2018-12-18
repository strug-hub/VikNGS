#pragma once
#include <string>

enum class Filter { NONE, VALID, INVALID, IGNORE, NOT_SNP, NO_PASS, MISSING_DATA, NO_VARIATION, MAF };

inline std::string filterToString(Filter f){
    switch(f) {
        case Filter::VALID: return "Valid";
        case Filter::INVALID: return "Invalid information";
        case Filter::IGNORE: return "Ignored";
        case Filter::NOT_SNP: return "Not SNP";
        case Filter::NO_PASS: return "PASS fail";
        case Filter::MISSING_DATA: return "Missing";
        case Filter::NO_VARIATION: return "No variation";
        case Filter::MAF: return "MAF";
        default: return "?";
    }
}
