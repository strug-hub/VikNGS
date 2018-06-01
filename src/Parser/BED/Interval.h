#pragma once
#include "../../Variant.h"

struct Interval{
    std::string id;
    int index;
    int exon;

    std::string chr;
    int start;
    int end;
};

struct BEDInterval {
    bool valid = false;
    std::string id;

    int txStart = -1;
    int txEnd = -1;
    int cdsStart = -1;
    int cdsEnd = -1;
    std::vector<int> exonStart;
    std::vector<int> exonEnd;

    std::vector<int> indexes;
    std::string chr;

    void print() {
        std::cout << "Interval " + id + " (" + std::to_string(indexes.size()) + " variants)";
        std::cout << "\n\tvalid = ";
        std::cout << valid;
        std::cout << "\n\tTranscription: ";
        std::cout << std::to_string(txStart) + " to " + std::to_string(txEnd);

        if (cdsStart > -1) {
            std::cout << "\n\tCoding: ";
            std::cout << std::to_string(cdsStart) + " to " + std::to_string(cdsEnd);
        }
        for (int i = 0; i < exonStart.size(); i++) {
            std::cout << "\n\tExon " + std::to_string(i) + ": ";
            std::cout << std::to_string(exonStart[i]) + " to " + std::to_string(exonEnd[i]);
        }
        std::cout << "\n";
    }
};


/*

inline void txAddIfIn(Variant variant, int index) {
    if (variant.pos >= txStart && variant.pos <= txEnd && chr.compare(variant.chr) == 0)
        indexes.push_back(index);
}

inline void exonAddIfIn(Variant variant, int index) {
    if (variant.pos >= txStart && variant.pos <= txEnd && chr.compare(variant.chr) == 0) {
        for (int i = 0; i < exonStart.size(); i++) {
            if (variant.pos >= exonStart[i] && variant.pos <= exonEnd[i]) {
                indexes.push_back(index);
                return;
            }
        }
    }
}

*/
