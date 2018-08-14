#pragma once

#include <string>
#include <map>
#include <vector>
#include <algorithm>

struct Interval{
    std::string id;
    int index;
    int exon;

    std::string chr;
    int start;
    int end;

    inline bool isIn (std::string chrom, int pos){
        return chr == chrom && pos >= start && pos <= end;
    }
    inline bool isIn (int pos){
        return pos >= start && pos <= end;
    }

    inline bool isSmaller (int pos){ return pos < start; }
    inline bool isBigger (int pos){ return pos > end; }

    inline bool operator < (Interval& inv) {
        if (this->chr == inv.chr)
            return this->start < inv.start;
        else
            return this->chr < inv.chr;
    }
};

inline bool intervalCompare(Interval lhs, Interval rhs) { return lhs < rhs; }

struct IntervalSet{

private:
    std::map<std::string, std::vector<Interval>> chrMap;
public:
    void addInterval(Interval &inv){ chrMap[inv.chr].push_back(inv); }
    void sort(){
        for ( const auto &keyValue : chrMap )
            std::sort(keyValue.second.begin(), keyValue.second.end(), intervalCompare);
    }

    std::vector<Interval> * get(std::string key) { return &chrMap[key]; }
    bool isEmpty(){ return chrMap.empty(); }
};
