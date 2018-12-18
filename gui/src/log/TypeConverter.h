#pragma once

#include <string>
#include <iostream>
#include <limits>
#include "../src/Log.h"

static const int infinity = std::numeric_limits<int>::max();
static const int neginfinity = std::numeric_limits<int>::min();

inline double toDouble(std::string valueName, std::string value){
    try{
        return std::stod(value);
    }catch(...){
        throwError("TypeConverter", valueName + " is expected to be a decimal value.", value );
    }
    return -1;
}

inline int toInt(std::string valueName, std::string value){
    try{
        return std::stoi(value);
    }catch(...){
        throwError("TypeConverter", valueName + " is expected to be an integer value.", value );
    }
    return -1;
}

inline std::string toString(double value){
    std::string v = std::to_string(value);
    v = v.erase ( v.find_last_not_of('0') + 1, std::string::npos );
    if(v.back() == '.')
        v = v + "0";
    return v;
}


inline std::string toString(int value){
    return std::to_string(value);
}
