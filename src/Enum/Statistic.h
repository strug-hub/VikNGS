#pragma once

enum class Statistic { NONE, COMMON, CAST, SKAT, CALPHA };

inline bool isRare(Statistic s) {
    return s == Statistic::CAST      ||
            s == Statistic::SKAT     ||
            s == Statistic::CALPHA;}
