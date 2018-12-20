#include "Parser.h"
#include <regex>
#include <fstream>

/**
Removes whitespace from a string

@param str String to remove whitespace from.
@return str without whitespace.
*/
std::string trim(std::string str){
    str.erase(std::remove(str.begin(), str.end(), '\t'), str.end());
    str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

/**
Separates a string into a vector, splitting at every postion with
the sep character

@param s String to split.
@param sep Character to split the string at.
@return Split string.
*/
std::vector<std::string> splitString(std::string &s, char sep) {
    std::vector<std::string> split;
    size_t start = 0;
    for (size_t i = 0; i <= s.length(); i++) {
        if (s[i] == sep || i == s.length()) {
            split.emplace_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    return split;
}

/**
Separates a string into a vector, splitting at every postion with
the sep character, up to a given position

@param s String to split.
@param sep Character to split the string at.
@param stop Stop splitting after this index is reached (inclusive).
@return Split string.
*/
std::vector<std::string> splitString(std::string &s, char sep, int stop) {
    std::vector<std::string> split;
    size_t start = 0;
    size_t stp = static_cast<size_t>(stop);
    for (size_t i = 0; i <= s.length(); i++) {
        if (s[i] == sep || i == s.length()) {
            split.emplace_back(s.substr(start, i - start));
            start = i + 1;
            if(split.size() > stp)
                return split;
        }
    }
    return split;
}


/**
Checks if a string is a number

@param str string to check
@return true if str is numeric
*/
inline bool isNumeric(const std::string& str) {
    return (std::regex_match(str, std::regex("-?[1234567890]+(\\.[1234567890]+)?")));
}



//todo: move
/**
Create a vector of indexes used for collapsing by k

@param k number of variants to collapse.
@param n total number of variants.
@return Collapsed indexes.
*/
std::vector<std::vector<int>> collapseEveryK(int k, int n){
  std::vector<std::vector<int>> collapsedVariants;

  int index = 0;

  while(index < n){
    std::vector<int> group;

    for (int i = 0; i < k; i++) {
	group.push_back(index);
	index++;
	if(index >= n)
	  break;
    }

    collapsedVariants.push_back(group);
  }

  return collapsedVariants;
}
