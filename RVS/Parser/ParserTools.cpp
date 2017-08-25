#include "InputParser.h"
#include <regex>

/**
Extracts string from a MemoryMapped class

@param charArray MemoryMapped to extract string from.
@param start Index to start extracting string from.
@param end Index to stop extracting string from.
@return Extracted string.
*/
std::string extractString(MemoryMapped &charArray, int start, int end) {
	std::string ret;
	for (; start < end; start++) { ret += charArray[start]; }
	return ret;
}

/**
Removes whitespace from a string

@param str String to remove whitespace from.
@return str without whitespace.
*/
std::string trim(std::string str){
	str.erase(std::remove(str.begin(), str.end(), '\t'), str.end());
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
std::vector<std::string> split(std::string s, char sep) {
	std::vector<std::string> split;
	int start = 0;
	for (int i = 0; i <= s.length(); i++) {
		if (s[i] == sep || i == s.length()) {
			split.push_back(s.substr(start, i - start));
			start = i + 1;
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
	return (std::regex_match(str, std::regex("-?[1234567890]+(\.[1234567890]+)?")));
}
