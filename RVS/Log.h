#pragma once

#include <string>
#include <iostream>  

inline void printInfo(std::string message) {
	std::cout << "[INFO] ";
	std::cout << message;
	std::cout << "\n";
}

inline void printWarning(std::string message) {
	std::cout << "[WARNING] ";
	std::cout << message;
	std::cout << "\n";
}

inline void printError(std::string message) {
	std::cout << "[ERROR] ";
	std::cout << message;
	std::cout << "\n";
}
