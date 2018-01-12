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

inline void throwError(std::string source, std::string message) {
	printError(message);

	throw std::runtime_error(message);
}

inline void throwError(std::string source, std::string message, std::string valueGiven) {
	printError(message);

	throw std::runtime_error(message);
}

inline void printWarning(std::string source, std::string message, std::string valueGiven) {
	printWarning(message + " ||| " + valueGiven);
}

inline void printWarning(std::string source, std::string message) {
	printWarning(message);
}