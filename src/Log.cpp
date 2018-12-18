#include "Log.h"
#include <iostream>

std::string newline = "\n";

void printInfo(std::string message) {
	std::cout << "[INFO] " << message << newline;
}

void printWarning(std::string message) {
	std::cout << "[WARNING] " << message << newline;
}

void printError(std::string message) {
	std::cout << "[ERROR] " << message << newline;
}

void throwError(std::string source, std::string message) {
	printError(message);
	throw std::runtime_error(message);
}

void throwError(std::string source, std::string message, std::string valueGiven) {
	printError(message);
	throw std::runtime_error(message);
}

void printWarning(std::string source, std::string message, std::string valueGiven) {
	printWarning(message + " ||| " + valueGiven);
}

void printWarning(std::string source, std::string message) {
	printWarning(message);
}

bool BED_WARNING_PRINTED=false;
void printBedIDWarning(std::string chr){
    if(!BED_WARNING_PRINTED)
        printWarning("There was a chromosome name in the VCF file (" + chr + ") which was not found in the BED file provided. Names should match exactly.");
    BED_WARNING_PRINTED = true;
}

