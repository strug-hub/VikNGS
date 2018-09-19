#pragma once
#include <string>

void printInfo(std::string message);
void printWarning(std::string message);
void printError(std::string message);
void throwError(std::string source, std::string message);
void throwError(std::string source, std::string message, std::string valueGiven);
void printWarning(std::string source, std::string message, std::string valueGiven);
void printWarning(std::string source, std::string message);

void printBedIDWarning(std::string chr);
