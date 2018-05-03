#pragma once

#include <string>
#include <iostream>  

void printInfo(std::string message);
void printWarning(std::string message);
void printError(std::string message);
void throwError(std::string source, std::string message);
void throwError(std::string source, std::string message, std::string valueGiven);
void printWarning(std::string source, std::string message, std::string valueGiven);
void printWarning(std::string source, std::string message);

class variant_exception: public std::exception{
    std::string message;

public:
    variant_exception(std::string msg){
        this->message= msg;
    }

  virtual const char* what() const throw()
  {
        std::string what= "Error in single variant p-value computation: " + this->message;
        return what.c_str();
  }
};
