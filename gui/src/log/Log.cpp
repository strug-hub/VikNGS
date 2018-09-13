#include "../src/Log.h"
#include "qlog.h"

#include <QObject>
#include <QString>
#include <QColor>

#include <thread>
#include <chrono>

static int chill = 0;

static QLog* qlog = new QLog();
QLog* getQLog(){
    return qlog;
}

static QString newl = "\n";
static QColor info = QColor::fromRgb(89, 94, 109);
static QColor error = QColor::fromRgb(201, 42, 42);
static QColor warning = QColor::fromRgb(163, 86, 58);

static int sleepMilli = 10;

void printInfo(std::string message) {\
   chill = 0;
   QString s = "[INFO] ";
   s.append(QString::fromStdString(message));
   qlog->out(s, info);

   std::this_thread::sleep_for(std::chrono::milliseconds(sleepMilli));
}

void printWarning(std::string message) {
    chill++;

    if(chill < 3){
        QString s = "[WARNING] ";
        s.append(QString::fromStdString(message));
        qlog->out(s, warning);

        std::this_thread::sleep_for(std::chrono::milliseconds(sleepMilli));
    }
}

void printError(std::string message) {
    chill++;

    if(chill < 3){
        QString s = "[ERROR] ";
        s.append(QString::fromStdString(message));
        qlog->out(s, error);

        std::this_thread::sleep_for(std::chrono::milliseconds(sleepMilli));
    }
}

void throwError(std::string source, std::string message) {
	printError(message);
	throw std::runtime_error(message);
}

void throwError(std::string source, std::string message, std::string valueGiven) {
    printError(message + " Value given: " + valueGiven);
	throw std::runtime_error(message);
}

void printWarning(std::string source, std::string message, std::string valueGiven) {
	printWarning(message + " ||| " + valueGiven);
}

void printWarning(std::string source, std::string message) {
	printWarning(message);
}
