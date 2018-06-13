#include "RareTestObject.h"

double RareTestObject::getYmHigh(std::vector<VectorXd> & ycenter) {
    double ym = 0;
    for (int i = 0; i < size(); i++)
        if (isHRG(i))
            ym += ycenter[i].array().pow(2).sum();

    return ym;
}

double RareTestObject::getYmLow(std::vector<VectorXd> & ycenter) {
    double ym = 0;
    for (int i = 0; i < size(); i++)
        if (!isHRG(i))
            ym += ycenter[i].array().pow(2).sum();

    return ym;
}

double RareTestObject::getYm(std::vector<VectorXd> & ycenter) {
    double ym = 0;
    for (int i = 0; i < size(); i++)
        ym += ycenter[i].array().pow(2).sum();

    return ym;
}

void RareTestObject::bootstrap() {

    x = groupwiseShuffleWithReplacement(xcenter);
}

void RareTestObject::regularBootstrap() {

    x = shuffleWithoutReplacement(xcenter);
}
