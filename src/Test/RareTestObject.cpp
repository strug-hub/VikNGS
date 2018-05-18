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


void RareTestObject::normalBootstrap(bool covariates) {

    /*std::vector<VectorXd> xboot;
	int i, j, length;
	int c = 0;

    int n = 0;
    for (i = 0; i < size(); i++)
        n += y[i].rows();
    VectorXd Xnew(n);

	for (i = 0; i < size(); i++) {
		length = xcenter[i].rows();
		VectorXd xrand(length);
		for (j = 0; j < length; j++) {
			xrand[j] = xcenter[i][randomInt(0, length - 1)];
			Xnew[c] = xrand[j];
			c++;
		}

		xboot.push_back(xrand);
	}

	this->x = xboot;

    //todo: FIX!!
    //VectorXd beta = getBeta(Xnew, Y, Z, family);
    //ycenter = fitModel(beta, y_, z_, family);

    */
}

void RareTestObject::bootstrap() {

    /*
	std::vector<VectorXd> xboot;
	int i, j, length;

	for (i = 0; i < size(); i++) {
		length = xcenter[i].rows();
		VectorXd xrand(length);
		for (j = 0; j < length; j++)
			xrand[j] = xcenter[i][randomInt(0, length - 1)];

		xboot.push_back(xrand);
	}

	this->x = xboot;

    //todo fix bootstrappign!!!
    double ybar = average(y);
	
	for (int i = 0; i < size(); i++) 
        ycenter[i] = y[i].array() - ybar;

        */
}
