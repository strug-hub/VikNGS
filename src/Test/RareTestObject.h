#pragma once
#include "../Math/MathHelper.h"

class RareTestObject {
private:

    //xcenter gets randomized for bootstrap
    std::vector<VectorXd> xcenter;
    std::vector<VectorXd> x;

	std::vector<int> readDepth;
	double robustVar;

public:
    RareTestObject(std::vector<VectorXd> &x, std::vector<int> &readDepth, VectorXd &p) {

		this->x = x;
		this->readDepth = readDepth;
        this->robustVar = calcRobustVar(p);

		for (int i = 0; i < size(); i++)
            this->xcenter.push_back(x[i].array() - x[i].mean());
	}

	inline double getRobustVar() { return robustVar; }
    double getYmHigh(std::vector<VectorXd> & ycenter);
    double getYmLow(std::vector<VectorXd> & ycenter);
    double getYm(std::vector<VectorXd> & ycenter);

    inline double getScore(std::vector<VectorXd> & ycenter) {
		double score = 0;
        for (int i = 0; i < size(); i++)
            score += (ycenter[i].array() * x[i].array()).sum();

        return score;
	}

	inline bool isHRG(int i) { return readDepth[i] == 1; }
	inline VectorXd getX(int i) { return x[i]; }
	inline int size() { return x.size(); }
    inline int xSize(int i) { return x[i].size(); }

    void bootstrap();

};
