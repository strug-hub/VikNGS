#include "../RVS.h"
#include "../VectorHelper.h"

class RareTestObject {
private:
	VectorXd Y;
	MatrixXd Z;
	bool hasCov;

	//filter NAN from Y and Z only
	//this x stays the same after every bootstrap!
	std::vector<VectorXd> x;
	std::vector<VectorXd> y;
	std::vector<MatrixXd> z;

	std::vector<VectorXd> x_;
	std::vector<VectorXd> y_;
	std::vector<MatrixXd> z_;

	std::vector<int> readDepth;
	std::vector<VectorXd> xcenter;
	std::vector<VectorXd> ycenter;
	double robustVar;
	double nhrd;
	double nlrd;
	double nhrd_;
	double nlrd_;

	void filterNAN();
	void filterNAN_yz();
	void filterNAN_xy();
	void filterNAN_YZ();

	void countRD();

public:
	RareTestObject(VectorXd &X, VectorXd &Y, MatrixXd &Z,
		std::vector<VectorXd> &x, std::vector<VectorXd> &y, std::vector<MatrixXd> &z,
		std::vector<int> &readDepth, VectorXd &P) {
		hasCov = true;
		this->Y = Y;
		this->Z = Z;
		filterNAN_YZ();

		this->x = x;
		this->y = y;
		this->z = z;
		filterNAN_yz();

		this->readDepth = readDepth;
		this->robustVar = calcRobustVar(P);

		filterNAN();
		countRD();

		VectorXd beta = getBeta(X, Y, Z);
		ycenter = fitModel(beta, y_, z_, "norm");

		for (int i = 0; i < size(); i++)
			this->xcenter.push_back(this->x[i].array() - x_[i].mean());
	}

	RareTestObject(std::vector<VectorXd> &x, std::vector<VectorXd> &y,
		std::vector<int> &readDepth, VectorXd &P) {
		hasCov = false;
		this->x = x;
		this->y = y;

		this->readDepth = readDepth;
		this->robustVar = calcRobustVar(P);

		filterNAN_xy();
		countRD();

		double ybar = average(y_);

		for (int i = 0; i < size(); i++) {
			xcenter.push_back(this->x[i].array() - x_[i].mean());
			ycenter.push_back(y_[i].array() - ybar);
		}

	}

	inline double getRobustVar() { return robustVar; }
	inline double getYm(int depth) {
		double ym = 0;
		for (int i = 0; i < size(); i++)
			if (readDepth[i] == depth)
				ym += ycenter[i].array().pow(2).sum();

		if (depth == 1)  ym = ym / nhrd_ * nhrd;
		else  ym = ym / nlrd_ * nlrd;

		return sqrt(ym);
	}

	inline double getScore() {
		double score = 0;
		for (int i = 0; i < size(); i++) 
			score += (ycenter[i].array() * x_[i].array()).sum();
		return score;
	}

	inline bool isHRG(int i) { return readDepth[i] == 1; }
	inline VectorXd getX(int i) { return x[i]; }
	inline int size() { return x.size(); }

	//bootstrapping ------------------------------

	void covBootstrap();
	void noCovBootstrap();

	inline void bootstrap() {
		if (hasCov)
			covBootstrap();
		else
			noCovBootstrap();
	}

	//bootstrapping ------------------------------

};
