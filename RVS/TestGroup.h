#pragma once

class TestGroup {
	protected:
		VectorXd Xboot;
		VectorXd Xcenter;

		bool isCentered = false;
		bool isCentered_filterz = false;
		void centerX();
		void centerX_filterz();

		void filterNAN();
		void filterNAN_z();

		VectorXd X_nofilter;
		VectorXd Y_nofilter;
		MatrixXd Z_nofilter;

		VectorXd X;
		VectorXd Y;
		MatrixXd Z;

		VectorXd X_filterz;
		VectorXd Y_filterz;
		MatrixXd Z_filterz;
		VectorXd Xcenter_filterz;

	public:

		VectorXd Ycenter;
		double robustVar;

		TestGroup() {};
		TestGroup(SNP &snp, std::vector<Sample> &sample, Group &group, bool rare);
	
		void bootstrapX();
		void fitModel(VectorXd &beta, std::string distribution);

		inline double score() { return (Ycenter.array() * X.array()).sum(); }

		//may be a lot more work than anticipated....  :(
		//TODO!!!
		//inline double bootstrapScore() { return (Ycenter.array() * X.array()).sum(); }
		
		virtual double variance(bool rvs) = 0;
		inline double bootScore() { return (Ycenter.array() * Xboot.array()).sum(); };
		double bootVariance();

		inline int length() { return X.size(); }
		inline int length_filterz() { return X_filterz.size(); }
		inline VectorXd getX_filterz() { return X_filterz; }
		inline VectorXd getBootstrapX_filterz();

};

class TestHRG : public TestGroup {
	public:
		inline TestHRG(SNP &snp, std::vector<Sample> &sample, Group &group, bool rare)
			: TestGroup(snp, sample, group, rare) {}	
		inline bool isHRG() { return true; }
		double variance(bool rvs);
};

class TestLRG : public TestGroup {
	public:
		inline TestLRG(SNP &snp, std::vector<Sample> &sample, Group &group, bool rare)
			: TestGroup(snp, sample, group, rare) {}
		inline bool isHRG() { return false; }
		double variance(bool rvs);
};