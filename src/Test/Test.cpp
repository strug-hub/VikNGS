#include "Test.h"
#include "../Eigen/src/MatrixFunctions/MatrixSquareRoot.h"
#include "../Math/CompQuadForm.h"
#include "ScoreTestFunctions.h"
#include "TestObject.h"

/*
Calculates test statistic.

@param TestObject Test object containing data.
@param test Indicates which test to use.
@param family Statistical distribution family.

@return test statistic
*/
double calculateTestStatistic(TestObject& o, Test& test, Family family, bool print) {

    VectorXd score = getScoreVector(*o.getYcenter(), *o.getX());
    MatrixXd variance = getVarianceMatrix(o, test, family);

    Statistic s = test.getStatistic();

    if(s == Statistic::COMMON || s == Statistic::CAST){

        double testStat = std::pow(score.sum(), 2) / variance.sum();
        return chiSquareOneDOF(testStat);

    }

    if(s == Statistic::SKAT || s == Statistic::CALPHA){

        auto g = variance.eigenvalues();
        VectorXd f = variance.eigenvalues().real();

        //todo: temporary solution to eigenvalue issue
        bool useCalpha = f.sum() < 1e-4;
        for(int e = 0; e < f.rows(); e++)
            if(f[e] < 0)
                useCalpha = true;

        if (s == Statistic::CALPHA)
            useCalpha = true;

        if(useCalpha){
            VectorXd e = variance.eigenvalues().real();
            std::vector<double> eigenvals(e.data(), e.data() + e.size());
            CQF pval;
            return pval.qfc(eigenvals, score.array().pow(2).sum(), score.rows());
        }

        //skat-Z
        VectorXd A = o.mafWeightVector();

        double quad = 0;
        for(int i = 0; i < score.rows(); i++)
            quad += score[i]*A[i]*score[i];

        MatrixXd sigma = variance.sqrt() * A.asDiagonal() * variance.sqrt();
        VectorXd e = sigma.eigenvalues().real();
        std::vector<double> eigenvalues(e.data(), e.data() + e.size());
        CQF pval;
        return pval.qfc(eigenvalues, quad, score.rows());
    }

    throwError("Test", "Unsure which score test to use. This should not happen.");
    return NAN;
}

double bootstrapTest(double testStatistic, TestObject& o, Test bootTest, Family bootFam, int nboot, bool stopEarly){

    bootTest.setRVSFalse();

    double bootCount = 0;
    double tcount = 0;
    double tsamp;

    for (int h = 0; h < nboot; h++) {

        if(STOP_RUNNING_THREAD)
            return NAN;

        o.bootstrap(bootTest, bootFam);

        tsamp = calculateTestStatistic(o, bootTest, bootFam, false);

        if (std::abs(tsamp) <= std::abs(testStatistic))
            tcount++;

        bootCount++;

        if (stopEarly && bootCount > 10) {
            double pstar = 5 / ((bootCount + 0.0737224)*(1 + 0.0737224));
            if (tcount / (1.0*h) > pstar)
                break;
        }
    }

    return (tcount + 1) / (bootCount + 1) ;
}

double runTest(SampleInfo* sampleInfo, VariantSet* variant, Test test, int nboot, bool stopEarly){

    if(STOP_RUNNING_THREAD)
        return NAN;
    if(variant->validSize() < 1)
        return NAN;

    MatrixXd X = variant->getX(test.getGenotype());
    VectorXd Y = sampleInfo->getY();
    MatrixXd Z = sampleInfo->getZ();
    if(Z.cols()<1)
        Z = MatrixXd::Constant(Y.rows(), 1, 1);
    MatrixXd P = variant->getP(test.getGenotype());
    VectorXi G = sampleInfo->getG();
    std::map<int, Depth> map = sampleInfo->getGroupDepthMap();

    if(test.getSampleSize() > 0){
        int size = test.getSampleSize();
        MatrixXd x = X.block(0, 0, size, X.cols()); X = x;
        VectorXd y = Y.block(0, 0, size, Y.cols()); Y = y;
        if(Z.rows() > 0){
            MatrixXd z = Z.block(0, 0, size, Z.cols()); Z = z;
        }
        VectorXi g = G.block(0, 0, size, G.cols()); G = g;
    }

    Group group(G, map);
    Genotype geno(X, P, test.getGenotype());
    Phenotype pheno(Y, Z, sampleInfo->getFamily());

    //if( test.isRVS() && readDepthFlip(Y, group, sampleInfo->getFamily()))
   //     test.setRegularVariance();

    TestObject o(geno, pheno, group, test.isRareTest());

    double testStatistic = calculateTestStatistic(o, test, sampleInfo->getFamily(), true);

    if(nboot > 1)
        return bootstrapTest(testStatistic, o, test, sampleInfo->getFamily(), nboot, stopEarly);

    return testStatistic;
}
