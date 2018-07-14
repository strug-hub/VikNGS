#include "simulation.h"
#include "../global.h"

/*
Simulates a dataset which can be used for an association test.
@param simReq Parameters required to simulate a dataset.
@return TestInput with information for running association test.
*/
std::vector<TestInput> simulate(SimulationRequest& simReq) {

	printInfo("Setting up simulation parameters.");

    int nsnp = simReq.nsnp;
    double oddsRatio = simReq.oddsRatio;
    double mafMin = simReq.mafMin;
    double mafMax = simReq.mafMax;

    printInfo("Simulating data.");

    VectorXd mafs = simulateMinorAlleleFrequency(nsnp, mafMin, mafMax);

    std::vector<Variant> variants;
    for (int i = 0; i < nsnp; i++)
        variants.push_back(randomVariant());

    std::map<int, SimulationRequestGroup> group;
    std::map<int, int> readGroup;

    for (SimulationRequestGroup srg : simReq.groups) {
        group[srg.index] = srg;
        readGroup[srg.index] = srg.isHrg;
    }

    std::vector<VectorXd> G = simulateG(simReq);
    std::vector<VectorXd> Y = simulateY(simReq);
    //note: each matrix is nsnp x nsamp
    std::vector<MatrixXd> X_ = simulateX(simReq, oddsRatio, mafs, simReq.collapse);

    std::vector<TestInput> inputs;

    VectorXd g;
    MatrixXd x;
    VectorXd y;
    //empty variable
    MatrixXd z;

    std::vector<std::vector<int>> collapse;
    if(simReq.isRare()){
        std::vector<int> c;
        for(int i = 0; i < nsnp; i++){
            if(i%simReq.collapse == 0 && i > 0){
                collapse.push_back(c);
                c.clear();
            }
            c.push_back(i);
        }

        collapse.push_back(c);
    }

    for(int i = 0; i < simReq.steps; i++){

        printInfo("Simulating step " + std::to_string(i+1));
    
        for (int j = 0; j < X_[i].cols(); j++) {
            
            if(STOP_RUNNING_THREAD)
                throw std::runtime_error("terminate thread");

            int groupID = G[i][j];
            VectorXd x_i = X_[i].col(j);
            std::vector<GenotypeLikelihood> likelihoods = generateSeqData(x_i, group.at(groupID));
            for(int k = 0; k < nsnp; k++)
                variants[k].likelihood.emplace_back(likelihoods[k]);
        }

        if(i == 0){
            x = X_[0].transpose();
            y = Y[0];
            g = G[0];

        }
        else{
            MatrixXd newX(x.rows() + X_[i].cols(), x.cols());
            newX << x, X_[i].transpose();
            x = newX;

            VectorXd newY(y.rows() + Y[i].rows());
            newY << y, Y[i];
            y = newY;

            VectorXd newG(g.rows() + G[i].rows());
            newG << g, G[i];
            g = newG;
        }

        for(int k = 0; k < nsnp; k++){
            variants[k].trueGenotype = x.col(k);
            variants[k].calculateGenotypeFrequency();
            variants[k].expectedGenotype = calculateExpectedGenotypes(variants[k].likelihood, variants[k].P);
            variants[k].genotypeCalls = calculateGenotypeCalls(variants[k].likelihood);
        }

        TestInput t = buildTestInput(y, z, g, readGroup, variants, "binomial");

        if(simReq.isRare())
            t.addCollapse(collapse);

        inputs.push_back(t);
    }

    return inputs;
}


