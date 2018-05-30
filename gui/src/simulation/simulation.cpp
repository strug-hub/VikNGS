#include "simulation.h"

/*
Simulates a dataset which can be used for an association test.
@param simReq Parameters required to simulate a dataset.
@return TestInput with information for running association test.
*/
std::vector<TestInput> simulate(std::vector<SimulationRequest> simReqs) {

    SimulationRequest sr = simReqs[0];

	printInfo("Setting up simulation parameters.");

    int nsnp = sr.nsnp;
    double oddsRatio = sr.oddsRatio;
    double mafMin = sr.mafMin;
    double mafMax = sr.mafMax;

    VectorXd mafs = simulateMinorAlleleFrequency(nsnp, mafMin, mafMax);

    std::vector<Variant> variantInfo;

    for (int i = 0; i < nsnp; i++)
        variantInfo.push_back(randomVariant());

    printInfo("Simulating data.");

    std::vector<TestInput> inputs;

    for(SimulationRequest simReq : simReqs){

        int nsamp;
        int ncase = simReq.ncase();
        int ncont = simReq.ncontrol();

        nsamp = ncase + ncont;

        std::map<int, SimulationRequestGroup> group;
        std::map<int, int> readGroup;
        VectorXd g(nsamp);

        int counter = 0;
        for (SimulationRequestGroup srg : simReq.groups) {
            group[srg.index] = srg;
            readGroup[srg.index] = srg.isHrg;
            for(int i = counter; i < counter + srg.n; i++ )
                g[i] = srg.index;
            counter += srg.n;
        }

        MatrixXd x;

        if(simReq.isRare())
            x = simulateX(simReq, oddsRatio, mafs, simReq.collapse);

        else
            x = simulateX(simReq, oddsRatio, mafs);

        VectorXd y = simulateY(simReq);

        MatrixXd EG(nsamp, nsnp);
        MatrixXd p(nsnp, 3);

        std::vector<Variant> variants;

        for (int i = 0; i < EG.cols(); i++) {

            Variant info = variantInfo[i];
            Variant variant = generateSeqData(x.col(i), y, g, group, info);

            variants.push_back(variant);
            EG.col(i) = variant.expectedGenotype;
            p.row(i) = variant.P;
        }

        //empty variables
        MatrixXd z;
        std::vector<std::vector<int>> interval;

        if(simReq.isRare()){

            std::vector<int> inv;
            for(int i = 0; i < nsnp; i++){
                if(i%simReq.collapse == 0 && i > 0){
                    interval.push_back(inv);
                    inv.clear();
                }
                inv.push_back(i);
            }
        }

        TestInput t = buildTestInput(EG, y, z, g, p, readGroup, interval, variants, "binomial");
        inputs.push_back(t);
    }

    return inputs;
}


