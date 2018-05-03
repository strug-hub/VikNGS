#include "simulation.h"

/*
Simulates a dataset which can be used for an association test.
@param simReq Parameters required to simulate a dataset.
@return TestInput with information for running association test.
*/
std::vector<TestInput> simulate(std::vector<SimulationRequest> simReqs) {

    SimulationRequest sr = simReqs[0];

	printInfo("Setting up simulation parameters.");

    int npop = sr.npop;
    double prevalence = sr.prevalence;

	int ncase_pop = floor(npop * prevalence);

    int nsnp = sr.nsnp;

    double me = sr.me;
    double sde = sr.sde;

    double oddsRatio = sr.oddsRatio;
    double maf = sr.maf;

    VectorXd mafs = simulateMinorAlleleFrequency(nsnp, maf, maf);
    //todo:?
    /* or we can determine the minor allele frequency fixed for each collapsed SNPs(5 SNPs in the current setting)
        njoint = 5
        loopno = nsnp / njoint

        mafco_5 = c(0.040856639, 0.021548447, 0.015053696, 0.005911941, 0.022559459)
        mafco = rep(mafco_5, loopno)
    */

    std::vector<Variant> variantInfo;

    for (int i = 0; i < nsnp; i++)
        variantInfo.push_back(randomVariant());

    printInfo("Simulating population data.");

    MatrixXd Xpop = simulatePopulationX(npop, ncase_pop, oddsRatio, mafs);
    VectorXd Ypop = simulatePopulationY(npop, ncase_pop);

    printInfo("Simulating sample data.");

    std::vector<TestInput> inputs;

    for(SimulationRequest simReq : simReqs){

        int nsamp;
        int ncase = 0;
        int ncont = 0;

        for (SimulationRequestGroup srg : simReq.groups) {
            if (srg.isCase)
                ncase += srg.n;
            else
                ncont += srg.n;
        }

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

        MatrixXd x = sampleX(simReq, Xpop, nsamp, ncase, ncase_pop);
        VectorXd y = sampleY(simReq, Ypop, nsamp, ncase, ncase_pop);

        MatrixXd EG(nsamp, nsnp);
        MatrixXd p(nsnp, 3);

        std::vector<Variant> variants;

        for (int i = 0; i < EG.cols(); i++) {

            Variant info = variantInfo[i];
            Variant variant = generateSeqData(x.col(i), y, g, group, me, sde, info);

            variants.push_back(variant);
            EG.col(i) = variant.expectedGenotype;
            p.row(i) = variant.P;
        }

        //dummy variables
        MatrixXd z;
        std::vector<std::vector<int>> interval;

        TestInput t = buildTestInput(EG, y, z, g, p, readGroup, interval, variants, "binomial");
        inputs.push_back(t);
    }

    return inputs;
}


