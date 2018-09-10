#include "Simulation.h"


/*
Simulates sample information used for an association test.
@param simReq Parameters required to simulate a dataset.
@return SampleInfo with information regarding simulated samples.
*/
SampleInfo simulateSampleInfo(SimulationRequest& simReq) {

    SampleInfo info;
    info.setG(simulateG(simReq));
    info.setY(simulateY(simReq));

    std::map<int, Depth> groupDepth;
    for (SimulationRequestGroup srg : simReq.groups)
        groupDepth[srg.index] = srg.readDepth;

    info.setGroupDepthMap(groupDepth);
    return info;
}

/*
Simulates variant information used for an association test.
@param simReq Parameters required to simulate a dataset.
@return VariantSets with information regarding simulated variants.
*/
std::vector<VariantSet> simulateVariants(SimulationRequest& simReq) {

    std::vector<VariantSet> variants;

    int nsnp = simReq.nsnp;
    int collapseSize = simReq.collapse;
    double oddsRatio = simReq.oddsRatio;
    double minMaf = simReq.mafMin;
    double maxMaf = simReq.mafMax;

    VectorXd maf(collapseSize);

    for (int i = 0; i < nsnp; i+=collapseSize){

        for (int i = 0; i < collapseSize; i++)
            maf[i] = randomDouble(minMaf, maxMaf);

        //note: each matrix is nsnp x nsamp
        MatrixXd X = simulateX(simReq, oddsRatio, maf);

        VariantSet vs;
        variants.push_back(vs);
        VariantSet* setPointer = &variants.back();

        for (int j = 0; j < collapseSize; j++){
            Variant v = randomVariant();
            v.setTrueMaf(maf[j]);

            VectorXd trueX = X.col(j);
            v.setTrueGenotypes(trueX);

            VectorXi readDepths = generateReadDepths(simReq);
            std::vector<std::vector<int>> baseCalls = generateBaseCalls(simReq, trueX, readDepths);
            std::vector<Vector3d> likelihoods = generateLikelihoods(simReq, baseCalls);

            v.setExpectedGenotypes(likelihoods);
            v.setCallGenotypes(likelihoods);

            setPointer->addVariant(v);
        }
    }

    return variants;
}
