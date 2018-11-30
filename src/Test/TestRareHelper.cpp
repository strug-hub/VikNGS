#include "Test.h"
#include "TestObject.h"
#include "ScoreTestFunctions.h"

MatrixXd getVarianceMatrix(TestObject& o, Test& test, Family family){
    if(test.isExpectedGenotypes()){

        if(family == Family::BINOMIAL)
            return getRobustVarianceBinomial(*o.getYcenter(), *o.getX(), *o.getG(),
                                       *o.getDepths(), o.robustVarVector(), test.isRVS());
        if(family == Family::NORMAL)
            return getRobustVarianceNormal(*o.getYcenter(), *o.getX(), *o.getG(),
                                     *o.getDepths(), o.robustVarVector(), test.isRVS());
    }
    else
            return getRegularVariance(*o.getYcenter(), *o.getX(), *o.getZ(), *o.getMU(), family);


    throwError("TEST", "Don't know how to compute variance during test, this error should not happen.");
}

