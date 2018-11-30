#pragma once
#include "../Math/Math.h"
#include "../vikNGS.h"
#include "TestObject.h"


//Test.cpp
double runTest(SampleInfo* sampleInfo, VariantSet* variant, Test test, int nboot, bool stopEarly);

//TestRareHelper.cpp
MatrixXd getVarianceMatrix(TestObject& o, Test& test, Family family);
