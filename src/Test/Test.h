#pragma once
#include "../Math/EigenStructures.h"

enum class Family;
struct SampleInfo;
struct VariantSet;

class TestSettings;
class TestObject;

//Test.cpp
double runTest(SampleInfo* sampleInfo, VariantSet* variant, TestSettings test, int nboot, bool stopEarly);

//TestRareHelper.cpp
MatrixXd getVarianceMatrix(TestObject& o, TestSettings& test, Family family);
