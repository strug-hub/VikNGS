#pragma once

#include "../Enums.h"
struct SampleInfo;
struct VariantSet;

#include "../Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::VectorXi;

class TestObject;

//Test.cpp
double runTest(SampleInfo* sampleInfo, VariantSet* variant, Test test, int nboot, bool stopEarly);

//TestRareHelper.cpp
MatrixXd getVarianceMatrix(TestObject& o, Test& test, Family family);
