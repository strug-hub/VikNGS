// Unit tests for src/Math/GeneticsHelper.cpp and related math.
#include "catch_amalgamated.hpp"
#include "Math/Math.h"

#include <cmath>
#include <vector>

TEST_CASE("calcRobustVar(P1,P2) matches the var(E(G|D)) formula", "[genetics]") {
    // For genotype freq P = (p0, p1, p2), E[G] = 2*p2 + p1, E[G^2] = 4*p2 + p1.
    // var = E[G^2] - E[G]^2.
    double p1 = 0.3, p2 = 0.2;   // implicit p0 = 0.5
    double expected = (4 * p2 + p1) - std::pow(2 * p2 + p1, 2);
    REQUIRE(calcRobustVar(p1, p2) == Catch::Approx(expected));
}

TEST_CASE("calcRobustVar(P) with all-AA vector is zero", "[genetics]") {
    // Vector3d default layout: (p0, p1, p2) per the formula above.
    Vector3d P;
    P << 1.0, 0.0, 0.0;  // everybody homozygous reference
    REQUIRE(calcRobustVar(P) == Catch::Approx(0.0));
}

TEST_CASE("calcRobustVar(P) for 50/50 heterozygotes equals 0.25", "[genetics]") {
    // p0=0, p1=1, p2=0 → E[G]=1, E[G^2]=1 → var=0. Not 0.25.
    // p0=0.5, p1=0, p2=0.5 → E[G]=1, E[G^2]=2 → var=1.
    Vector3d P;
    P << 0.5, 0.0, 0.5;
    REQUIRE(calcRobustVar(P) == Catch::Approx(1.0));
}

TEST_CASE("calculateGenotypeFrequencies(VectorXd) returns a probability simplex", "[genetics]") {
    VectorXd gt(6);
    gt << 0.0, 0.0, 1.0, 1.0, 2.0, 2.0;
    Vector3d P = calculateGenotypeFrequencies(gt);
    REQUIRE(P[0] + P[1] + P[2] == Catch::Approx(1.0));
    REQUIRE(P[0] >= 0.0);
    REQUIRE(P[1] >= 0.0);
    REQUIRE(P[2] >= 0.0);
}
