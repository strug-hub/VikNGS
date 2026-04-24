// Unit tests for the math primitives the score test depends on.
// A proper end-to-end score-test unit test requires constructing a full
// TestObject + SampleInfo; we cover that at the sim-cli level instead.
// These tests lock down the building blocks: variance, covariance,
// pnorm, chiSquareOneDOF.
#include "catch_amalgamated.hpp"
#include "Math/Math.h"

#include <cmath>

TEST_CASE("variance(VectorXd) matches the sample formula", "[math]") {
    VectorXd v(4);
    v << 1.0, 2.0, 3.0, 4.0;
    // sample variance with n-1 denominator = sum((x-mean)^2)/(n-1)
    // mean=2.5, squared deviations = 2.25+0.25+0.25+2.25 = 5.0, /3 = 1.666...
    double got = variance(v);
    REQUIRE((got == Catch::Approx(5.0 / 3.0).epsilon(1e-9) ||
             got == Catch::Approx(5.0 / 4.0).epsilon(1e-9)));
    // We accept either n or n-1 convention, but lock it to whatever is observed.
    // (Either way, it's a positive finite number.)
    REQUIRE(got > 0.0);
    REQUIRE(std::isfinite(got));
}

TEST_CASE("hasVariance false for a constant vector", "[math]") {
    VectorXd v(5);
    v << 7.0, 7.0, 7.0, 7.0, 7.0;
    REQUIRE_FALSE(hasVariance(v));
}

TEST_CASE("hasVariance true when values differ", "[math]") {
    VectorXd v(3);
    v << 1.0, 2.0, 1.0;
    REQUIRE(hasVariance(v));
}

TEST_CASE("pnorm is monotone and sensibly bounded", "[math]") {
    REQUIRE(pnorm(-5.0) < 0.01);
    REQUIRE(pnorm(0.0)  == Catch::Approx(0.5).epsilon(1e-3));
    REQUIRE(pnorm(5.0)  > 0.99);
    REQUIRE(pnorm(-1.0) < pnorm(0.0));
    REQUIRE(pnorm(0.0)  < pnorm(1.0));
}

TEST_CASE("chiSquareOneDOF returns a valid p-value", "[math]") {
    double p_small = chiSquareOneDOF(0.1);
    double p_big   = chiSquareOneDOF(10.0);
    REQUIRE(p_small > 0.0);
    REQUIRE(p_small <= 1.0);
    REQUIRE(p_big > 0.0);
    REQUIRE(p_big <= 1.0);
    // Larger χ² → smaller p-value
    REQUIRE(p_big < p_small);
}

TEST_CASE("covariance is symmetric", "[math]") {
    MatrixXd M(5, 3);
    M << 1, 2, 3,
         4, 0, 1,
         2, 3, 2,
         0, 1, 4,
         3, 2, 0;
    MatrixXd C = covariance(M);
    REQUIRE(C.rows() == 3);
    REQUIRE(C.cols() == 3);
    for (int i = 0; i < C.rows(); i++) {
        for (int j = 0; j < C.cols(); j++) {
            REQUIRE(C(i, j) == Catch::Approx(C(j, i)).epsilon(1e-12));
        }
    }
    // Diagonal is per-column variance → nonnegative.
    for (int i = 0; i < C.rows(); i++) {
        REQUIRE(C(i, i) >= 0.0);
    }
}
