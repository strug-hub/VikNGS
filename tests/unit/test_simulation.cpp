// Unit tests for the simulation helpers (src/gui/src/simulation/Simulation*.cpp).
#include "catch_amalgamated.hpp"
#include "gui/src/simulation/Simulation.h"

#include <algorithm>

static SimulationRequest makeBinomialRequest(int nsnp, double mafMin, double mafMax) {
    SimulationRequest req;
    req.nsnp = nsnp;
    req.effectSize = 1.0;   // null
    req.mafMin = mafMin;
    req.mafMax = mafMax;
    req.steps = 1;
    req.collapse = 1;
    req.nboot = 1;
    req.useBootstrap = false;
    req.stopEarly = false;
    req.nthreads = 1;
    req.family = Family::BINOMIAL;
    req.testStatistic = Statistic::COMMON;

    SimulationRequestGroup caseG;
    caseG.index = 0; caseG.n = 100; caseG.family = Family::BINOMIAL;
    caseG.isCase = true; caseG.meanDepth = 20; caseG.sdDepth = 2;
    caseG.errorRate = 0.01; caseG.readDepth = Depth::HIGH;
    req.groups.push_back(caseG);

    SimulationRequestGroup ctrlG;
    ctrlG.index = 1; ctrlG.n = 100; ctrlG.family = Family::BINOMIAL;
    ctrlG.isCase = false; ctrlG.meanDepth = 20; ctrlG.sdDepth = 2;
    ctrlG.errorRate = 0.01; ctrlG.readDepth = Depth::HIGH;
    req.groups.push_back(ctrlG);

    return req;
}

TEST_CASE("generateMafs returns values in [mafMin, mafMax]", "[simulation]") {
    auto req = makeBinomialRequest(500, 0.05, 0.30);
    VectorXd mafs = generateMafs(req);
    REQUIRE(mafs.size() == 500);
    double mn = mafs.minCoeff();
    double mx = mafs.maxCoeff();
    REQUIRE(mn >= 0.05);
    REQUIRE(mx <= 0.30);
}

TEST_CASE("simulateYCaseControl returns 0/1 phenotypes of correct length", "[simulation]") {
    auto req = makeBinomialRequest(10, 0.05, 0.5);
    MatrixXd Y = simulateYCaseControl(req);
    REQUIRE(Y.rows() == 200);  // 100 case + 100 control
    for (int i = 0; i < Y.rows(); i++) {
        double v = Y(i, 0);
        REQUIRE((v == 0.0 || v == 1.0));
    }
    // First 100 are cases (= 1), next 100 are controls (= 0) per group layout.
    REQUIRE(Y.block(0, 0, 100, 1).sum() == Catch::Approx(100.0));
    REQUIRE(Y.block(100, 0, 100, 1).sum() == Catch::Approx(0.0));
}

TEST_CASE("SimulationRequest::underNull respects family semantics", "[simulation]") {
    auto req = makeBinomialRequest(10, 0.05, 0.5);
    req.effectSize = 1.0;
    REQUIRE(req.underNull());
    req.effectSize = 1.5;
    REQUIRE_FALSE(req.underNull());

    req.family = Family::NORMAL;
    req.effectSize = 0.0;
    REQUIRE(req.underNull());
    req.effectSize = 0.1;
    REQUIRE_FALSE(req.underNull());
}
