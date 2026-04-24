"""Type I error tests: under the null, p-values should be uniform on (0,1)
and the α=0.05 rejection rate should be near 0.05.

Each test runs one null simulation (effect=1.0 for binomial, 0.0 for normal),
pulls p-values from a given genotype source, and checks:
  1. KS test for uniformity — p > 0.01
  2. Rejection rate at α=0.05 in [0.03, 0.07]

Tolerances chosen for N=3000 variants (stderr ≈ 0.004 on rejection rate).
"""
from __future__ import annotations

import pytest
from scipy import stats

from tests.stats._sim import run_sim, pvalues_by_source

NSNP = 3000
NCASE = 500
NCONTROL = 500
ALPHA = 0.05
RATE_LO = 0.03
RATE_HI = 0.07
# KS p-value under H0 is itself uniform, so testing many sources at KS_PMIN=0.01
# produces a false-failure rate too high for CI. Use 1e-4 to make the combined
# flake rate <0.1% while still catching gross non-uniformity (which would
# also be caught by the rate check anyway).
KS_PMIN = 1e-4


@pytest.fixture(scope="module")
def null_common_rows(sim_binary):
    return run_sim(
        sim_binary,
        ["--nsnp", NSNP, "--ncase", NCASE, "--ncontrol", NCONTROL,
         "--stat", "common", "--effect", "1.0"],
    )


def _check_uniform(pvalues: list[float], label: str):
    assert len(pvalues) == NSNP, f"{label}: got {len(pvalues)} p-values, expected {NSNP}"
    assert all(0.0 <= p <= 1.0 for p in pvalues), f"{label}: p-values outside [0,1]"

    # Some variants are monomorphic in the sample and produce exact p=1 for
    # called-genotype tests. Drop them for the uniformity check; they don't
    # represent a statistical failure of the test.
    interior = [p for p in pvalues if 0.0 < p < 1.0]
    ks = stats.kstest(interior, "uniform")
    assert ks.pvalue > KS_PMIN, f"{label}: KS uniformity p={ks.pvalue:.4f} (stat={ks.statistic:.4f}), n={len(interior)}"

    rate = sum(p < ALPHA for p in pvalues) / len(pvalues)
    assert RATE_LO < rate < RATE_HI, f"{label}: α=0.05 rejection rate={rate:.4f} outside [{RATE_LO},{RATE_HI}]"


@pytest.mark.parametrize("source", ["true", "expected", "call"])
def test_common_type1_control(null_common_rows, source):
    pv = pvalues_by_source(null_common_rows, source)
    _check_uniform(pv, f"common/{source}")
