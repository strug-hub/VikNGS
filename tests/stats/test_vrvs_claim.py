"""Validates the paper's central claim: under depth-imbalanced null data,
the vRVS (expected-genotype) test maintains α=0.05 calibration while the
naive called-genotype test inflates.

Scenario: case cohort has mean depth 50, error 0.001 (near-perfect sequencing).
Control cohort has mean depth 3, error 0.05 (shallow + noisy). No true effect.

Expected behavior (verified 2026-04):
  expected (vRVS): ≈0.047  (calibrated)
  call:            ≈0.170  (3-4× inflated)
  true:            ≈0.048  (calibrated; true genotypes are depth-invariant)

If 'expected' inflates above ~0.07 here, the vRVS correction is broken.
If 'call' fails to inflate, the simulation isn't actually creating depth-driven
bias and the other vRVS tests are untrustworthy.
"""
from __future__ import annotations

import pytest

from tests.stats._sim import run_sim, pvalues_by_source

NSNP = 3000
NCASE = 500
NCONTROL = 500
ALPHA = 0.05


@pytest.fixture(scope="module")
def imbalanced_null_rows(sim_binary):
    return run_sim(
        sim_binary,
        ["--nsnp", NSNP, "--ncase", NCASE, "--ncontrol", NCONTROL,
         "--stat", "common", "--effect", "1.0",
         "--case-depth", "50", "--case-depth-sd", "5",
         "--case-high", "1", "--case-error", "0.001",
         "--control-depth", "3", "--control-depth-sd", "1",
         "--control-high", "0", "--control-error", "0.05"],
    )


def _rate(rows, source):
    pv = pvalues_by_source(rows, source)
    return sum(p < ALPHA for p in pv) / len(pv), len(pv)


def test_vrvs_maintains_calibration_under_depth_imbalance(imbalanced_null_rows):
    rate, n = _rate(imbalanced_null_rows, "expected")
    assert n == NSNP
    assert rate < 0.07, (
        f"vRVS (expected) rejection rate={rate:.4f} exceeds 0.07 under depth imbalance "
        "— vRVS variance correction may be broken"
    )


def test_called_genotype_inflates_under_depth_imbalance(imbalanced_null_rows):
    rate, n = _rate(imbalanced_null_rows, "call")
    assert n == NSNP
    assert rate > 0.10, (
        f"Called-genotype rejection rate={rate:.4f} should inflate above 0.10 under this "
        "extreme depth imbalance. If it doesn't, the simulation isn't producing "
        "depth-driven bias and other vRVS claims become untestable."
    )


def test_vrvs_beats_call_under_depth_imbalance(imbalanced_null_rows):
    rate_expected, _ = _rate(imbalanced_null_rows, "expected")
    rate_call, _ = _rate(imbalanced_null_rows, "call")
    assert rate_expected < rate_call - 0.05, (
        f"vRVS ({rate_expected:.4f}) should be substantially better calibrated "
        f"than call ({rate_call:.4f}) under depth imbalance."
    )
