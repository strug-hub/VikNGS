"""Power tests: under a fixed alternative, rejection rate should exceed a floor.

Floors were pinned on 2026-04 with N=1000, OR=1.5, ncase=ncontrol=500:
  common/expected: ~0.95
  common/true:     ~0.95

The floor is set well below observed so occasional RNG variation doesn't flake
the suite, while any real power regression (e.g. a score statistic bug that
halves sensitivity) is caught immediately.
"""
from __future__ import annotations

import pytest

from tests.stats._sim import run_sim, pvalues_by_source

NSNP = 1000
NCASE = 500
NCONTROL = 500
OR_ALT = 1.5
ALPHA = 0.05


@pytest.fixture(scope="module")
def alt_common_rows(sim_binary):
    return run_sim(
        sim_binary,
        ["--nsnp", NSNP, "--ncase", NCASE, "--ncontrol", NCONTROL,
         "--stat", "common", "--effect", OR_ALT],
    )


@pytest.mark.parametrize("source,floor", [
    ("true",     0.85),
    ("expected", 0.85),
])
def test_common_power(alt_common_rows, source, floor):
    pv = pvalues_by_source(alt_common_rows, source)
    assert len(pv) == NSNP, f"got {len(pv)} p-values for {source}"
    rate = sum(p < ALPHA for p in pv) / len(pv)
    assert rate > floor, f"{source} power={rate:.3f} below floor {floor}"
