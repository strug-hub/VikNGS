"""Smoke test for sim-cli: plumbing works, p-values in (0,1)."""
from __future__ import annotations

from tests.stats._sim import run_sim, pvalues_by_source


def test_sim_cli_runs_and_emits_pvalues(sim_binary):
    rows = run_sim(sim_binary, ["--nsnp", "10", "--ncase", "100", "--ncontrol", "100", "--stat", "common"])
    # 10 variants, up to 3 tests per variant (true, expected, call under null)
    assert len(rows) > 0
    assert len({r.variant_idx for r in rows}) == 10
    for r in rows:
        assert 0.0 <= r.pvalue <= 1.0, f"p-value out of range: {r}"
        assert r.test_source in {"true", "expected", "call"}
        assert r.test_stat == "common"


def test_sim_cli_expected_source_present(sim_binary):
    rows = run_sim(sim_binary, ["--nsnp", "5", "--ncase", "50", "--ncontrol", "50", "--stat", "common"])
    pv = pvalues_by_source(rows, "expected")
    assert len(pv) == 5
