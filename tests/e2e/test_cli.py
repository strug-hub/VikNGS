"""End-to-end CLI regression tests: run each case and diff against committed goldens."""
from __future__ import annotations

from pathlib import Path

import pytest

from tests.e2e._diff import diff_pvalues, format_errors
from tests.e2e._harness import run_cli
from tests.e2e.cases import CASES

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"


@pytest.mark.parametrize("case_id,flags", CASES, ids=[c[0] for c in CASES])
def test_cli_matches_golden(case_id, flags, cli_binary, example_vcf, example_info):
    golden_pv = (GOLDEN_DIR / f"{case_id}.pvalues.txt").read_text()
    golden_filt = (GOLDEN_DIR / f"{case_id}.filtered.txt").read_text()

    r = run_cli(cli_binary, example_vcf, example_info, extra_args=flags)
    assert r.returncode == 0, f"CLI failed (rc={r.returncode}): {r.stderr[:500]}"

    errors = diff_pvalues(r.pvalues_text, golden_pv)
    assert not errors, "pvalues drift:\n" + format_errors(errors)

    errors = diff_pvalues(r.filtered_text, golden_filt)
    assert not errors, "filtered drift:\n" + format_errors(errors)
