"""Regenerate golden output files for all e2e cases.

Usage:
    python -m tests.e2e.regen

Writes tests/e2e/golden/<case_id>.pvalues.txt and .filtered.txt.
Only run this when you intentionally want to accept current CLI behavior
as the new baseline.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO))  # so `tests` is importable when run as script

from tests.e2e._harness import run_cli
from tests.e2e.cases import CASES

GOLDEN_DIR = HERE / "golden"
BIN = REPO / "bin" / "vikNGS"
VCF = REPO / "example" / "example.vcf"
INFO = REPO / "example" / "example_info.txt"


def main() -> int:
    if not BIN.exists():
        print(f"ERROR: CLI not built at {BIN}. Run `make` in bin/ first.", file=sys.stderr)
        return 1
    GOLDEN_DIR.mkdir(exist_ok=True)
    for case_id, flags in CASES:
        print(f"[{case_id}] flags={flags}")
        r = run_cli(BIN, VCF, INFO, extra_args=flags)
        if r.returncode != 0:
            print(f"  FAILED: rc={r.returncode}\n  stderr: {r.stderr[:500]}", file=sys.stderr)
            return 2
        (GOLDEN_DIR / f"{case_id}.pvalues.txt").write_text(r.pvalues_text)
        (GOLDEN_DIR / f"{case_id}.filtered.txt").write_text(r.filtered_text)
        print(f"  wrote {len(r.pvalues_text.splitlines())} pvalue lines, {len(r.filtered_text.splitlines())} filtered lines")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
