"""Test matrix: (case_id, CLI flag list) pairs used by both regen and test_cli."""
from __future__ import annotations

CASES: list[tuple[str, list[str]]] = [
    ("common_default",        []),
    ("common_call",           ["-g", "call"]),
    ("common_vcf",            ["-g", "vcf"]),
    ("rare_cast",             ["-r", "cast"]),
    ("rare_cast_call",        ["-r", "cast", "-g", "call"]),
    ("rare_skat",             ["-r", "skat"]),
    ("rare_skat_call",        ["-r", "skat", "-g", "call"]),
    ("maf_tight",             ["-m", "0.01"]),
    ("depth_low",             ["-d", "10"]),
    ("missing_strict",        ["-x", "0.05"]),
    ("region_filter",         ["--chr", "11", "--from", "193000", "--to", "195000"]),
]
