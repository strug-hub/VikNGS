"""Helpers for running sim-cli and parsing its TSV output."""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass
class SimRow:
    variant_idx: int
    test_source: str
    test_stat: str
    pvalue: float


def run_sim(
    sim_binary: Path,
    extra_args: Sequence[str],
    timeout: int = 600,
) -> list[SimRow]:
    """Run sim-cli, return parsed rows. Raises if rc != 0."""
    cmd = [str(sim_binary), "--header", *map(str, extra_args)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(f"sim-cli rc={result.returncode}\nstderr: {result.stderr[:500]}\ncmd: {cmd}")
    rows: list[SimRow] = []
    # printInfo writes to stdout with "[INFO] " prefix, interleaved with our TSV.
    tsv_lines = [ln for ln in result.stdout.splitlines() if ln and not ln.startswith("[")]
    if not tsv_lines:
        return rows
    header = tsv_lines[0].split("\t")
    expected = ["variant_idx", "test_source", "test_stat", "pvalue"]
    if header != expected:
        raise RuntimeError(f"unexpected sim-cli header: {header}")
    for line in tsv_lines[1:]:
        parts = line.split("\t")
        rows.append(SimRow(
            variant_idx=int(parts[0]),
            test_source=parts[1],
            test_stat=parts[2],
            pvalue=float(parts[3]),
        ))
    return rows


def pvalues_by_source(rows: list[SimRow], source: str) -> list[float]:
    return [r.pvalue for r in rows if r.test_source == source]
