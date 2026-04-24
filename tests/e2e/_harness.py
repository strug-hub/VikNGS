"""Run the VikNGS CLI and capture its output files.

CLI writes timestamped files (pvalues_YYYY-MM-DD_HH-MM-SS-NN.txt, filtered_...txt)
to the directory passed via -o. This module runs the binary in a fresh tmpdir
and returns the contents by file-type.
"""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass
class CliRun:
    returncode: int
    stdout: str
    stderr: str
    pvalues_text: str
    filtered_text: str


def run_cli(
    cli_binary: Path,
    vcf: Path,
    sample: Path,
    extra_args: Sequence[str] = (),
    timeout: int = 120,
) -> CliRun:
    with tempfile.TemporaryDirectory(prefix="vikngs-e2e-") as tmp:
        tmp_path = Path(tmp)
        cmd = [
            str(cli_binary),
            "--vcf", str(vcf),
            "--sample", str(sample),
            "-o", str(tmp_path),
            *extra_args,
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        pvalues = _read_unique(tmp_path, "pvalues_*.txt")
        filtered = _read_unique(tmp_path, "filtered_*.txt")
        return CliRun(
            returncode=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            pvalues_text=pvalues,
            filtered_text=filtered,
        )


def _read_unique(directory: Path, pattern: str) -> str:
    hits = sorted(directory.glob(pattern))
    if not hits:
        return ""
    if len(hits) > 1:
        raise AssertionError(f"Expected 1 file matching {pattern} in {directory}, got {len(hits)}")
    return hits[0].read_text()
