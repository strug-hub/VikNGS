"""Smoke test: verify CLI builds and --help exits 0. Validates pytest wiring."""
from __future__ import annotations

import subprocess


def test_cli_help(cli_binary):
    result = subprocess.run(
        [str(cli_binary), "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert "vikNGS" in result.stdout
    assert "--vcf" in result.stdout


def test_example_files_exist(example_vcf, example_info):
    assert example_vcf.exists()
    assert example_info.exists()
