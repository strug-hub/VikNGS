"""Shared pytest fixtures for VikNGS test suite."""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN_DIR = REPO_ROOT / "bin"
EXAMPLE_DIR = REPO_ROOT / "example"


def _build(target: str) -> None:
    subprocess.run(
        ["make", target] if target != "all" else ["make"],
        cwd=BIN_DIR,
        check=True,
    )


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def example_vcf() -> Path:
    return EXAMPLE_DIR / "example.vcf"


@pytest.fixture(scope="session")
def example_info() -> Path:
    return EXAMPLE_DIR / "example_info.txt"


@pytest.fixture(scope="session")
def cli_binary() -> Path:
    """Path to the compiled vikNGS CLI binary. Builds it if missing."""
    binary = BIN_DIR / "vikNGS"
    if not binary.exists():
        _build("all")
    assert binary.exists(), f"CLI build did not produce {binary}"
    return binary


@pytest.fixture(scope="session")
def sim_binary() -> Path:
    """Path to the compiled sim-cli binary. Builds it if missing."""
    binary = BIN_DIR / "sim-cli"
    if not binary.exists():
        _build("sim-cli")
    assert binary.exists(), f"sim-cli build did not produce {binary}"
    return binary
