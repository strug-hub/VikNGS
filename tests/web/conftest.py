"""Fixtures for the browser test layer.

Ensures the WASM module is built, boots a Vite preview server on a free
port, and tears it all down at end of session.
"""
from __future__ import annotations

import os
import socket
import subprocess
import time
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent.parent
WEB  = REPO / "web"


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _wait_for(url: str, timeout: float = 30.0) -> None:
    import urllib.request
    deadline = time.monotonic() + timeout
    last_err: Exception | None = None
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=1) as resp:
                if resp.status < 500:
                    return
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(0.25)
    raise RuntimeError(f"Vite preview at {url} never came up. Last error: {last_err}")


@pytest.fixture(scope="session")
def wasm_built() -> None:
    need = [WEB / "public" / "vikngs-core.cjs", WEB / "public" / "vikngs-core.wasm"]
    if not all(p.exists() for p in need):
        subprocess.run(["./emcmake-build.sh"], cwd=REPO, check=True)
        # emcmake writes vikngs-core.js; rename to .cjs for the UMD loader
        js = WEB / "public" / "vikngs-core.js"
        cjs = WEB / "public" / "vikngs-core.cjs"
        if js.exists() and not cjs.exists():
            js.rename(cjs)
    # Build the UI bundle once so vite preview has something to serve
    if not (WEB / "dist" / "index.html").exists():
        if not (WEB / "node_modules").exists():
            subprocess.run(["npm", "install"], cwd=WEB, check=True)
        subprocess.run(["npm", "run", "build"], cwd=WEB, check=True)
        # vite build doesn't copy public/ into dist automatically for SPA — it does.
        # If the cjs loader references are off, manual copy:
        dist_pub = WEB / "dist"
        for name in ("vikngs-core.cjs", "vikngs-core.wasm"):
            src = WEB / "public" / name
            dst = dist_pub / name
            if src.exists() and not dst.exists():
                dst.write_bytes(src.read_bytes())


@pytest.fixture(scope="session")
def vite_server(wasm_built) -> str:
    """Boot `vite preview` and return its base URL."""
    port = _free_port()
    env = os.environ.copy()
    proc = subprocess.Popen(
        ["npx", "vite", "preview", "--port", str(port), "--strictPort"],
        cwd=WEB,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    url = f"http://127.0.0.1:{port}"
    try:
        _wait_for(url)
        yield url
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
