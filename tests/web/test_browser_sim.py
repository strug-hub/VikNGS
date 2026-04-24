"""Browser regression test for the simulation tab. Drives the page,
switches to the Simulation tab, runs a small seeded simulation, and
sanity-checks the derived power-curve data exposed on window.__simData.

The WASM sim's RNG diverges from sim-cli's (libc++ vs libstdc++), so we
can't bit-diff p-values. Instead we assert determinism (same seed
reproduces identical __simData) and loose null-behavior (power < alpha
within sampling noise at alpha=0.05 with NSNP=200).
"""
from __future__ import annotations

import pytest

playwright_api = pytest.importorskip("playwright.sync_api")
sync_playwright = playwright_api.sync_playwright


def _run_sim(page, seed: int):
    page.locator('nav.tabs .tab[data-tab="simulation"]').click()
    page.wait_for_selector("#sim-form")

    # Expand every <details> so the seed input is visible to Playwright's
    # visibility check. The form has several collapsed groups by default.
    page.evaluate("document.querySelectorAll('#sim-form details').forEach(d => d.open = true)")

    # Default sim form already has nsnp=200, effectSize=1, 2 groups of 500.
    # Override seed so runs are reproducible.
    page.fill("#sf-seed", str(seed))

    page.locator("#sim-run-btn").click()
    # Power plot appears when render completes.
    page.wait_for_selector("#sim-power .u-over", timeout=120_000)

    return page.evaluate("() => window.__simData")


def test_sim_runs_and_produces_power_curve(vite_server: str):
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.on("pageerror", lambda e: print(f"pageerror: {e}"))
        page.goto(vite_server, wait_until="domcontentloaded")

        data = _run_sim(page, seed=1337)

        assert data is not None, "window.__simData was not exposed"
        assert len(data["keys"]) >= 2, f"expected >=2 test keys, got {data['keys']}"
        assert data["steps"] == 1
        assert data["rowCount"] > 0

        # Under the null (effectSize=1), alpha=0.05 → type-I rate close to
        # 0.05 per test. With nsnp=200 the 99% CI is roughly [0.01, 0.10].
        for key, powers in zip(data["keys"], data["power"]):
            assert len(powers) == 1
            rate = powers[0]
            assert 0.0 <= rate <= 0.25, (
                f"type-I rate for {key} = {rate} looks far from alpha=0.05"
            )

        browser.close()


def test_sim_is_deterministic_with_seed(vite_server: str):
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(vite_server, wait_until="domcontentloaded")

        a = _run_sim(page, seed=42)
        # Second run with same seed.
        b = _run_sim(page, seed=42)

        assert a["keys"] == b["keys"]
        assert a["power"] == b["power"], (
            f"seeded reruns diverged: {a['power']} vs {b['power']}"
        )

        browser.close()
