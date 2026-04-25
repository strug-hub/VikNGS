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


def test_family_toggle_shows_phenotype_and_covariate(vite_server: str):
    """Switching family=normal should reveal the Phenotype and Covariate
    groups (hidden under binomial), and the cohort dropdown should switch
    to a single 'normal' option."""
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(vite_server, wait_until="domcontentloaded")
        page.locator('nav.tabs .tab[data-tab="simulation"]').click()
        page.wait_for_selector("#sim-form")

        # Binomial (default): phenotype + covariate groups are hidden.
        assert page.evaluate("() => document.querySelector('#sf-phenotypeMean').closest('details').hidden") is True
        assert page.evaluate("() => document.querySelector('#sf-covariate').closest('details').hidden") is True

        # Switch family to normal.
        page.select_option("#sf-family", "normal")

        # Now both groups should be visible.
        assert page.evaluate("() => document.querySelector('#sf-phenotypeMean').closest('details').hidden") is False
        assert page.evaluate("() => document.querySelector('#sf-covariate').closest('details').hidden") is False

        # Cohort dropdowns in the groups table should now only offer 'normal'.
        cohorts = page.evaluate("""() => Array.from(document.querySelectorAll('table.sim-groups tbody tr td:nth-child(2) select'))
            .map(s => Array.from(s.options).map(o => o.value))""")
        assert cohorts, "no cohort selects found"
        for opts in cohorts:
            assert opts == ["normal"], f"unexpected cohort options under normal family: {opts}"

        browser.close()


def test_normal_preset_runs(vite_server: str):
    """The 'Quant null' preset switches to normal family and produces a
    valid power curve."""
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.on("pageerror", lambda e: print(f"pageerror: {e}"))
        page.goto(vite_server, wait_until="domcontentloaded")
        page.locator('nav.tabs .tab[data-tab="simulation"]').click()
        page.wait_for_selector("#sim-form")

        # Click the "Quant null" preset.
        page.locator("button[title^='Quantitative null']").click()
        # Make seed input visible to override.
        page.evaluate("document.querySelectorAll('#sim-form details').forEach(d => d.open = true)")
        page.fill("#sf-seed", "7")
        # Reduce work for a fast test.
        page.fill("#sf-nsnp", "100")

        page.locator("#sim-run-btn").click()
        page.wait_for_selector("#sim-power .u-over", timeout=180_000)

        data = page.evaluate("() => window.__simData")
        assert data["rowCount"] > 0
        # Sample-info table should populate with one 'normal' cohort.
        cohort_text = page.locator("#sim-sample-info table tbody tr").first.inner_text()
        assert "normal" in cohort_text.lower(), f"sample table missing normal cohort: {cohort_text!r}"

        browser.close()


def test_sample_info_renders_after_run(vite_server: str):
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(vite_server, wait_until="domcontentloaded")

        _run_sim(page, seed=99)
        # First (default-binomial) sample table should have header + ≥2 rows.
        n_rows = page.locator("#sim-sample-info table tbody tr").count()
        assert n_rows >= 2, f"expected ≥2 group rows in sample table, got {n_rows}"

        browser.close()


def test_html_export_is_self_contained_interactive(vite_server: str):
    """Exported report should be self-contained: inlines uPlot, the
    viewer script, and the raw p-value rows so the file is fully
    interactive offline."""
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(vite_server, wait_until="domcontentloaded")

        _run_sim(page, seed=11)
        assert page.locator("#sim-export-btn").is_enabled(), "Export button should be enabled after run"

        with page.expect_download(timeout=15_000) as download_info:
            page.locator("#sim-export-btn").click()
        download = download_info.value
        assert download.suggested_filename.endswith(".html"), download.suggested_filename
        path = download.path()
        assert path is not None, "download has no path"
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()

        assert "<!doctype html>" in content.lower()
        assert "VikNGS simulation report" in content
        # uPlot inlined.
        assert "uPlot" in content
        # Viewer module inlined and exposes VikngsReport.render.
        assert "VikngsReport" in content
        # Raw rows are embedded so the file actually works offline.
        assert '"rows":' in content
        assert '"pvalue":' in content
        # Reproducibility block.
        assert "Run configuration" in content

        # Open the exported file in the browser and verify it actually
        # produces an interactive plot — uPlot's container class appears.
        # Use file:// load via context.new_page so the file is offline-rendered.
        # Copy the downloaded file to a stable path; the playwright
        # artifacts dir is locked to the originating context.
        import tempfile, shutil
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as tf:
            stable_path = tf.name
        shutil.copy(str(path), stable_path)

        ctx = browser.new_context()
        offline_page = ctx.new_page()
        offline_page.goto("file://" + stable_path, wait_until="load")
        offline_page.wait_for_selector("#r-power .u-over", timeout=15_000)
        # Sample table rendered.
        n_rows = offline_page.locator("#r-sample table tbody tr").count()
        assert n_rows >= 1, f"expected sample rows in offline report, got {n_rows}"
        ctx.close()

        browser.close()
