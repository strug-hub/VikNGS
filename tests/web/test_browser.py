"""End-to-end browser test: drives the Vite-served VikNGS page, uploads the
example files, runs the analysis, scrapes the result rows, and diffs them
against the committed CLI golden (common_default). Proves WASM+UI output
stays numerically equivalent to the native CLI.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

playwright_api = pytest.importorskip("playwright.sync_api")
sync_playwright = playwright_api.sync_playwright


REPO = Path(__file__).resolve().parent.parent.parent
EXAMPLE_VCF   = REPO / "example" / "example.vcf"
EXAMPLE_INFO  = REPO / "example" / "example_info.txt"
GOLDEN_COMMON = REPO / "tests" / "e2e" / "golden" / "common_default.pvalues.txt"


def _parse_golden(text: str) -> list[tuple[str, int, float]]:
    rows = []
    for line in text.strip().splitlines():
        parts = line.split("\t")
        rows.append((parts[0], int(parts[1]), float(parts[4])))
    return rows


def test_example_run_matches_cli_golden(vite_server: str):
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        # Surface console errors in pytest output to aid debugging
        errors: list[str] = []
        page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))
        page.on("console", lambda m: errors.append(f"console[{m.type}]: {m.text}") if m.type == "error" else None)

        page.goto(vite_server, wait_until="domcontentloaded")

        # Upload files via the file inputs.
        page.locator("#f-vcf").set_input_files(str(EXAMPLE_VCF))
        page.locator("#f-sample").set_input_files(str(EXAMPLE_INFO))

        # Default form values match common_default (maf=0.05, depth=30,
        # missing=0.1, mustPass=false, stat=common, gt=expected).
        page.locator("#run-btn").click()

        # Wait for the results table to populate. "Done." appears in the log
        # when the run finishes successfully.
        page.wait_for_selector("#results-table table", timeout=90_000)

        rows = page.locator("#results-table tbody tr").all()
        actual = []
        for r in rows:
            tds = r.locator("td").all_inner_texts()
            actual.append((tds[0], int(tds[1]), float(tds[4])))

        expected = _parse_golden(GOLDEN_COMMON.read_text())

        assert len(actual) == len(expected), (
            f"row count mismatch: actual={len(actual)} expected={len(expected)}. "
            f"page errors so far: {errors[:5]}"
        )
        for i, (a, e) in enumerate(zip(actual, expected)):
            assert a[0] == e[0], f"row {i} chrom: {a[0]!r} vs {e[0]!r}"
            assert a[1] == e[1], f"row {i} pos: {a[1]} vs {e[1]}"
            # Golden p-values are truncated to 6 decimals; allow 5e-6 rtol.
            drift = abs(a[2] - e[2]) / max(abs(e[2]), 1e-30)
            assert drift < 5e-6, f"row {i} pvalue drift {drift}: {a[2]} vs {e[2]}"

        browser.close()

        if errors:
            # Non-fatal but surfaced
            print("page events:\n" + "\n".join(errors))


def test_preload_example_button(vite_server: str):
    """The 'Preload example' header button should fetch the bundled files and
    populate both file inputs; running after that should produce the same
    CLI-golden-matching 16 rows as when uploaded manually.
    """
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(vite_server, wait_until="domcontentloaded")

        page.locator("#preload-example-btn").click()
        # Wait for the ok log line confirming the load finished.
        page.wait_for_selector("#log .ok", timeout=10_000)

        # File inputs should now report a selected file.
        vcf_name = page.eval_on_selector("#f-vcf", "el => el.files[0]?.name")
        sample_name = page.eval_on_selector("#f-sample", "el => el.files[0]?.name")
        assert vcf_name == "example.vcf", f"vcf input: {vcf_name!r}"
        assert sample_name == "example_info.txt", f"sample input: {sample_name!r}"

        page.locator("#run-btn").click()
        page.wait_for_selector("#results-table table", timeout=90_000)
        row_count = page.locator("#results-table tbody tr").count()
        assert row_count == 16, f"expected 16 rows, got {row_count}"

        browser.close()
