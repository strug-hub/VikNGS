# VikNGS Test Suite

Three layers of tests, chosen for different jobs:

| Layer | Location | Runner | Purpose | Speed |
|-------|----------|--------|---------|-------|
| **e2e** | `tests/e2e/` | pytest | Regression safety: diff CLI output against golden files | ~1s |
| **stats** | `tests/stats/` | pytest + scipy | Statistical correctness: simulation-driven Type I / power / vRVS-claim checks | ~5s |
| **unit** | `tests/unit/` | Catch2 (vendored) | Fast micro-tests on individual math primitives | <1s |

## Quick start

```bash
# One-shot: builds everything, runs everything
make -C tests build
make -C tests test
```

Individual layers:

```bash
make -C tests e2e     # golden-file regression
make -C tests stats   # simulation-driven
make -C tests unit    # Catch2
```

## Requirements

- `g++` with C++17 support
- Python 3 with `pytest`, `numpy`, `scipy` (install via `pip install -r tests/requirements.txt`)

## Binaries built

- `bin/vikNGS` — main CLI (existing)
- `bin/sim-cli` — headless simulation driver (built from `src/cmd/simcmd.cpp`)
- `bin/test-unit` — Catch2 unit-test binary

All three go in `bin/` and are gitignored.

## E2E: regression tests

The test matrix is defined in `tests/e2e/cases.py`. Each case runs the CLI with a flag combination against `example/example.vcf` and diffs the result against a committed golden file in `tests/e2e/golden/`.

The diff (`tests/e2e/_diff.py`) is tolerant on numeric columns (`rtol=1e-9`) and exact on chr/pos/ref/alt/test-name.

### Add a new e2e case

1. Add `(case_id, [flags...])` to `CASES` in `tests/e2e/cases.py`.
2. Generate its golden: `make -C tests regen CONFIRM=1`.
3. Commit both the new case and its goldens.

### Regenerate all goldens

Only do this when you've intentionally changed CLI behavior and want to accept the new output as the baseline. It **overwrites** every file in `tests/e2e/golden/`.

```bash
make -C tests regen CONFIRM=1
```

## Stats: statistical correctness

These run `bin/sim-cli` (a headless wrapper around the existing `startSimulation` engine) and apply statistical checks to the resulting p-values.

- `test_sim_smoke.py` — plumbing
- `test_type1.py` — null simulations; p-values should be uniform and the α=0.05 rejection rate in `[0.03, 0.07]`
- `test_power.py` — alternative-hypothesis simulations; rejection rate exceeds a pinned floor
- `test_vrvs_claim.py` — under depth-imbalanced null data, `expected` (vRVS) stays calibrated while `call` inflates. This test validates the paper's central claim.

Statistical tests use fresh random seeds on each run. Tolerances were chosen so false-failure probability per run is <0.1%.

## Unit: Catch2 micro-tests

Tests for individual math and genetics primitives. Runs sub-second.

`tests/unit/catch_amalgamated.{hpp,cpp}` is vendored Catch2 v3.5.2 (Boost license). To update: replace both files from a newer release.

Adding a test: drop a new `test_*.cpp` in `tests/unit/`, it'll be picked up automatically by the Makefile.
