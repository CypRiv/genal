# genal-python Test Suite

## Overview
This folder contains the pytest suite for `genal`:
- **Unit tests** (`unit/`): fast, isolated tests of individual functions and edge cases.
- **Integration tests** (`integration/`): multi-component workflows (may touch files, PLINK2, larger fixtures).
- **Baselines** (`baselines/`): regression “golden” outputs to detect unintended behavior changes.
- **Test data** (`tests_data/`): real GWAS examples + StudyA genetic/phenotype files used by tests.

The suite is designed to be **side-effect isolated**:
- `tests/conftest.py` forces `HOME` under `tests/.genal_test_home/<worker>/` (so `~/.genal` stays in the repo).
- Each test runs from its own temporary working directory (so `tmp_GENAL/` and PLINK outputs don’t pollute the repo).
- Temporary test artifacts are cleaned between runs (see “Temporary Files & Cleanup”).

## Quick Start
Run from the git repo root (`genal/`). If you're at the workspace root, `cd genal` first.

- Run everything:
  - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests -v`
- Skip slow tests:
  - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests -v -m "not slow"`
- Run network tests (real GWAS Catalog calls):
  - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests -v --run-network -m network`
- Full run including network tests:
  - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests -v --run-network`
- Run only PLINK2 tests:
  - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests -v -m plink`
- Generate a run report (JUnit + stdout + markdown summary under `tests/reports/`):
  - `conda run -n genal_test python tests/run_suite.py --run-network`

## Test Structure
```
tests/
├── unit/           # fast, isolated function tests
├── integration/    # multi-component workflow tests
├── baselines/      # persisted regression reference data (goldens)
├── tests_data/     # large GWAS + StudyA genetic test files
├── reports/        # test run artifacts (junit, stdout)
├── conftest.py     # shared fixtures + environment isolation
├── pytest.ini      # pytest configuration + markers
└── run_suite.py    # local runner that writes a report under reports/
```

## Test Data
Primary files live under `tests/tests_data/`:
- `GWAS_example1_b37.txt` (GRCh37, whitespace-delimited)
- `GWAS_example2_b38.csv` (GRCh38, whitespace-delimited despite extension)
- `GWAS_example3_b38.parquet` (GRCh38, parquet)
- `GWAS_example4_b37.tsv` (GRCh37, TSV)
- `StudyA_genetic_files/` (PLINK2-compatible reference data used for clumping/extraction/PRS/association)
- `StudyA_phenotype_example.csv` (phenotype + covariates)

Helpers in `tests/gwas_examples.py` provide standardized, row-limited loaders for speed (they never load the full multi-GB files in unit tests).

## Pytest Markers
- `slow`: longer-running tests (still offline).
- `network`: uses the network (real GWAS Catalog API). Skipped by default unless `--run-network`.
- `plink`: requires a `plink2` executable. Locally, tests will also pick up `tests/.genal_test_home/.genal/plink2/plink2` if present.

`tests/pytest.ini` enables xdist by default (`-n auto`). For debugging, disable parallelism with `-n 0`:
- `conda run -n genal_test python -m pytest -c tests/pytest.ini tests -n 0 -vv -s`

## Baseline Truth Testing
Baselines are JSON files capturing “known good” outputs for key workflows (MR, liftover, PRS, coloc, preprocessing, MR-PRESSO, …).

### What counts as “drift”?
Any change in results compared to `tests/baselines/*.json` is considered drift. Drift can be:
- a real regression (bug),
- an intentional algorithmic change (baseline should be updated), or
- a non-determinism issue (same code yields different results due to randomness/parallelism).

### Run baseline checks
- Run just the baseline assertions:
  - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests/integration/test_baselines.py -v`
- Include MR ground-truth (MR baseline lives in `tests/baselines/mr_baseline_v1.json` but is asserted in `tests/integration/test_mr_ground_truth.py`):
  - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests/integration/test_mr_ground_truth.py -v`

### Investigate a baseline failure (drift workflow)
1) Re-run the failing test in serial with full output:
   - `conda run -n genal_test python -m pytest -c tests/pytest.ini tests/integration/test_baselines.py -n 0 -vv -s`
2) Open the referenced baseline JSON under `tests/baselines/` and check:
   - `source`: inputs/parameters used to generate the baseline
   - `expected`: the stored “golden” outputs
   - `tolerances`: numeric tolerances used in assertions
3) Decide whether the change is expected:
   - If the change is **not expected**, treat it as a regression: fix the code, then re-run the baseline test.
   - If the change is **expected**, regenerate the baseline (next section), review diffs, and bump the baseline version if appropriate.

### Regenerate baselines (intentional changes only)
Baselines are generated by `tests/generate_baselines.py`.

- Regenerate all baselines:
  - `conda run -n genal_test python tests/generate_baselines.py --all`
- Regenerate a single baseline:
  - `conda run -n genal_test python tests/generate_baselines.py --mrpresso`

After regenerating:
- Re-run the baseline tests.
- Review the JSON diff (`git diff tests/baselines/`) and validate the changes are intended.

### Handling randomness in baseline comparisons
Some algorithms have inherent randomness (bootstrapping, MR-PRESSO permutations) or can behave non-deterministically under parallelism.

Baseline strategy in this repo:
- Prefer **deterministic baselines**: store explicit seeds in baseline `source` and force deterministic RNG in tests when needed.
  - Example: MR-PRESSO baseline tests monkeypatch the RNG and run inline (no process pool) to remove parallel non-determinism.
- Use **tolerances** for floating-point comparisons (`numeric_rel` / `numeric_abs`) instead of exact equality.
- For highly stochastic methods, baseline only **stable summaries** (e.g., min p-value, outlier indices) rather than full per-iteration traces.
- For bootstrap-heavy MR methods, prefer “property tests” (e.g., estimate is within a wide tolerance of the known causal effect) rather than baseline exact coefficients.

## Adding New Tests
### Where to place tests
- Prefer `unit/` when the test can be pure, small, and fast (no PLINK, no network, minimal I/O).
- Use `integration/` when calling `Geno` methods, touching disk, running PLINK2, or exercising multi-module workflows.

### Preserve test isolation
- Don’t write to the repo root. The suite changes CWD to a per-test temp dir; rely on that (or use `tmp_path`) for any file outputs.
- If you call functions that expect `tmp_GENAL/`, call `genal.tools.create_tmp()` after CWD isolation.
- If you need to inspect artifacts during debugging, run with `--keep-tmp` (see next section).

### Markers (required conventions)
- Mark expensive tests as `@pytest.mark.slow`.
- Mark real network/API tests as `@pytest.mark.network` (and ensure they run only with `--run-network`).
- Mark PLINK2-backed tests as `@pytest.mark.plink` (tests will skip automatically if `plink2` can’t be found).

### Adding a new regression baseline
When a module’s output should not drift silently:
1) Add a baseline generator to `tests/generate_baselines.py`.
2) Write the baseline JSON to `tests/baselines/`.
3) Add an assertion test in `tests/integration/test_baselines.py`.

Guidelines for stable baselines:
- Store only a **small subset** of rows/outputs (first N rows, key summary metrics).
- Include all parameters and seeds under `source`.
- For stochastic algorithms, enforce determinism (fixed seeds / monkeypatch RNG / avoid nondeterministic parallelism) and/or compare only stable summaries with tolerances.

## Fixtures Reference
See `tests/conftest.py` for the authoritative list. Common fixtures include:
- `tests_data_dir`, `studyA_genetic_dir`, `phenotype_csv_path`
- `reset_genal_config` (isolated `~/.genal/config.json` under the test HOME)

## Temporary Files & Cleanup
Temporary outputs are written only under:
- `tests/.tmp/` (pytest temp dirs; includes any `tmp_GENAL/` created during tests)
- `tests/.genal_test_home/` (isolated `~/.genal` config + small caches)

By default, `tests/.tmp/` and any repo-root `tmp_GENAL/` are cleaned at the start/end of each pytest run.
To keep them for debugging:
- `conda run -n genal_test python -m pytest -c tests/pytest.ini tests --keep-tmp -n 0 -vv -s`

## Troubleshooting
- **Network tests skipped**: run with `--run-network`.
- **PLINK2 tests skipped**: ensure `plink2` is in `PATH` or install it (CI should provide it).
- **HDF5 outcome test skipped**: install `tables`/PyTables in the test env if you want to exercise `.h5` I/O.
- **Unexpected files in repo**: tests should not write outside their temp dir; check for missing `monkeypatch.chdir()` usage.

## CI/CD Notes
- CI should run `conda run -n genal_test ...` (no new env creation).
- Artifacts (JUnit + stdout) are written by `tests/run_suite.py` into `tests/reports/`.
