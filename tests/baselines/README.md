# Baseline Truth Files

This folder contains persisted “golden” baselines used to detect unintended behavioral changes
after refactors, dependency upgrades, or algorithm changes.

- `mr_baseline_v1.json`: deterministic MR baseline built from a sampled subset of the real GWAS examples (see file for details).

If you intentionally change an algorithm and the baseline should be updated:
1) Regenerate with `genal/tests/generate_baselines.py`
2) Review diffs and document why the baseline changed
