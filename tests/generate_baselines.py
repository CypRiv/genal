"""
Generate baseline truth files for regression testing.

Usage:
  python genal/tests/generate_baselines.py --all
  python genal/tests/generate_baselines.py --liftover
  python genal/tests/generate_baselines.py --prs
  python genal/tests/generate_baselines.py --coloc
  python genal/tests/generate_baselines.py --preprocessing
  python genal/tests/generate_baselines.py --mrpresso

IMPORTANT: Only regenerate baselines after validating outputs are correct.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any


TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
BASELINES_DIR = TESTS_DIR / "baselines"
TESTS_DATA_DIR = TESTS_DIR / "tests_data"


def _ensure_import_paths() -> None:
    # Mirror genal/tests/conftest.py behavior for running as a standalone script.
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    if str(TESTS_DIR) not in sys.path:
        sys.path.insert(0, str(TESTS_DIR))


def _set_isolated_test_home() -> Path:
    test_home = TESTS_DIR / ".genal_test_home"
    test_home.mkdir(parents=True, exist_ok=True)
    os.environ["HOME"] = str(test_home)
    os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    sys.dont_write_bytecode = True
    return test_home


def _write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def _find_plink2(tests_dir: Path) -> str | None:
    plink2 = shutil.which("plink2")
    if plink2:
        return plink2
    candidate = tests_dir / ".genal_test_home" / ".genal" / "plink2" / "plink2"
    if candidate.exists():
        return str(candidate)
    return None


def generate_liftover_baseline() -> None:
    _ensure_import_paths()
    from genal.lift import lift_coordinates_python
    import pandas as pd

    chain_text = (
        "chain 1 chr1 200000 + 0 100000 chr1 200000 + 1000 101000 1\n"
        "100000\n"
    )
    inputs = [
        {"CHR": 1, "POS": 0},
        {"CHR": 1, "POS": 10},
        {"CHR": 1, "POS": 99999},
        {"CHR": 1, "POS": 100000},  # unmapped (outside the block)
        {"CHR": 2, "POS": 10},      # unmapped (no chr2 chain)
    ]

    with tempfile.TemporaryDirectory() as td:
        chain_path = Path(td) / "synthetic.chain"
        chain_path.write_text(chain_text)
        df = pd.DataFrame(inputs)
        out = lift_coordinates_python(df, chain_path=str(chain_path))
        expected = out[["CHR", "POS"]].astype(int).to_dict(orient="records")

    baseline = {
        "version": 1,
        "description": "Synthetic liftover regression baseline (chr1 positions +1000).",
        "source": {"chain_text": chain_text, "input": inputs},
        "expected": expected,
        "tolerances": {"numeric_rel": 0.0, "numeric_abs": 0.0},
    }
    _write_json(BASELINES_DIR / "liftover_baseline_v1.json", baseline)


def generate_prs_baseline() -> None:
    _ensure_import_paths()
    _set_isolated_test_home()

    import numpy as np
    import pandas as pd

    from genal.extract_prs import prs_func
    from genal.tools import default_config, read_config, write_config, load_reference_panel

    plink2 = _find_plink2(TESTS_DIR)
    if not plink2:
        raise RuntimeError("plink2 not found in PATH or under tests/.genal_test_home; cannot generate PRS baseline.")

    cfg = default_config()
    cfg["paths"]["plink2_path"] = plink2
    write_config(cfg)

    geno_prefix = TESTS_DATA_DIR / "StudyA_genetic_files" / "plink_sampled_9"
    ref = load_reference_panel(str(geno_prefix))
    ref = ref[(ref["A1"].str.len() == 1) & (ref["A2"].str.len() == 1)].head(30).reset_index(drop=True)

    weights = pd.DataFrame(
        {
            "SNP": ref["SNP"].astype(str),
            "EA": ref["A1"].astype(str),
            "BETA": np.linspace(-0.2, 0.2, len(ref)),
        }
    )

    old_cwd = Path.cwd()
    with tempfile.TemporaryDirectory() as td:
        os.chdir(td)
        scores = prs_func(weights, weighted=True, path=str(geno_prefix), name="prs_baseline", cpus=1, ram=4000)
    os.chdir(old_cwd)

    # Keep only a small, stable subset in the baseline file.
    scores = scores.sort_values("IID").reset_index(drop=True)
    score_cols = [c for c in scores.columns if c.startswith("SCORE")]
    baseline_rows = scores[["FID", "IID", *score_cols]].head(15).to_dict(orient="records")

    baseline = {
        "version": 1,
        "description": "PRS regression baseline on StudyA chr9 (first 30 biallelic variants).",
        "source": {
            "geno_prefix": "StudyA_genetic_files/plink_sampled_9",
            "n_variants": int(len(weights)),
            "weights": weights.to_dict(orient="records"),
        },
        "expected": {"rows": baseline_rows, "score_columns": score_cols},
        "tolerances": {"numeric_rel": 1e-6, "numeric_abs": 1e-6},
    }
    _write_json(BASELINES_DIR / "prs_baseline_v1.json", baseline)


def generate_colocalization_baseline() -> None:
    _ensure_import_paths()
    import numpy as np
    import pandas as pd

    from genal.colocalization import coloc_abf_func

    rng = np.random.default_rng(0)
    n = 200
    pos = np.arange(1_000_000, 1_000_000 + n)
    beta1 = rng.normal(0.0, 0.05, size=n)
    beta2 = beta1 * 0.8 + rng.normal(0.0, 0.01, size=n)

    data1 = pd.DataFrame(
        {
            "CHR": 1,
            "POS": pos,
            "SNP": [f"rs{i}" for i in range(n)],
            "EA": ["A"] * n,
            "NEA": ["G"] * n,
            "BETA": beta1,
            "SE": 0.1,
            "EAF": rng.uniform(0.05, 0.95, size=n),
        }
    )
    data2 = pd.DataFrame(
        {
            "CHR": 1,
            "POS": pos,
            "SNP": [f"rs{i}" for i in range(n)],
            "EA": ["A"] * n,
            "NEA": ["G"] * n,
            "BETA": beta2,
            "SE": 0.1,
            "EAF": rng.uniform(0.05, 0.95, size=n),
        }
    )

    res = coloc_abf_func(data1.copy(), data2.copy(), trait1_type="quant", trait2_type="quant", sdY1=1, sdY2=1)
    baseline = {
        "version": 1,
        "description": "Colocalization ABF baseline on synthetic correlated traits.",
        "source": {"seed": 0, "n": n},
        "expected": res,
        "tolerances": {"numeric_rel": 1e-10, "numeric_abs": 1e-10},
    }
    _write_json(BASELINES_DIR / "colocalization_baseline_v1.json", baseline)


def generate_preprocessing_baseline() -> None:
    _ensure_import_paths()
    _set_isolated_test_home()

    import numpy as np
    import pandas as pd

    from genal.Geno import Geno
    from genal.tools import default_config, write_config, load_reference_panel

    # Keep config deterministic for baseline generation.
    write_config(default_config())

    ref_prefix = TESTS_DATA_DIR / "StudyA_genetic_files" / "plink_sampled_9"
    ref = load_reference_panel(str(ref_prefix))
    ref = ref[(ref["A1"].str.len() == 1) & (ref["A2"].str.len() == 1)].head(25).reset_index(drop=True)

    raw = pd.DataFrame(
        {
            "CHR": ref["CHR"].astype(int),
            "POS": ref["POS"].astype(int),
            "BETA": np.linspace(-0.1, 0.1, len(ref)),
            "SE": 0.1,
            "P": 0.5,
        }
    )
    g = Geno(raw)
    g.preprocess_data(preprocessing="Fill", reference_panel=str(ref_prefix), fill_snpids=True)

    cols = ["CHR", "POS", "SNP", "EA", "NEA"]
    expected_rows = g.data[cols].head(15).to_dict(orient="records")

    baseline = {
        "version": 1,
        "description": "Preprocessing baseline: Fill SNP/alleles from StudyA chr9 reference.",
        "source": {"reference_panel": "StudyA_genetic_files/plink_sampled_9", "n": int(len(ref))},
        "expected": {"columns": cols, "rows": expected_rows},
        "tolerances": {"numeric_rel": 0.0, "numeric_abs": 0.0},
    }
    _write_json(BASELINES_DIR / "preprocessing_baseline_v1.json", baseline)


def generate_mrpresso_baseline() -> None:
    _ensure_import_paths()
    import numpy as np
    import pandas as pd

    import genal.MRpresso as mrp

    # Run inline to keep deterministic outputs and avoid process spawning in baseline generation.
    class _InlineProcessPool:
        def __init__(self, max_workers: int | None = None) -> None:  # noqa: ARG002
            return

        def __enter__(self):  # noqa: ANN001
            return self

        def __exit__(self, exc_type, exc, tb) -> bool:  # noqa: ANN001
            return False

        def map(self, func, iterable):  # noqa: ANN001
            return list(map(func, iterable))

    mrp.ProcessPoolExecutor = _InlineProcessPool  # type: ignore[attr-defined]
    mrp.tqdm = lambda it, **kwargs: it  # type: ignore[assignment]

    seed_base = 12345
    counter = {"n": 0}

    def _seeded_default_rng():  # noqa: ANN001
        counter["n"] += 1
        return np.random.default_rng(seed_base + counter["n"])

    mrp.default_rng = _seeded_default_rng  # type: ignore[assignment]

    rng = np.random.default_rng(0)
    n = 20
    outlier_idx = 7
    beta_e = rng.normal(0.0, 0.1, size=n)
    beta_o = 0.5 * beta_e + rng.normal(0.0, 0.01, size=n)
    beta_o[outlier_idx] = 5.0
    df = pd.DataFrame({"BETA_o": beta_o, "BETA_e": beta_e, "SE_o": 0.05, "SE_e": 0.05})

    mod_table, global_test, outlier_test, bias_test, _subset = mrp.mr_presso(
        df,
        n_iterations=400,
        outlier_test=True,
        distortion_test=True,
        significance_p=0.05,
        cpus=1,
    )

    baseline = {
        "version": 1,
        "description": "MR-PRESSO baseline on synthetic data with one extreme outlier.",
        "source": {"seed": 0, "n": n, "outlier_index": outlier_idx, "n_iterations": 400},
        "expected": {
            "global_test": global_test,
            "bias_test": bias_test,
            "outlier_test_min_p": float(outlier_test["Pvalue"].min()) if not outlier_test.empty else None,
            "mod_table": mod_table.to_dict(orient="records"),
        },
        "tolerances": {"numeric_rel": 1e-6, "numeric_abs": 1e-6},
    }
    _write_json(BASELINES_DIR / "mrpresso_baseline_v1.json", baseline)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="Regenerate all baselines.")
    parser.add_argument("--liftover", action="store_true", help="Regenerate liftover baseline.")
    parser.add_argument("--prs", action="store_true", help="Regenerate PRS baseline.")
    parser.add_argument("--coloc", action="store_true", help="Regenerate colocalization baseline.")
    parser.add_argument("--preprocessing", action="store_true", help="Regenerate preprocessing baseline.")
    parser.add_argument("--mrpresso", action="store_true", help="Regenerate MR-PRESSO baseline.")

    args = parser.parse_args()
    if not any(vars(args).values()):
        parser.error("No baseline selected. Use --all or one of the specific flags.")

    BASELINES_DIR.mkdir(parents=True, exist_ok=True)

    if args.all or args.liftover:
        generate_liftover_baseline()
    if args.all or args.prs:
        generate_prs_baseline()
    if args.all or args.coloc:
        generate_colocalization_baseline()
    if args.all or args.preprocessing:
        generate_preprocessing_baseline()
    if args.all or args.mrpresso:
        generate_mrpresso_baseline()

    print(f"Wrote baselines under: {BASELINES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
