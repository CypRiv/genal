from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _baselines_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "baselines"


def test_liftover_baseline_matches() -> None:
    from genal.lift import lift_coordinates_python

    baseline = json.loads((_baselines_dir() / "liftover_baseline_v1.json").read_text())
    chain_text = baseline["source"]["chain_text"]
    inputs = baseline["source"]["input"]

    with tempfile.TemporaryDirectory() as td:
        chain_path = Path(td) / "synthetic.chain"
        chain_path.write_text(chain_text)
        df = pd.DataFrame(inputs)
        out = lift_coordinates_python(df, chain_path=str(chain_path))

    expected = baseline["expected"]
    got = out[["CHR", "POS"]].astype(int).to_dict(orient="records")
    assert got == expected


@pytest.mark.slow
@pytest.mark.plink
def test_prs_baseline_matches(reset_genal_config, tests_dir: Path, tests_data_dir: Path) -> None:
    import os

    from genal.extract_prs import prs_func
    from genal.tools import read_config, write_config, create_tmp

    baseline = json.loads((_baselines_dir() / "prs_baseline_v1.json").read_text())

    plink2 = shutil.which("plink2")
    if not plink2:
        candidate = tests_dir / ".genal_test_home" / ".genal" / "plink2" / "plink2"
        if candidate.exists():
            plink2 = str(candidate)
    if not plink2:
        pytest.skip("plink2 not available")

    cfg = read_config()
    cfg["paths"]["plink2_path"] = plink2
    write_config(cfg)

    geno_prefix = tests_data_dir / baseline["source"]["geno_prefix"]
    weights = pd.DataFrame(baseline["source"]["weights"])

    with tempfile.TemporaryDirectory() as td:
        old_cwd = Path.cwd()
        try:
            # Keep tmp_GENAL and PLINK outputs out of the repo tree.
            os.chdir(td)
            create_tmp()
            scores = prs_func(weights, weighted=True, path=str(geno_prefix), name="prs_baseline_test", cpus=1, ram=4000)
        finally:
            os.chdir(old_cwd)

    scores = scores.sort_values("IID").reset_index(drop=True)
    score_cols = baseline["expected"]["score_columns"]
    got_rows = scores[["FID", "IID", *score_cols]].head(len(baseline["expected"]["rows"])).to_dict(orient="records")

    expected_rows = baseline["expected"]["rows"]
    tol = baseline["tolerances"]
    for got, exp in zip(got_rows, expected_rows):
        assert got["FID"] == exp["FID"]
        assert got["IID"] == exp["IID"]
        for col in score_cols:
            assert float(got[col]) == pytest.approx(float(exp[col]), rel=tol["numeric_rel"], abs=tol["numeric_abs"])


def test_colocalization_baseline_matches() -> None:
    from genal.colocalization import coloc_abf_func

    baseline = json.loads((_baselines_dir() / "colocalization_baseline_v1.json").read_text())
    seed = int(baseline["source"]["seed"])
    n = int(baseline["source"]["n"])

    rng = np.random.default_rng(seed)
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

    got = coloc_abf_func(data1.copy(), data2.copy(), trait1_type="quant", trait2_type="quant", sdY1=1, sdY2=1)
    expected = baseline["expected"]
    tol = baseline["tolerances"]
    for k, v in expected.items():
        if isinstance(v, (int, float)):
            assert float(got[k]) == pytest.approx(float(v), rel=tol["numeric_rel"], abs=tol["numeric_abs"])
        else:
            assert got[k] == v


def test_preprocessing_baseline_matches(reset_genal_config, tests_data_dir: Path) -> None:
    from genal.Geno import Geno

    baseline = json.loads((_baselines_dir() / "preprocessing_baseline_v1.json").read_text())
    ref_prefix = tests_data_dir / baseline["source"]["reference_panel"]

    # Recreate the same raw input (head of the reference panel).
    from genal.tools import load_reference_panel

    ref = load_reference_panel(str(ref_prefix))
    ref = ref[(ref["A1"].str.len() == 1) & (ref["A2"].str.len() == 1)].head(int(baseline["source"]["n"])).reset_index(drop=True)

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

    cols = baseline["expected"]["columns"]
    got_rows = g.data[cols].head(len(baseline["expected"]["rows"])).to_dict(orient="records")
    assert got_rows == baseline["expected"]["rows"]


def test_mrpresso_baseline_matches(monkeypatch: pytest.MonkeyPatch) -> None:
    import genal.MRpresso as mrp

    baseline = json.loads((_baselines_dir() / "mrpresso_baseline_v1.json").read_text())

    # Match the baseline generation strategy (inline execution + seeded RNG per call).
    class _InlineProcessPool:
        def __init__(self, max_workers: int | None = None) -> None:  # noqa: ARG002
            return

        def __enter__(self):  # noqa: ANN001
            return self

        def __exit__(self, exc_type, exc, tb) -> bool:  # noqa: ANN001
            return False

        def map(self, func, iterable):  # noqa: ANN001
            return list(map(func, iterable))

    monkeypatch.setattr(mrp, "ProcessPoolExecutor", _InlineProcessPool)
    monkeypatch.setattr(mrp, "tqdm", lambda it, **kwargs: it)

    seed_base = 12345
    counter = {"n": 0}

    def _seeded_default_rng():  # noqa: ANN001
        counter["n"] += 1
        return np.random.default_rng(seed_base + counter["n"])

    monkeypatch.setattr(mrp, "default_rng", _seeded_default_rng)

    rng = np.random.default_rng(int(baseline["source"]["seed"]))
    n = int(baseline["source"]["n"])
    outlier_idx = int(baseline["source"]["outlier_index"])
    beta_e = rng.normal(0.0, 0.1, size=n)
    beta_o = 0.5 * beta_e + rng.normal(0.0, 0.01, size=n)
    beta_o[outlier_idx] = 5.0
    df = pd.DataFrame({"BETA_o": beta_o, "BETA_e": beta_e, "SE_o": 0.05, "SE_e": 0.05})

    mod_table, global_test, outlier_test, bias_test, _subset = mrp.mr_presso(
        df,
        n_iterations=int(baseline["source"]["n_iterations"]),
        outlier_test=True,
        distortion_test=True,
        significance_p=0.05,
        cpus=1,
    )

    expected = baseline["expected"]
    tol = baseline["tolerances"]

    assert float(global_test["global_test_p"]) == pytest.approx(float(expected["global_test"]["global_test_p"]), rel=tol["numeric_rel"], abs=tol["numeric_abs"])
    assert bias_test.get("outliers_indices") == expected["bias_test"].get("outliers_indices")
    got_rows = mod_table.to_dict(orient="records")
    exp_rows = expected["mod_table"]
    assert len(got_rows) == len(exp_rows)
    for got, exp in zip(got_rows, exp_rows):
        for k, v in exp.items():
            if isinstance(v, (int, float)) and v == v:  # not NaN
                assert float(got[k]) == pytest.approx(float(v), rel=tol["numeric_rel"], abs=tol["numeric_abs"])
            else:
                assert got[k] == v
