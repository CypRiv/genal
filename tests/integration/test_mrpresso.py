from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


class _InlineProcessPool:
    def __init__(self, max_workers: int | None = None) -> None:  # noqa: ARG002
        return

    def __enter__(self) -> "_InlineProcessPool":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False

    def map(self, func, iterable):  # noqa: ANN001
        return list(map(func, iterable))


def _patch_mrpresso_inline(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Make MR-PRESSO deterministic and fast for tests by avoiding process-based parallelism.
    """
    # Avoid spawning processes during tests.
    monkeypatch.setattr("genal.MRpresso.ProcessPoolExecutor", _InlineProcessPool)

    # Silence tqdm overhead.
    monkeypatch.setattr("genal.MRpresso.tqdm", lambda it, **kwargs: it)

    # Deterministic-but-varying RNG per call.
    seed_base = 12345
    counter = {"n": 0}

    def _seeded_default_rng():  # noqa: ANN001
        counter["n"] += 1
        return np.random.default_rng(seed_base + counter["n"])

    monkeypatch.setattr("genal.MRpresso.default_rng", _seeded_default_rng)


def _synthetic_mrpresso_df(n: int = 10, *, outlier_index: int | None = None) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    beta_e = rng.normal(0.0, 0.1, size=n)
    se_e = np.full(n, 0.05)
    se_o = np.full(n, 0.05)

    beta_o = 0.5 * beta_e + rng.normal(0.0, 0.01, size=n)
    if outlier_index is not None:
        beta_o[outlier_index] = 5.0  # extreme outlier

    return pd.DataFrame({"BETA_o": beta_o, "BETA_e": beta_e, "SE_o": se_o, "SE_e": se_e})


def test_mr_presso_validates_minimum_ivs() -> None:
    from genal.MRpresso import mr_presso

    df = _synthetic_mrpresso_df(n=3)
    with pytest.raises(Exception, match="Not enough instrumental variables"):
        mr_presso(df, n_iterations=100, cpus=1)


def test_mr_presso_validates_n_iterations() -> None:
    from genal.MRpresso import mr_presso

    df = _synthetic_mrpresso_df(n=10)
    with pytest.raises(Exception, match="increase n_iterations"):
        mr_presso(df, n_iterations=10, cpus=1)


def test_power_eigen_basic() -> None:
    from genal.MRpresso import power_eigen

    x = np.diag([2.0, 3.0])
    out = power_eigen(x, 2)
    assert out[0, 0] == pytest.approx(4.0)
    assert out[1, 1] == pytest.approx(9.0)


def test_getRSS_LOO_single_exposure() -> None:
    from genal.MRpresso import getRSS_LOO

    beta_e = np.array([0.1, 0.2, 0.3, 0.4])
    df = pd.DataFrame(
        {
            "BETA_o": 2.0 * beta_e,
            "BETA_e": beta_e,
            "Weights": np.ones_like(beta_e),
        }
    )
    rss, estimates = getRSS_LOO(df, ["BETA_e"], returnIV=True)
    assert rss == pytest.approx(0.0, abs=1e-10)
    assert np.allclose(estimates, 2.0, atol=1e-10)


def test_mr_presso_global_test_structure(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.MRpresso import mr_presso

    _patch_mrpresso_inline(monkeypatch)
    df = _synthetic_mrpresso_df(n=12)
    mod_table, global_test, outlier_test, bias_test, subset_data = mr_presso(
        df,
        n_iterations=200,
        outlier_test=True,
        distortion_test=True,
        cpus=1,
    )

    assert {"method", "nSNP", "b", "se", "pval"}.issubset(mod_table.columns)
    assert isinstance(global_test, dict)
    assert "global_test_p" in global_test
    assert isinstance(outlier_test, pd.DataFrame)
    assert isinstance(bias_test, dict)
    assert subset_data is None or isinstance(subset_data, pd.DataFrame)


def test_mr_presso_outlier_test_disabled(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.MRpresso import mr_presso

    _patch_mrpresso_inline(monkeypatch)
    df = _synthetic_mrpresso_df(n=12)
    _, global_test, outlier_test, bias_test, subset_data = mr_presso(
        df,
        n_iterations=200,
        outlier_test=False,
        distortion_test=True,
        cpus=1,
    )

    assert "global_test_p" in global_test
    assert outlier_test.empty
    assert bias_test == {}
    assert subset_data is None


@pytest.mark.slow
def test_mr_presso_with_synthetic_outliers(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.MRpresso import mr_presso

    _patch_mrpresso_inline(monkeypatch)

    outlier_idx = 7
    df = _synthetic_mrpresso_df(n=20, outlier_index=outlier_idx)
    _, global_test, outlier_test, bias_test, subset_data = mr_presso(
        df,
        n_iterations=400,
        outlier_test=True,
        distortion_test=True,
        significance_p=0.05,
        cpus=1,
    )

    assert global_test["global_test_p"] < 0.05
    assert not outlier_test.empty
    assert isinstance(bias_test.get("outliers_indices"), list)
    assert outlier_idx in bias_test["outliers_indices"]
    assert subset_data is not None
    assert len(subset_data) < len(df)

