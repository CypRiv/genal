from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_weighted_median_matches_expected() -> None:
    from genal.MR import weighted_median

    b_iv = np.array([1.0, 2.0, 3.0, 4.0])
    weights = np.array([0.1, 0.2, 0.3, 0.4])
    b = weighted_median(b_iv, weights)
    assert b == pytest.approx(3.142857142857143, rel=1e-12)


def test_mr_ivw_matches_closed_form_slope() -> None:
    from genal.MR import mr_ivw

    beta_e = pd.Series([0.10, 0.20, 0.30, 0.40])
    beta_o = pd.Series([0.12, 0.18, 0.33, 0.41])
    se_o = pd.Series([0.05, 0.05, 0.05, 0.05])
    se_e = pd.Series([0.01, 0.01, 0.01, 0.01])

    res = mr_ivw(beta_e, se_e, beta_o, se_o)[0]
    w = 1 / (se_o**2)
    expected_b = float((w * beta_e * beta_o).sum() / (w * (beta_e**2)).sum())
    assert res["b"] == pytest.approx(expected_b, rel=1e-12)
    assert res["se"] > 0
    assert 0 <= res["pval"] <= 1
    assert res["nSNP"] == len(beta_e)


def test_mr_egger_regression_agrees_with_statsmodels() -> None:
    import statsmodels.api as sm
    from scipy.stats import t, chi2
    from genal.MR import mr_egger_regression

    beta_e = pd.Series([0.10, -0.20, 0.30, 0.40])
    beta_o = pd.Series([0.12, -0.10, 0.33, 0.41])
    se_o = pd.Series([0.05, 0.05, 0.05, 0.05])
    se_e = pd.Series([0.01, 0.01, 0.01, 0.01])

    out = mr_egger_regression(beta_e, se_e, beta_o, se_o)
    slope = out[0]
    intercept = out[1]

    sign0 = np.sign(beta_e.replace(0, 1))
    y = (beta_o * sign0).copy()
    x = np.abs(beta_e.copy())
    X = sm.add_constant(x)
    w = 1 / (se_o**2)
    mod = sm.WLS(y, X, weights=w).fit()

    expected_b = float(mod.params.iloc[1])
    expected_se = float(mod.bse.iloc[1] / min(1, np.sqrt(mod.mse_resid)))
    expected_p = float(2 * t.sf(abs(expected_b / expected_se), len(beta_e) - 2))

    expected_bi = float(mod.params.iloc[0])
    expected_sei = float(mod.bse.iloc[0] / min(1, np.sqrt(mod.mse_resid)))
    expected_pi = float(2 * t.sf(abs(expected_bi / expected_sei), len(beta_e) - 2))

    expected_Q = float(mod.mse_resid * (len(beta_e) - 2))
    expected_Qp = float(chi2.sf(expected_Q, len(beta_e) - 2))

    assert slope["b"] == pytest.approx(expected_b, rel=1e-12)
    assert slope["se"] == pytest.approx(expected_se, rel=1e-12)
    assert slope["pval"] == pytest.approx(expected_p, rel=1e-12)
    assert slope["Q"] == pytest.approx(expected_Q, rel=1e-12)
    assert slope["Q_pval"] == pytest.approx(expected_Qp, rel=1e-12)

    assert intercept["b"] == pytest.approx(expected_bi, rel=1e-12)
    assert intercept["se"] == pytest.approx(expected_sei, rel=1e-12)
    assert intercept["pval"] == pytest.approx(expected_pi, rel=1e-12)


def test_mr_sign_basic() -> None:
    from scipy.stats import binomtest
    from genal.MR import mr_sign

    beta_e = np.array([1, 1, -1, -1, 1, -1], dtype=float)
    beta_o = np.array([1, -1, -1, 1, 1, -1], dtype=float)

    res = mr_sign(beta_e, beta_o)[0]
    assert res["nSNP"] == 6
    assert res["b"] == pytest.approx((4 / 6 - 0.5) * 2)
    assert res["pval"] == pytest.approx(binomtest(4, 6, p=0.5).pvalue)


def test_mode_helpers_validate_method() -> None:
    from genal.MR import bootstrap_mode_iteration

    with pytest.raises(ValueError):
        bootstrap_mode_iteration(
            0,
            BETA_IV=np.array([1.0, 2.0, 3.0]),
            SE_IV=np.array([0.1, 0.1, 0.1]),
            phi=1,
            method="BAD",
        )

