from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gwas_examples import (
    _filter_reasonable_rows,
    load_gwas_example1_b37,
    make_mr_outcome_from_exposure,
)


@pytest.mark.slow
def test_stochastic_mr_methods_close_to_ground_truth() -> None:
    """
    Egger-boot and mode methods involve bootstrap sampling (and multiprocessing in the
    current implementation). We can't assert exact equality, but we can assert they are
    reasonably close to the known causal effect used to generate the outcome.
    """
    baseline_path = Path(__file__).resolve().parents[1] / "baselines" / "mr_baseline_v1.json"
    baseline = json.loads(baseline_path.read_text())
    scenario = baseline["scenario"]

    causal_beta = float(scenario["causal_beta"])

    # Reconstruct instruments.
    sample = load_gwas_example1_b37(nrows=baseline["source"]["nrows_read"]).df
    filtered = _filter_reasonable_rows(sample).drop_duplicates(subset=["SNP"], keep="first")
    instrument_snps = baseline["source"]["instrument_snps"]
    inst = filtered[filtered["SNP"].isin(instrument_snps)].copy()
    inst = inst.set_index("SNP").loc[instrument_snps].reset_index()

    exposure_df = inst[["SNP", "BETA", "SE", "EA", "NEA", "EAF"]].copy()
    outcome_df = make_mr_outcome_from_exposure(
        exposure_df,
        causal_beta=causal_beta,
        noise_seed=int(scenario["noise_seed"]),
        noise_sd=float(scenario["noise_sd"]),
    )

    from genal.Geno import Geno

    exposure = Geno(exposure_df)
    outcome = Geno(outcome_df)

    # Reduce worker count to keep runtime manageable in CI / laptops.
    exposure.cpus = 1

    exposure.query_outcome(outcome, proxy=False)
    res = exposure.MR(
        methods=["Egger-boot", "Simple-mode", "Weighted-mode"],
        action=int(scenario["action"]),
        nboot=50,
        phi=float(scenario["phi"]),
        heterogeneity=False,
        exposure_name="exposure",
        outcome_name="outcome",
    )

    methods = set(res["method"].tolist())
    assert "MR Egger bootstrap" in methods
    assert "Egger Intercept bootstrap" in methods
    assert "Simple mode" in methods
    assert "Weighted mode" in methods

    slope = float(res[res["method"] == "MR Egger bootstrap"]["b"].iloc[0])
    intercept = float(res[res["method"] == "Egger Intercept bootstrap"]["b"].iloc[0])
    assert slope == pytest.approx(causal_beta, abs=1.0)
    assert intercept == pytest.approx(0.0, abs=0.5)

    for m in ["Simple mode", "Weighted mode"]:
        b = float(res[res["method"] == m]["b"].iloc[0])
        assert b == pytest.approx(causal_beta, abs=1.0)
