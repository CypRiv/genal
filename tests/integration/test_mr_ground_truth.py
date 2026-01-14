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


def test_mr_baseline_matches_ground_truth() -> None:
    """
    Compares core MR outputs against a persisted baseline built from real GWAS example data.
    """
    baseline_path = Path(__file__).resolve().parents[1] / "baselines" / "mr_baseline_v1.json"
    baseline = json.loads(baseline_path.read_text())

    scenario = baseline["scenario"]
    methods = scenario["methods"]

    # Reconstruct the exact instrument set from the real GWAS fixture.
    sample = load_gwas_example1_b37(nrows=baseline["source"]["nrows_read"]).df
    filtered = _filter_reasonable_rows(sample).drop_duplicates(subset=["SNP"], keep="first")

    instrument_snps = baseline["source"]["instrument_snps"]
    inst = filtered[filtered["SNP"].isin(instrument_snps)].copy()
    assert inst.shape[0] == baseline["source"]["n_instruments"]

    inst = inst.set_index("SNP").loc[instrument_snps].reset_index()
    exposure_df = inst[["SNP", "BETA", "SE", "EA", "NEA", "EAF"]].copy()
    outcome_df = make_mr_outcome_from_exposure(
        exposure_df,
        causal_beta=float(scenario["causal_beta"]),
        noise_seed=int(scenario["noise_seed"]),
        noise_sd=float(scenario["noise_sd"]),
    )

    from genal.Geno import Geno

    # Make bootstrap-based methods deterministic.
    np.random.seed(int(scenario["bootstrap_seed"]))

    exposure = Geno(exposure_df)
    outcome = Geno(outcome_df)
    exposure.query_outcome(outcome, proxy=False)
    res = exposure.MR(
        methods=methods,
        action=int(scenario["action"]),
        nboot=int(scenario["nboot"]),
        penk=float(scenario["penk"]),
        phi=float(scenario["phi"]),
        heterogeneity=False,
        exposure_name="exposure",
        outcome_name="outcome",
    )

    expected = pd.DataFrame(baseline["expected"]).sort_values("method").reset_index(drop=True)
    got = res.sort_values("method").reset_index(drop=True)

    assert got["method"].tolist() == expected["method"].tolist()
    assert got["nSNP"].tolist() == expected["nSNP"].tolist()

    # Compare effect estimates and SEs with tolerances (SE can vary slightly between versions).
    for col in ["b", "se"]:
        for i in range(len(expected)):
            exp = expected.loc[i, col]
            gotv = got.loc[i, col]
            if pd.isna(exp):
                assert pd.isna(gotv)
            else:
                assert gotv == pytest.approx(float(exp), rel=1e-6, abs=1e-6)
