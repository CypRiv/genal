from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


def _minimal_gwas_df(snps: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "SNP": snps,
            "BETA": [0.1] * len(snps),
            "SE": [0.01] * len(snps),
            "EA": ["A"] * len(snps),
            "NEA": ["G"] * len(snps),
        }
    )


def test_query_outcome_func_validates_exposure_columns() -> None:
    from genal.MR_tools import query_outcome_func

    exposure_missing = pd.DataFrame({"SNP": ["rs1"], "BETA": [0.1], "SE": [0.01]})
    with pytest.raises(ValueError, match="EA"):
        query_outcome_func(
            exposure_missing,
            outcome="not_used.h5",
            name=None,
            proxy=False,
            reference_panel="EUR_37",
            kb=5000,
            r2=0.8,
            window_snps=1000000,
            cpus=1,
        )


def test_query_outcome_func_validates_outcome_type() -> None:
    from genal.MR_tools import query_outcome_func

    exposure = _minimal_gwas_df(["rs1"])
    with pytest.raises(ValueError, match="Geno object or filepath"):
        query_outcome_func(
            exposure,
            outcome=123,
            name=None,
            proxy=False,
            reference_panel="EUR_37",
            kb=5000,
            r2=0.8,
            window_snps=1000000,
            cpus=1,
        )


def test_load_outcome_from_filepath_validation(tmp_path: Path) -> None:
    from genal.MR_tools import load_outcome_from_filepath

    with pytest.raises(ValueError, match="doesn't lead to a file"):
        load_outcome_from_filepath(str(tmp_path / "missing.h5"))

    wrong_ext = tmp_path / "outcome.txt"
    wrong_ext.write_text("nope")
    with pytest.raises(ValueError, match="needs to be in .h5"):
        load_outcome_from_filepath(str(wrong_ext))


def test_query_outcome_func_no_overlap_returns_empty_geno_outcome() -> None:
    from genal.Geno import Geno
    from genal.MR_tools import query_outcome_func

    exposure = Geno(_minimal_gwas_df(["rs1", "rs2"]))
    outcome = Geno(_minimal_gwas_df(["rsX", "rsY"]))

    exp_df, out_df, out_name = query_outcome_func(
        exposure.data,
        outcome=outcome,
        name=None,
        proxy=False,
        reference_panel="EUR_37",
        kb=5000,
        r2=0.8,
        window_snps=1000000,
        cpus=1,
    )

    assert exp_df.empty
    assert out_df.empty
    assert out_name == outcome.name


def test_query_outcome_func_path_outcome(tmp_path: Path) -> None:
    from genal.MR_tools import query_outcome_func

    pytest.importorskip("tables", reason="pandas HDF5 support requires optional dependency 'tables'")

    exposure = _minimal_gwas_df(["rs1", "rs2"])
    outcome = _minimal_gwas_df(["rs2", "rs3"])

    outcome_path = tmp_path / "outcome.h5"
    outcome.to_hdf(outcome_path, key="data", mode="w")

    exp_df, out_df, out_name = query_outcome_func(
        exposure,
        outcome=str(outcome_path),
        name=None,
        proxy=False,
        reference_panel="EUR_37",
        kb=5000,
        r2=0.8,
        window_snps=1000000,
        cpus=1,
    )

    assert set(exp_df["SNP"]) == {"rs2"}
    assert set(out_df["SNP"]) == {"rs2"}
    assert out_name == "outcome"


def test_geno_mr_after_no_overlap_should_return_empty_not_error() -> None:
    from genal.Geno import Geno

    exposure = Geno(_minimal_gwas_df(["rs1", "rs2"]))
    outcome = Geno(_minimal_gwas_df(["rsX", "rsY"]))

    exposure.query_outcome(outcome, proxy=False)
    res = exposure.MR(methods=["IVW"], action=3, heterogeneity=False)
    assert res.empty
