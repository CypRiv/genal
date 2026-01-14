from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


def test_set_phenotype_infers_binary_and_recodes(phenotype_csv_path: Path) -> None:
    from genal.association import set_phenotype_func

    df = pd.read_csv(phenotype_csv_path, nrows=200)
    out, pheno_type = set_phenotype_func(
        df, PHENO="binary_column", PHENO_type=None, IID="IID", FID="FID"
    )

    assert pheno_type == "binary"
    assert {"IID", "FID", "PHENO"}.issubset(out.columns)
    assert set(out["PHENO"].dropna().unique()).issubset({0, 1})


def test_set_phenotype_infers_quantitative(phenotype_csv_path: Path) -> None:
    from genal.association import set_phenotype_func

    df = pd.read_csv(phenotype_csv_path, nrows=200)
    out, pheno_type = set_phenotype_func(
        df, PHENO="continuous_column", PHENO_type=None, IID="IID", FID="FID"
    )
    assert pheno_type == "quant"
    assert out["PHENO"].dtype.kind in {"f", "i", "u"}


def test_set_phenotype_binary_requires_two_levels() -> None:
    from genal.association import set_phenotype_func

    df = pd.DataFrame({"IID": [1, 2, 3], "PH": [0, 1, 2]})
    with pytest.raises(ValueError):
        set_phenotype_func(df, PHENO="PH", PHENO_type="binary", IID="IID")


def test_set_phenotype_quant_requires_numeric() -> None:
    from genal.association import set_phenotype_func

    df = pd.DataFrame({"IID": [1, 2, 3], "PH": ["a", "b", "c"]})
    with pytest.raises(ValueError):
        set_phenotype_func(df, PHENO="PH", PHENO_type="quant", IID="IID")


def test_set_phenotype_alternate_control_flips_coding() -> None:
    from genal.association import set_phenotype_func

    df = pd.DataFrame({"IID": [1, 2, 3, 4, 5], "PH": [3, 3, 4, 4, 4]})
    out_default, _ = set_phenotype_func(df, PHENO="PH", PHENO_type="binary", IID="IID")
    out_alt, _ = set_phenotype_func(
        df, PHENO="PH", PHENO_type="binary", IID="IID", alternate_control=True
    )

    assert out_default["PHENO"].sum() != out_alt["PHENO"].sum()

