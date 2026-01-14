from __future__ import annotations

import pandas as pd

from gwas_examples import _filter_reasonable_rows, corrupt_gwas_df, load_gwas_example4_b37


def test_preprocess_handles_corrupted_rows_by_sanitizing_or_dropping() -> None:
    """
    Realistic path: load real GWAS rows, inject common data issues, then run preprocessing.
    """
    from genal.Geno import Geno

    base = _filter_reasonable_rows(load_gwas_example4_b37(nrows=50_000).df).head(2_000)
    assert base.shape[0] > 100

    corrupted = corrupt_gwas_df(base, seed=123, frac=0.05)
    g = Geno(corrupted, keep_columns=True)

    # Should not error; should delete invalid rows under Fill_delete.
    before = g.data.shape[0]
    g.preprocess_data(preprocessing="Fill_delete", keep_indel=False, keep_dups=False)
    after = g.data.shape[0]

    assert after < before
    # After Fill_delete, standard columns should contain no missing values.
    assert g.data[["CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P"]].isna().sum().sum() == 0


def test_preprocess_warns_when_effect_without_allele() -> None:
    from genal.geno_tools import fill_ea_nea
    from genal.tools import load_reference_panel
    from pathlib import Path
    import pytest

    base = _filter_reasonable_rows(load_gwas_example4_b37(nrows=20_000).df).head(200)
    df = base[["CHR", "POS", "BETA"]].copy()

    # Use a local reference panel (no downloads) and ensure the warning triggers.
    tests_dir = Path(__file__).resolve().parents[1]
    ref_prefix = str(tests_dir / "tests_data" / "StudyA_genetic_files" / "plink_sampled_9")
    ref = load_reference_panel(ref_prefix)

    with pytest.warns(UserWarning):
        out = fill_ea_nea(df, ref)
    assert {"EA", "NEA"}.issubset(out.columns)
