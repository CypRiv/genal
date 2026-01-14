from __future__ import annotations

import pandas as pd

from gwas_examples import (
    _filter_reasonable_rows,
    load_gwas_example1_b37,
    load_gwas_example2_b38,
    load_gwas_example3_b38_parquet,
    load_gwas_example4_b37,
)


def _assert_usable(df: pd.DataFrame) -> None:
    assert df.shape[0] > 0
    assert {"CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P"}.issubset(df.columns)
    assert df["SNP"].notna().any()


def test_can_load_all_gwas_examples_small_sample() -> None:
    _assert_usable(load_gwas_example1_b37(nrows=20_000).df)
    _assert_usable(load_gwas_example2_b38(nrows=50_000).df)
    _assert_usable(load_gwas_example3_b38_parquet(nrows=50_000).df)
    _assert_usable(load_gwas_example4_b37(nrows=50_000).df)


def test_can_build_geno_from_each_example_after_standardization() -> None:
    from genal.Geno import Geno

    for sample in [
        load_gwas_example1_b37(nrows=20_000),
        load_gwas_example2_b38(nrows=50_000),
        load_gwas_example3_b38_parquet(nrows=20_000),
        load_gwas_example4_b37(nrows=20_000),
    ]:
        df = _filter_reasonable_rows(sample.df).head(1_000)
        _assert_usable(df)
        g = Geno(df)
        assert g.data.shape[0] == df.shape[0]

