from __future__ import annotations

from pathlib import Path

import pandas as pd


def test_geno_preprocess_fill_adds_ids_and_alleles(
    reset_genal_config, studyA_genetic_dir: Path
) -> None:
    from genal.Geno import Geno
    from genal.tools import load_reference_panel

    ref_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(ref_prefix)

    subset = ref[(ref["A1"].str.len() == 1) & (ref["A2"].str.len() == 1)].head(25)
    df = subset[["CHR", "POS"]].copy()
    df["BETA"] = 0.1
    df["SE"] = 0.05
    df["P"] = 0.05

    g = Geno(df)
    g.preprocess_data(preprocessing="Fill", reference_panel=ref_prefix, fill_snpids=True)

    assert g.data.shape[0] == subset.shape[0]
    assert g.data["SNP"].notna().all()
    assert g.data["EA"].notna().all()
    assert g.data["NEA"].notna().all()

