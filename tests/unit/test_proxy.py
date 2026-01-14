from __future__ import annotations

import pandas as pd
import pytest


def _ld_row(**kwargs):
    defaults = {
        "SNP_A": "rs_base",
        "SNP_B": "rs_proxy",
        "R": 0.8,
        "BP_A": 123,
        "CHR_A": 1,
        "MAJ_A": "G",
        "NONMAJ_A": "A",
        "NONMAJ_FREQ_A": 0.2,
        "BP_B": 456,
        "CHR_B": 1,
        "MAJ_B": "C",
        "NONMAJ_B": "T",
        "NONMAJ_FREQ_B": 0.3,
    }
    defaults.update(kwargs)
    return defaults


def test_query_outcome_proxy_rewrites_to_base_snp() -> None:
    from genal.proxy import query_outcome_proxy

    df_outcome = pd.DataFrame(
        {
            "SNP": ["rs_proxy"],
            "CHR": [1],
            "POS": [456],
            "EA": ["T"],  # NONMAJ_B
            "NEA": ["C"],  # MAJ_B
            "BETA": [0.5],
            "SE": [0.1],
            "EAF": [0.3],
        }
    )
    ld = pd.DataFrame([_ld_row()])

    out = query_outcome_proxy(df_outcome, ld, snps_to_extract=[], snps_df=df_outcome.SNP.values)
    assert out.shape[0] == 1
    assert out.loc[0, "SNP"] == "rs_base"  # proxy rewritten to original SNP
    assert out.loc[0, "EA"] == "A"  # NONMAJ_A
    assert out.loc[0, "NEA"] == "G"  # MAJ_A
    assert out.loc[0, "BETA"] == pytest.approx(0.5)


def test_apply_proxies_replaces_snp_ids() -> None:
    from genal.proxy import apply_proxies

    df = pd.DataFrame(
        {
            "SNP": ["rs_base", "rs_other"],
            "EA": ["A", "N"],  # rs_other will be dropped due to allele mismatch
            "BETA": [1.0, 2.0],
        }
    )
    ld = pd.DataFrame(
        [
            _ld_row(SNP_A="rs_base", SNP_B="rs_proxy", R=0.9),
            _ld_row(SNP_A="rs_other", SNP_B="rs_proxy2", R=0.9, MAJ_A="A", NONMAJ_A="G"),
        ]
    )

    out = apply_proxies(df, ld)
    assert out["SNP"].tolist() == ["rs_proxy"]
    assert out.loc[0, "EA"] == "T"  # NONMAJ_B
    assert out.loc[0, "BETA"] == pytest.approx(1.0)


def test_proxy_functions_require_ld_dataframe() -> None:
    from genal.proxy import apply_proxies, query_outcome_proxy

    df = pd.DataFrame({"SNP": ["rs1"], "EA": ["A"], "BETA": [0.1]})
    with pytest.raises(ValueError):
        apply_proxies(df, ld=None)

    with pytest.raises(ValueError):
        query_outcome_proxy(df, ld=None, snps_to_extract=[])

