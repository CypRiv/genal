from __future__ import annotations

import pandas as pd
import pytest


def test_harmonize_mr_flips_and_filters() -> None:
    from genal.MR_tools import harmonize_MR

    df_exposure = pd.DataFrame(
        {
            "SNP": ["rs_aligned", "rs_inverted", "rs_to_flip", "rs_mismatch", "rs_pal"],
            "EA": ["A", "A", "A", "A", "A"],
            "NEA": ["C", "G", "C", "C", "T"],
            "BETA": [0.10, 0.20, 0.30, 0.40, 0.50],
            "SE": [0.01, 0.01, 0.01, 0.01, 0.01],
            "EAF": [0.10, 0.20, 0.30, 0.40, 0.10],
        }
    )
    df_outcome = pd.DataFrame(
        {
            "SNP": ["rs_aligned", "rs_inverted", "rs_to_flip", "rs_mismatch", "rs_pal"],
            "EA": ["A", "G", "T", "A", "T"],
            "NEA": ["C", "A", "G", "T", "A"],
            "BETA": [1.0, 2.0, 3.0, 4.0, 5.0],
            "SE": [0.02, 0.02, 0.02, 0.02, 0.02],
            "EAF": [0.11, 0.21, 0.31, 0.41, 0.91],
        }
    )

    out = harmonize_MR(df_exposure, df_outcome, action=3)
    assert set(out["SNP"]) == {"rs_aligned", "rs_inverted", "rs_to_flip"}

    inv = out[out["SNP"] == "rs_inverted"].iloc[0]
    assert inv["BETA_o"] == pytest.approx(-2.0)

    flipped = out[out["SNP"] == "rs_to_flip"].iloc[0]
    assert flipped["BETA_o"] == pytest.approx(3.0)


def test_apply_action_2_removes_ambiguous_and_flips() -> None:
    from genal.MR_tools import apply_action_2

    df = pd.DataFrame(
        {
            "SNP": ["rs_ambiguous", "rs_flip", "rs_ok"],
            "palindrome": [True, True, False],
            "EAF_e": [0.5, 0.1, 0.2],
            "EAF_o": [0.1, 0.9, 0.2],
            "BETA_o": [1.0, 2.0, 3.0],
        }
    )
    out = apply_action_2(df, eaf_threshold=0.42)
    assert "rs_ambiguous" not in set(out["SNP"])

    flip = out[out["SNP"] == "rs_flip"].iloc[0]
    assert flip["BETA_o"] == pytest.approx(-2.0)
    assert flip["EAF_o"] == pytest.approx(0.1)


def test_flip_alleles_complements() -> None:
    from genal.MR_tools import flip_alleles

    s = pd.Series(["A", "C", "G", "T"])
    assert flip_alleles(s).tolist() == ["T", "G", "C", "A"]


def test_check_required_columns_errors() -> None:
    from genal.MR_tools import check_required_columns

    df = pd.DataFrame({"SNP": ["rs1"]})
    with pytest.raises(ValueError, match="BETA"):
        check_required_columns(df, ["SNP", "BETA"])

