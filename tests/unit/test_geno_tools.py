from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def test_remove_na_drops_rows_in_standard_columns() -> None:
    from genal.geno_tools import remove_na

    df = pd.DataFrame(
        {
            "CHR": [1, 1, 1],
            "POS": [10, 20, 30],
            "SNP": ["rs1", "rs2", None],
            "BETA": [0.1, None, 0.3],
        }
    )
    remove_na(df)
    assert df.shape[0] == 1
    assert df.iloc[0]["SNP"] == "rs1"


def test_check_snp_column_deduplicates() -> None:
    from genal.geno_tools import check_snp_column

    df = pd.DataFrame({"SNP": ["rs1", "rs1", "rs2"], "BETA": [0.1, 0.2, 0.3]})
    check_snp_column(df)
    assert df["SNP"].tolist() == ["rs1", "rs2"]


def test_check_allele_column_validates_and_handles_indels() -> None:
    from genal.geno_tools import check_allele_column

    df = pd.DataFrame({"EA": ["a", "t", "N", "AT", 1, None]})
    check_allele_column(df, "EA", keep_indel=False)
    assert df["EA"].tolist()[:2] == ["A", "T"]
    assert pd.isna(df.loc[2, "EA"])  # invalid nucleotide
    assert pd.isna(df.loc[3, "EA"])  # indel removed


def test_fill_se_p_creates_missing_column() -> None:
    from genal.geno_tools import fill_se_p

    df = pd.DataFrame({"BETA": [0.2, -0.2], "P": [0.05, 0.05]})
    fill_se_p(df)
    assert "SE" in df.columns
    assert df["SE"].notna().all()

    df2 = pd.DataFrame({"BETA": [0.2, -0.2], "SE": [0.1, 0.1]})
    fill_se_p(df2)
    assert "P" in df2.columns
    assert (df2["P"] > 0).all() and (df2["P"] <= 1).all()


def test_fill_fstatistic_prefers_beta_se_and_falls_back_to_p() -> None:
    from genal.geno_tools import fill_fstatistic

    df = pd.DataFrame(
        {
            "BETA": [0.2, 0.2, np.nan],
            "SE": [0.1, 0.0, np.nan],
            "P": [np.nan, np.nan, 0.05],
        }
    )
    fill_fstatistic(df)
    assert "FSTAT" in df.columns
    assert df.loc[0, "FSTAT"] == pytest.approx((0.2 / 0.1) ** 2)
    assert np.isinf(df.loc[1, "FSTAT"])
    assert df.loc[2, "FSTAT"] > 0


def test_check_p_column_sanitizes_range() -> None:
    from genal.geno_tools import check_p_column

    df = pd.DataFrame({"P": ["0.1", "2", "-1", "nan"]})
    check_p_column(df)
    assert df["P"].tolist()[0] == pytest.approx(0.1)
    assert pd.isna(df.loc[1, "P"])
    assert pd.isna(df.loc[2, "P"])
    assert pd.isna(df.loc[3, "P"])


def test_check_beta_column_detects_or_and_transforms() -> None:
    from genal.geno_tools import check_beta_column

    df_or = pd.DataFrame({"BETA": [1.1, 0.9, 1.0], "SE": [0.1, 0.1, 0.1]})
    check_beta_column(df_or, effect_column=None)
    assert "SE" not in df_or.columns  # dropped for OR inputs
    assert df_or["BETA"].tolist()[0] == pytest.approx(np.log(1.1))

    df_beta = pd.DataFrame({"BETA": [0.1, -0.2, 0.3]})
    check_beta_column(df_beta, effect_column=None)
    assert df_beta["BETA"].tolist()[1] == pytest.approx(-0.2)

    with pytest.raises(ValueError):
        check_beta_column(pd.DataFrame({"BETA": [1.0]}), effect_column="BAD")


def test_fill_ea_nea_and_fill_nea(studyA_genetic_dir: Path) -> None:
    from genal.tools import load_reference_panel
    from genal.geno_tools import fill_ea_nea, fill_nea

    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9"))
    row = ref.iloc[0][["CHR", "POS"]].to_dict()

    df_missing = pd.DataFrame({**row, "BETA": [0.1]})
    with pytest.warns(UserWarning):
        out = fill_ea_nea(df_missing, ref)
    assert {"EA", "NEA"}.issubset(out.columns)

    df_nea = pd.DataFrame({**row, "EA": [ref.iloc[0]["A1"]]})
    out2 = fill_nea(df_nea, ref)
    assert out2.loc[0, "NEA"] == ref.iloc[0]["A2"]


def test_fill_coordinates_and_snpids(studyA_genetic_dir: Path) -> None:
    from genal.tools import load_reference_panel
    from genal.geno_tools import fill_coordinates_func, fill_snpids_func

    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9"))
    # Use a biallelic SNP (no indel) for keep_indel=False behavior.
    ref_row = ref[(ref["A1"].str.len() == 1) & (ref["A2"].str.len() == 1)].iloc[0]

    df = pd.DataFrame({"SNP": [ref_row["SNP"]]})
    out = fill_coordinates_func(df, ref)
    assert out.loc[0, "CHR"] == ref_row["CHR"]
    assert out.loc[0, "POS"] == ref_row["POS"]

    df2 = pd.DataFrame(
        {
            "CHR": [ref_row["CHR"]],
            "POS": [ref_row["POS"]],
            "EA": [ref_row["A1"]],
            "NEA": [ref_row["A2"]],
        }
    )
    out2 = fill_snpids_func(df2, ref, keep_indel=False)
    assert out2.loc[0, "SNP"] == ref_row["SNP"]

    df3 = pd.DataFrame({"CHR": [999], "POS": [999], "EA": ["A"], "NEA": ["C"]})
    with pytest.warns(UserWarning):
        out3 = fill_snpids_func(df3, ref, keep_indel=False)
    assert out3.loc[0, "SNP"] == "999:999:C:A"


def test_check_int_column_extracts_digits() -> None:
    from genal.geno_tools import check_int_column

    df = pd.DataFrame({"CHR": ["1", "chr2", "X", None, "3.7"]})
    check_int_column(df, "CHR")
    assert df.loc[0, "CHR"] == 1
    assert df.loc[1, "CHR"] == 2
    assert pd.isna(df.loc[2, "CHR"])
    assert pd.isna(df.loc[3, "CHR"])
    assert df.loc[4, "CHR"] == 3  # digits are extracted from strings


def test_adjust_column_names_validates_duplicates() -> None:
    from genal.geno_tools import adjust_column_names

    df = pd.DataFrame({"chrom": [1], "pos": [10], "CHR": [1]})
    with pytest.raises(ValueError):
        adjust_column_names(df, CHR="chrom", POS="pos", SNP="SNP", EA="EA", NEA="NEA", BETA="BETA", SE="SE", P="P", EAF="EAF", keep_columns=True, FSTAT="FSTAT")


def test_check_arguments_defaults_and_validation() -> None:
    from genal.geno_tools import check_arguments

    keep_indel, keep_dups, fill_snpids, fill_coordinates = check_arguments(
        preprocessing="None",
        effect_column=None,
        fill_snpids=None,
        fill_coordinates=None,
        keep_indel=None,
        keep_dups=None,
    )
    assert keep_indel is True
    assert keep_dups is True
    assert fill_snpids is False
    assert fill_coordinates is False

    with pytest.raises(ValueError):
        check_arguments("BAD", None, None, None, None, None)


def test_save_data_writes_csv(tmp_path: Path) -> None:
    from genal.geno_tools import save_data

    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    save_data(df, name="out", path=str(tmp_path), fmt="csv")
    assert (tmp_path / "out.csv").is_file()


def test_update_eaf_requires_plink_path(reset_genal_config, studyA_genetic_dir: Path) -> None:
    from genal.geno_tools import update_eaf_func

    df = pd.DataFrame({"CHR": [1], "POS": [1], "EA": ["A"], "NEA": ["C"]})
    with pytest.raises(ValueError):
        update_eaf_func(
            df,
            reference_panel=str(studyA_genetic_dir / "plink_sampled_9"),
            object_name="obj",
        )
