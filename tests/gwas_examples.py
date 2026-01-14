from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


TESTS_DIR = Path(__file__).resolve().parent
TESTS_DATA_DIR = TESTS_DIR / "tests_data"


@dataclass(frozen=True)
class GWASSample:
    name: str
    df: pd.DataFrame


def _ensure_standard_columns(df: pd.DataFrame) -> pd.DataFrame:
    required = {"CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    return df


def _filter_reasonable_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Conservative filters for downstream MR/QA tests.
    Keeps only biallelic A/C/G/T alleles, finite numeric effect columns, and valid P range.
    """
    df = df.copy()
    for col in ["BETA", "SE", "P", "EAF"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

    allele_ok = (
        df["EA"].astype(str).str.fullmatch(r"[ACGT]", na=False)
        & df["NEA"].astype(str).str.fullmatch(r"[ACGT]", na=False)
    )
    numeric_ok = (
        df["CHR"].notna()
        & df["POS"].notna()
        & df["BETA"].notna()
        & df["SE"].notna()
        & np.isfinite(df["BETA"])
        & np.isfinite(df["SE"])
        & (df["SE"] > 0)
        & df["P"].notna()
        & np.isfinite(df["P"])
        & (df["P"] >= 0)
        & (df["P"] <= 1)
        & df["SNP"].notna()
    )
    df = df[allele_ok & numeric_ok].copy()
    df["CHR"] = df["CHR"].astype(int)
    df["POS"] = df["POS"].astype(int)
    return df


def load_gwas_example1_b37(nrows: int = 20_000) -> GWASSample:
    """
    GWAS_example1_b37.txt
    Columns (whitespace-delimited):
      MarkerName Allele1 Allele2 Freq1 Effect StdErr P TotalSampleSize N_effective
    MarkerName format looks like: CHR:POS:SNP
    """
    path = TESTS_DATA_DIR / "GWAS_example1_b37.txt"
    df = pd.read_csv(path, sep=r"\s+", engine="python", nrows=nrows)

    # Parse CHR/POS from MarkerName.
    parts = df["MarkerName"].astype(str).str.split(":", n=2, expand=True)
    df["CHR"] = pd.to_numeric(parts[0], errors="coerce")
    df["POS"] = pd.to_numeric(parts[1], errors="coerce")

    out = pd.DataFrame(
        {
            "CHR": df["CHR"],
            "POS": df["POS"],
            "SNP": df["MarkerName"].astype(str),
            "EA": df["Allele1"].astype(str).str.upper(),
            "NEA": df["Allele2"].astype(str).str.upper(),
            "EAF": pd.to_numeric(df["Freq1"], errors="coerce"),
            "BETA": pd.to_numeric(df["Effect"], errors="coerce"),
            "SE": pd.to_numeric(df["StdErr"], errors="coerce"),
            "P": pd.to_numeric(df["P"], errors="coerce"),
        }
    )
    return GWASSample("GWAS_example1_b37", _ensure_standard_columns(out))


def load_gwas_example2_b38(nrows: int = 50_000) -> GWASSample:
    """
    GWAS_example2_b38.csv (actually whitespace-delimited in the fixture).
    Columns:
      CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ ... BETA SE CHISQ LOG10P ...
    We derive P from LOG10P: P = 10 ** (-LOG10P).
    """
    path = TESTS_DATA_DIR / "GWAS_example2_b38.csv"
    df = pd.read_csv(path, sep=r"\s+", engine="python", nrows=nrows)

    log10p = pd.to_numeric(df.get("LOG10P"), errors="coerce")
    p = np.power(10.0, -log10p)

    out = pd.DataFrame(
        {
            "CHR": pd.to_numeric(df.get("CHROM"), errors="coerce"),
            "POS": pd.to_numeric(df.get("GENPOS"), errors="coerce"),
            "SNP": df.get("ID").astype(str),
            "EA": df.get("ALLELE1").astype(str).str.upper(),
            "NEA": df.get("ALLELE0").astype(str).str.upper(),
            "EAF": pd.to_numeric(df.get("A1FREQ"), errors="coerce"),
            "BETA": pd.to_numeric(df.get("BETA"), errors="coerce"),
            "SE": pd.to_numeric(df.get("SE"), errors="coerce"),
            "P": pd.to_numeric(p, errors="coerce"),
        }
    )
    return GWASSample("GWAS_example2_b38", _ensure_standard_columns(out))


def load_gwas_example3_b38_parquet(
    nrows: int = 50_000, row_group: int = 0
) -> GWASSample:
    """
    GWAS_example3_b38.parquet
    Has already-standardized columns: CHR, POS, NEA, EA, EAF, BETA, P, SNP, SE.
    We read only one parquet row-group and slice to nrows for speed.
    """
    import pyarrow.parquet as pq

    path = TESTS_DATA_DIR / "GWAS_example3_b38.parquet"
    pf = pq.ParquetFile(path)
    cols = ["CHR", "POS", "SNP", "EA", "NEA", "EAF", "BETA", "SE", "P"]
    table = pf.read_row_group(row_group, columns=cols)
    df = table.to_pandas().head(nrows)
    out = df[cols].copy()
    return GWASSample("GWAS_example3_b38", _ensure_standard_columns(out))


def load_gwas_example4_b37(nrows: int = 50_000) -> GWASSample:
    """
    GWAS_example4_b37.tsv (tab-delimited).
    Columns:
      chromosome base_pair_location effect_allele_frequency beta standard_error p_value ... effect_allele other_allele
    No SNP ID column: we generate a deterministic one from CHR/POS/alleles.
    """
    path = TESTS_DATA_DIR / "GWAS_example4_b37.tsv"
    df = pd.read_csv(path, sep="\t", nrows=nrows)

    chr_ = pd.to_numeric(df.get("chromosome"), errors="coerce")
    pos_ = pd.to_numeric(df.get("base_pair_location"), errors="coerce")
    ea = df.get("effect_allele").astype(str).str.upper()
    nea = df.get("other_allele").astype(str).str.upper()

    snp = chr_.astype("Int64").astype(str) + ":" + pos_.astype("Int64").astype(str) + ":" + nea + ":" + ea

    out = pd.DataFrame(
        {
            "CHR": chr_,
            "POS": pos_,
            "SNP": snp,
            "EA": ea,
            "NEA": nea,
            "EAF": pd.to_numeric(df.get("effect_allele_frequency"), errors="coerce"),
            "BETA": pd.to_numeric(df.get("beta"), errors="coerce"),
            "SE": pd.to_numeric(df.get("standard_error"), errors="coerce"),
            "P": pd.to_numeric(df.get("p_value"), errors="coerce"),
        }
    )
    return GWASSample("GWAS_example4_b37", _ensure_standard_columns(out))


def corrupt_gwas_df(
    df: pd.DataFrame,
    *,
    seed: int = 0,
    frac: float = 0.01,
    modes: set[Literal["bad_p", "bad_se", "bad_alleles", "bad_types"]] | None = None,
) -> pd.DataFrame:
    """
    Return a copy with injected errors to exercise validation/sanitization code paths.
    """
    if modes is None:
        modes = {"bad_p", "bad_se", "bad_alleles", "bad_types"}

    rng = np.random.default_rng(seed)
    out = df.copy()
    n = len(out)
    if n == 0:
        return out
    k = max(1, int(n * frac))
    idx = rng.choice(out.index.to_numpy(), size=k, replace=False)

    if "bad_p" in modes and "P" in out.columns:
        out["P"] = out["P"].astype("object")
        out.loc[idx[: k // 2], "P"] = 2  # out of range
        out.loc[idx[k // 2 :], "P"] = "not_a_number"

    if "bad_se" in modes and "SE" in out.columns:
        out.loc[idx, "SE"] = -1

    if "bad_alleles" in modes:
        out.loc[idx, "EA"] = "N"
        out.loc[idx, "NEA"] = "AT"

    if "bad_types" in modes:
        out["CHR"] = out["CHR"].astype("object")
        out["POS"] = out["POS"].astype("object")
        out.loc[idx, "CHR"] = out.loc[idx, "CHR"].astype(str).radd("chr")
        out.loc[idx, "POS"] = out.loc[idx, "POS"].astype(str).radd("pos")

    return out


def make_mr_outcome_from_exposure(
    exposure: pd.DataFrame,
    *,
    causal_beta: float,
    noise_seed: int,
    noise_sd: float,
) -> pd.DataFrame:
    """
    Build a synthetic outcome from a real exposure (same SNPs/alleles) with known causal effect.
    """
    rng = np.random.default_rng(noise_seed)
    out = exposure.copy()
    out["BETA"] = causal_beta * out["BETA"].astype(float) + rng.normal(0.0, noise_sd, size=len(out))
    return out
