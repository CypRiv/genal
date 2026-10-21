import pandas as pd
import numpy as np
import scipy.stats as st
import os
import warnings
from collections import Counter

from .constants import STANDARD_COLUMNS, CHR_CODES, NULL_SNP_IDS
from .tools import get_reference_panel_path, run_plink_command, create_tmp, get_plink_path


def remove_na(data):
    """Identify the standard columns containing NA values. Delete rows with NA values."""
    nrows = data.shape[0]
    present_standard_columns = [col for col in STANDARD_COLUMNS if col in data.columns]
    columns_na = (
        data[present_standard_columns]
        .columns[data[present_standard_columns].isna().any()]
        .tolist()
    )
    data.dropna(subset=present_standard_columns, inplace=True)
    n_del = nrows - data.shape[0]
    if n_del > 0:
        print(
            f"Deleted {n_del}({n_del/nrows*100:.3f}%) rows containing NA values in columns {columns_na}. Use preprocessing = 'Fill' to keep the rows containing NA values."
        )
    return


def missing_snp_ids(snp):
    """Boolean mask of missing SNP IDs: NaN or placeholder values such as '', '.', 'NA' (case-insensitive, surrounding spaces ignored)."""
    return snp.isna() | snp.astype(str).str.strip().str.upper().isin(NULL_SNP_IDS)


def check_snp_ids(data):
    """Set placeholder SNP IDs (e.g. '', '.', 'NA', see :func:`missing_snp_ids`) to NaN so that they are handled as missing IDs."""
    placeholder = missing_snp_ids(data["SNP"]) & data["SNP"].notna()
    n_placeholder = placeholder.sum()
    if n_placeholder > 0:
        values = ", ".join(repr(v) for v in data.loc[placeholder, "SNP"].astype(str).unique()[:5])
        data.loc[placeholder, "SNP"] = np.nan
        print(
            f"{n_placeholder}({n_placeholder/data.shape[0]*100:.3f}%) SNP IDs are placeholders ({values}) and have been set to nan (missing ID)."
        )
    return


def check_eaf_column(data):
    """Verify that the EAF column contains numeric values in the range [0, 1]. Set inappropriate values to NaN."""
    nrows = data.shape[0]
    eaf = pd.to_numeric(data["EAF"], errors="coerce")
    n_non_numeric = (eaf.isna() & data["EAF"].notna()).sum()
    out_of_range = (eaf < 0) | (eaf > 1)
    n_out_of_range = out_of_range.sum()
    data["EAF"] = eaf.where(~out_of_range)
    if n_non_numeric + n_out_of_range > 0:
        print(
            f"{n_non_numeric + n_out_of_range}({(n_non_numeric + n_out_of_range)/nrows*100:.3f}%) values in the EAF column have been set to nan: "
            f"{n_out_of_range} outside [0, 1] and {n_non_numeric} non-numeric."
        )
    return


def check_snp_column(data):
    """
    Remove duplicates in the SNP column (in place), keeping the first row. Rows with a missing SNP ID (NaN or placeholder) are never considered duplicates.
    Rows are removed by position, so repeated index labels do not remove other rows; the kept rows keep their index labels.
    """
    duplicated = (~missing_snp_ids(data["SNP"]) & data.duplicated(subset=["SNP"], keep="first")).to_numpy()
    n_del = int(duplicated.sum())
    n_initial = data.shape[0]
    if n_del > 0:
        labels = data.index[~duplicated]
        data.reset_index(drop=True, inplace=True)
        data.drop(index=np.flatnonzero(duplicated), inplace=True)
        data.index = labels
        print(
            f"{n_del}({n_del/n_initial*100:.3f}%) SNPs with duplicated IDs (SNP column) have been removed. Use keep_dups=True to keep them."
        )
    return


def check_allele_column(data, allele_col, keep_indel):
    """
    Verify that the corresponding allele column is upper case strings. Set to nan if not formed with A, T, C, G letters.
    Set to nan if values are insertions/deletions unless keep_indel=True.
    """
    nrows = data.shape[0]
    data[allele_col] = data[allele_col].astype(str).str.upper()
    atcg_condition = data[allele_col].str.contains("^[ATCG]+$", na=False)
    atcg_count = nrows - atcg_condition.sum()
    if atcg_count > 0:
        data.loc[~atcg_condition, allele_col] = np.nan
        print(
            f"{atcg_count}({atcg_count/nrows*100:.3f}%) rows contain non A, T, C, G values in the {allele_col} column and are set to NA."
        )
    if not keep_indel:
        nrows = data.shape[0]
        indel_condition = data[allele_col].str.len() > 1
        indel_count = indel_condition.sum()
        if indel_count > 0:
            data.loc[indel_condition, allele_col] = np.nan
            print(
                f"{indel_count}({indel_count/nrows*100:.3f}%) rows containing insertions/deletions in the {allele_col} column are set to NA. Use keep_indel=True to keep them."
            )
    return


def fill_se_p(data):
    """If either P or SE is missing but the other and BETA are present, fill it."""
    # Ensure SE is numeric and non-negative
    if ("SE" in data.columns):
        data["SE"] = pd.to_numeric(data["SE"], errors="coerce")
        data.loc[data["SE"] < 0, "SE"] = np.nan
        n_missing = data["SE"].isna().sum()
        if n_missing > 0:
            print(
                f"{n_missing}({n_missing/data.shape[0]*100:.3f}%) values in the SE column have been set to nan for being missing, negative or non-numeric."
            )
    # If SE is missing: SE = |BETA| / z with z = isf(P/2). SE is only defined when 0 < SE < inf.
    if ("P" in data.columns) & ("BETA" in data.columns) & ("SE" not in data.columns):
        beta = pd.to_numeric(data["BETA"], errors="coerce")
        p = pd.to_numeric(data["P"], errors="coerce")
        z = pd.Series(st.norm.isf(p / 2), index=data.index)
        with np.errstate(divide="ignore", invalid="ignore"):
            se = beta.abs() / z
        valid = np.isfinite(se) & (se > 0)
        data["SE"] = se.where(valid)
        print("The SE (Standard Error) column has been created.")
        n_invalid = (~valid).sum()
        if n_invalid > 0:
            missing = ~(np.isfinite(beta) & (p >= 0) & (p <= 1))
            beta_zero = ~missing & (beta == 0)
            p_zero = ~missing & ~beta_zero & np.isinf(z)
            p_one = ~missing & ~beta_zero & (z == 0)
            print(
                f"{n_invalid}({n_invalid/data.shape[0]*100:.3f}%) SE values could not be computed and are set to nan: "
                f"{beta_zero.sum()} with BETA = 0, {p_zero.sum()} with P = 0 (or too small to be represented), "
                f"{p_one.sum()} with P = 1, {missing.sum()} with missing or invalid BETA or P."
            )
    # If P is missing: P = 2 * sf(|BETA / SE|), only defined when SE > 0.
    if ("SE" in data.columns) & ("BETA" in data.columns) & ("P" not in data.columns):
        beta = pd.to_numeric(data["BETA"], errors="coerce")
        valid = np.isfinite(beta) & (data["SE"] > 0)
        z = (beta / data["SE"]).where(valid)
        data["P"] = 2 * st.norm.sf(np.abs(z))
        print("The P (P-value) column has been created.")
        n_invalid = (~valid).sum()
        if n_invalid > 0:
            print(
                f"{n_invalid}({n_invalid/data.shape[0]*100:.3f}%) P values could not be computed (missing BETA, or missing or zero SE) and are set to nan."
            )
    return


def fill_fstatistic(data, overwrite=False):
    """
    Compute or fill the per-variant F-statistic (FSTAT) column.
    
    The F-statistic is computed as:
    - Primary: FSTAT = (BETA / SE)² when BETA and SE are available and SE > 0.
    - Fallback: FSTAT = χ²_isf(P, df=1) when the primary route is not computable (missing BETA,
      missing SE or SE = 0) but P is present. For P=0 (extremely significant variants), this produces inf.
    - Variants with SE = 0, a non-zero BETA and no usable P (extremely significant variants) get FSTAT = inf.
    
    Args:
        data (pd.DataFrame): SNP-level DataFrame.
        overwrite (bool): If False (default), only fill missing FSTAT values if the column 
            exists; if it doesn't exist, create it. If True, recompute FSTAT for all rows 
            where computable, overwriting existing values; non-computable rows become NaN.
    
    Returns:
        None: Modifies data in place.
    
    Note:
        FSTAT is NOT added to STANDARD_COLUMNS to avoid row deletion due to missing values.
    """
    nrows = data.shape[0]
    column_created = False
    
    # Determine which rows need computation
    if "FSTAT" not in data.columns:
        data["FSTAT"] = np.nan
        column_created = True
        rows_to_compute = pd.Series([True] * nrows, index=data.index)
    elif overwrite:
        # Clear FSTAT first so non-computable rows become NaN (not stale values)
        data["FSTAT"] = np.nan
        rows_to_compute = pd.Series([True] * nrows, index=data.index)
    else:
        # Only compute for rows with missing FSTAT
        rows_to_compute = data["FSTAT"].isna()
    
    if not rows_to_compute.any():
        return
    
    # Track how many values are assigned
    n_assigned = 0
    
    # Primary route: FSTAT = (BETA / SE)² when BETA and SE are available and SE > 0
    beta_se_computable = pd.Series([False] * nrows, index=data.index)
    if "BETA" in data.columns and "SE" in data.columns:
        beta_se_computable = (
            rows_to_compute &
            data["BETA"].notna() &
            data["SE"].notna() &
            (data["SE"] > 0)
        )
        if beta_se_computable.any():
            data.loc[beta_se_computable, "FSTAT"] = (
                data.loc[beta_se_computable, "BETA"] / data.loc[beta_se_computable, "SE"]
            ) ** 2
            n_assigned += beta_se_computable.sum()
    
    # Fallback route: FSTAT = χ²_isf(P, df=1) for remaining rows where P is present
    # Allow P=0 (produces inf for extremely significant variants)
    p_fallback_computable = pd.Series([False] * nrows, index=data.index)
    if "P" in data.columns:
        p_fallback_computable = (
            rows_to_compute &
            ~beta_se_computable &
            data["P"].notna() &
            (data["P"] >= 0) &
            (data["P"] <= 1)
        )
        if p_fallback_computable.any():
            data.loc[p_fallback_computable, "FSTAT"] = st.chi2.isf(
                data.loc[p_fallback_computable, "P"], df=1
            )
            n_assigned += p_fallback_computable.sum()

    # SE = 0 with a non-zero BETA and no usable P: extremely significant variants, FSTAT = inf
    if "BETA" in data.columns and "SE" in data.columns:
        se_zero = (
            rows_to_compute &
            ~beta_se_computable &
            ~p_fallback_computable &
            (data["SE"] == 0) &
            data["BETA"].notna() &
            (data["BETA"] != 0)
        )
        if se_zero.any():
            data.loc[se_zero, "FSTAT"] = np.inf
            n_assigned += se_zero.sum()

    # Logging
    if column_created:
        print(
            f"The FSTAT (F-statistic) column has been created. "
            f"{n_assigned}({n_assigned/nrows*100:.3f}%) values computed."
        )
    elif overwrite:
        print(
            f"The FSTAT (F-statistic) column has been re-created. "
            f"{n_assigned}({n_assigned/nrows*100:.3f}%) values computed."
        )
    else:
        if n_assigned > 0:
            print(
                f"The FSTAT (F-statistic) column: {n_assigned}({n_assigned/nrows*100:.3f}%) "
                f"missing values have been filled."
            )
    
    return


def check_p_column(data):
    """Verify that the P column contains numeric values in the range [0,1]. Set inappropriate values to NA."""
    nrows = data.shape[0]
    data["P"] = pd.to_numeric(data["P"], errors="coerce")
    data.loc[(data["P"] < 0) | (data["P"] > 1), "P"] = np.nan
    n_missing = data["P"].isna().sum()
    if n_missing > 0:
        print(
            f"{n_missing}({n_missing/nrows*100:.3f}%) values in the P column have been set to nan for being missing, non numeric or out of range [0,1]."
        )
    return


def _scale_check_rows(data):
    """
    Rows usable to compare the effect estimates with the P-values: positive BETA different from 1, finite SE > 0
    and 0 < P < 1. Returns the BETA, SE and z-statistic implied by P (z_P = isf(P/2)) of these rows as arrays.
    """
    beta = pd.to_numeric(data["BETA"], errors="coerce")
    se = pd.to_numeric(data["SE"], errors="coerce")
    p = pd.to_numeric(data["P"], errors="coerce")
    z_p = pd.Series(st.norm.isf(p / 2), index=data.index)
    valid = (
        np.isfinite(beta) & (beta > 0) & (beta != 1) & np.isfinite(se) & (se > 0)
        & (p > 0) & (p < 1) & np.isfinite(z_p) & (z_p > 0)
    )
    return beta[valid].to_numpy(), se[valid].to_numpy(), z_p[valid].to_numpy()


def _z_difference(z_p, z_candidate):
    """Median relative difference between the z-statistics implied by P and candidate z-statistics (exp(median |log ratio|) - 1)."""
    return float(np.exp(np.median(np.abs(np.log(z_candidate / z_p)))) - 1)


def check_beta_column(data, effect_column):
    """
    Make sure the BETA column contains Beta estimates: odds ratios are log-transformed (odds ratios <= 0 give nan).

    effect_column="BETA": BETA is used as is.
    effect_column="OR": BETA is log-transformed. If SE and P are available, SE is kept if |log(OR)|/SE reproduces the
        z-statistic implied by P (SE on the log-odds scale), otherwise it is removed (it is then recomputed from P).
        Without P, SE is kept (assumed to be the standard error of log(OR)). Without SE, it is computed from P later.
    effect_column=None: the scale is determined from the data:
        - BETA contains negative or zero values: Beta estimates.
        - all BETA > 0 with SE and P available: the scale for which |BETA|/SE (Beta estimates) or |log(BETA)|/SE
          (odds ratios) is closer to the z-statistic implied by P (median relative difference). For odds ratios,
          SE is then checked as with effect_column="OR" (it is kept when it reproduces P on the log-odds scale).
        - all BETA > 0 without SE or P: odds ratios if the median of BETA is between 0.5 and 1.5 (heuristic).
          SE, if present, is kept (assumed to be the standard error of log(OR)).
    """
    # Ensure the BETA column is numeric
    data["BETA"] = pd.to_numeric(data["BETA"], errors="coerce")
    force_message = "Use effect_column='BETA' or effect_column='OR' to force the choice."
    has_se_p = "SE" in data.columns and "P" in data.columns

    if effect_column is None:
        n_nonpositive = (data["BETA"] <= 0).sum()
        beta, se, z_p = _scale_check_rows(data) if has_se_p else (np.array([]), None, None)
        if n_nonpositive > 0:
            effect_column = "BETA"
            print(
                f"The BETA column contains {n_nonpositive} negative or zero values (odds ratios are positive): it is treated as Beta estimates. "
                f"Use effect_column='OR' if it is a column of Odds Ratios."
            )
        elif len(beta) > 0:
            # Which scale reproduces the reported P-values?
            diff_beta = _z_difference(z_p, beta / se)
            diff_or = _z_difference(z_p, np.abs(np.log(beta)) / se)
            effect_column = "OR" if diff_or < diff_beta else "BETA"
            print(
                f"All BETA values are positive. Comparison with the z-statistics implied by the P column ({len(beta)} variants): "
                f"median relative difference {diff_beta*100:.2f}% for |BETA|/SE (Beta estimates) and {diff_or*100:.2f}% for |log(BETA)|/SE (Odds Ratios). "
                f"The BETA column is treated as {'Odds Ratios' if effect_column == 'OR' else 'Beta estimates'}. {force_message}"
            )
            if min(diff_beta, diff_or) > 0.1:
                print(
                    "Warning: neither scale reproduces the P-values closely. Check the BETA, SE and P columns."
                )
        else:
            median = data["BETA"].median()
            effect_column = "OR" if 0.5 < median < 1.5 else "BETA"
            print(
                f"All BETA values are positive and the scale cannot be checked against P (SE and P are not both available): "
                f"heuristic based on the median of BETA ({median:.3f}, odds ratios if between 0.5 and 1.5). "
                f"The BETA column is treated as {'Odds Ratios' if effect_column == 'OR' else 'Beta estimates'}. {force_message}"
            )

    ## Log transform the effect column if appropriate
    if effect_column not in ["BETA", "OR"]:
        raise ValueError(
            "The argument effect_column accepts only 'BETA' or 'OR' as values."
        )
    if effect_column == "OR":
        # Check the scale of SE with the odds ratios before transforming them
        se_message = None
        if "SE" in data.columns:
            odds_ratio, se, z_p = _scale_check_rows(data) if has_se_p else (np.array([]), None, None)
            if len(odds_ratio) > 0:
                diff_log = _z_difference(z_p, np.abs(np.log(odds_ratio)) / se)
                diff_or_scale = _z_difference(z_p, np.abs(np.log(odds_ratio)) / (se / odds_ratio))
                if diff_log <= 0.1 and diff_log <= diff_or_scale:
                    se_message = (
                        f"The SE column is kept: |log(OR)|/SE reproduces the z-statistics implied by P (median relative difference {diff_log*100:.2f}%, {len(odds_ratio)} variants)."
                    )
                else:
                    data.drop(columns=["SE"], inplace=True)
                    se_message = (
                        f"The SE column is not on the log-odds scale (median relative difference between |log(OR)|/SE and the z-statistics implied by P: {diff_log*100:.2f}%"
                        f"{f', {diff_or_scale*100:.2f}% if SE is the standard error of the odds ratio' if diff_or_scale < diff_log else ''}). "
                        f"It has been removed and will be computed from BETA and P."
                    )
            else:
                se_message = (
                    "The SE column is kept and assumed to be the standard error of log(OR) (the standard GWAS convention)"
                    f"{'' if has_se_p else ': there is no P column to check it'}."
                )
        n_nonpositive = (data["BETA"] <= 0).sum()
        data["BETA"] = np.log(data["BETA"].where(data["BETA"] > 0))
        print("The BETA column has been log-transformed to obtain Beta estimates.")
        if n_nonpositive > 0:
            print(
                f"{n_nonpositive}({n_nonpositive/data.shape[0]*100:.3f}%) odds ratios were zero or negative and their BETA has been set to nan."
            )
        if se_message:
            print(se_message)
    return


def _reference_at_positions(data, reference_panel_df, columns):
    """Return the reference records (CHR, POS and the requested columns) located at the CHR/POS positions present in data."""
    positions = data[["CHR", "POS"]].dropna().drop_duplicates()
    reference = reference_panel_df[reference_panel_df["POS"].isin(positions["POS"].unique())]
    return reference[["CHR", "POS", *columns]].merge(positions, on=["CHR", "POS"], how="inner")


def _count_rows_in(data, keys, key_values):
    """Count the rows of data whose keys are present in the key_values DataFrame."""
    if key_values.empty:
        return 0
    flagged = data[keys].merge(key_values[keys].drop_duplicates().assign(_flag=True), on=keys, how="left")
    return int(flagged["_flag"].notna().sum())


def fill_ea_nea(data, reference_panel_df):
    """
    Fill in the EA and NEA columns based on reference data (only for data without effect estimates).
    Positions with several reference records (multi-allelic sites) are ambiguous and their alleles are left as nan.

    Raises:
        ValueError: If the data has a BETA column: the allele to which BETA refers cannot be inferred from the reference.
    """
    if "BETA" in data.columns:
        raise ValueError(
            "The data has an effect (BETA) column but no effect allele (EA) column. The allele to which the BETA estimates refer "
            "cannot be inferred from the reference panel, so the EA/NEA columns cannot be filled. Provide the EA column "
            "(the NEA column can then be filled from the reference), or remove the BETA column."
        )
    reference = _reference_at_positions(data, reference_panel_df, ["A1", "A2"])
    ambiguous = reference.duplicated(subset=["CHR", "POS"], keep=False)
    n_ambiguous = _count_rows_in(data, ["CHR", "POS"], reference[ambiguous])
    data = data.merge(reference[~ambiguous], on=["CHR", "POS"], how="left")
    n_missing = data["A1"].isna().sum()
    data.rename(columns={"A1": "EA", "A2": "NEA"}, inplace=True)

    perc_missing = n_missing / data.shape[0] * 100
    print(
        f"Alleles columns created: effect (EA) and non-effect allele (NEA). {n_missing}({perc_missing:.3f}%) values are set to nan because SNPs were not found in the reference data"
        f"{f' or are located at a position with several reference variants ({n_ambiguous} SNPs)' if n_ambiguous else ''}."
    )
    if perc_missing > 50:
        warnings.warn(
            f"The EA (Effect Allele) and NEA (Non-Effect Allele) for many SNPs could not been found. Make sure the CHR/POS coordinates are in the same build as the reference panel (GRCh37 (hg19) for the default one). If not, you can first use the .lift() method to lift them. For instance: .lift(start='hg38', end='hg19', replace=True) if they are in build GRCh38 (hg38)."
        )
    return data


def fill_nea(data, reference_panel_df):
    """
    Fill in the NEA column based on reference data.
    The NEA is taken from the reference record at the same position that contains the EA allele.
    If several reference records at the position contain the EA allele with different other alleles, the NEA is ambiguous and left as nan.
    """
    reference = _reference_at_positions(data, reference_panel_df, ["A1", "A2"])
    candidates = data[["CHR", "POS", "EA"]].reset_index(drop=True)
    candidates["_row"] = np.arange(len(candidates))
    candidates = candidates.merge(reference, on=["CHR", "POS"], how="inner")
    candidates["NEA"] = np.select(
        [candidates["EA"] == candidates["A1"], candidates["EA"] == candidates["A2"]],
        [candidates["A2"], candidates["A1"]],
        default=None,
    )
    candidates = candidates.dropna(subset=["NEA"]).drop_duplicates(subset=["_row", "NEA"])
    n_candidates = candidates.groupby("_row")["NEA"].transform("size")
    unique = candidates[n_candidates == 1]
    nea = np.full(len(data), np.nan, dtype=object)
    nea[unique["_row"].to_numpy()] = unique["NEA"].to_numpy()
    data["NEA"] = nea
    n_missing = data["NEA"].isna().sum()
    n_ambiguous = candidates.loc[n_candidates > 1, "_row"].nunique()

    perc_missing = n_missing / data.shape[0] * 100
    print(
        f"The NEA (Non Effect Allele) column has been created. {n_missing}({perc_missing:.3f}%) values are set to nan because SNPs were not found in the reference data (position absent or EA not among the reference alleles)"
        f"{f' or their EA matches several reference variants with different alleles ({n_ambiguous} SNPs)' if n_ambiguous else ''}."
    )
    if perc_missing > 50:
        warnings.warn(
            f"The NEA (Non Effect Allele) for many SNPs could not been found. Make sure the CHR/POS coordinates are in the same build as the reference panel (GRCh37 (hg19) for the default one). If not, you can first use the .lift() method to lift them. For instance: .lift(start='hg38', end='hg19', replace=True) if they are in build GRCh38 (hg38)."
        )
    return data


def fill_coordinates_func(data, reference_panel_df):
    """Fill in the CHR/POS columns based on reference data."""
    if not "SNP" in data.columns:
        raise ValueError(
            f"The SNP column is not found in the data and is mandatory to fill coordinates!"
        )
    data.drop(columns=["CHR", "POS"], inplace=True, errors="ignore")
    reference = reference_panel_df.loc[
        reference_panel_df["SNP"].isin(data["SNP"].dropna().unique()), ["CHR", "POS", "SNP"]
    ]
    # SNP IDs present several times in the reference (at different positions) are ambiguous
    ambiguous = reference.duplicated(subset=["SNP"], keep=False)
    n_ambiguous = data["SNP"].isin(reference.loc[ambiguous, "SNP"]).sum()
    data = data.merge(reference[~ambiguous], on="SNP", how="left")
    n_missing = data["CHR"].isna().sum()
    data["CHR"] = data["CHR"].astype("Int64")
    data["POS"] = data["POS"].astype("Int64")
    perc_missing = n_missing / data.shape[0] * 100
    print(
        f"The coordinates columns (CHR for chromosome and POS for position) have been created. {n_missing}({perc_missing:.3f}%) SNPs were not found in the reference data"
        f"{f' or have several positions in the reference data ({n_ambiguous} SNPs)' if n_ambiguous else ''} and their values set to nan."
    )
    if perc_missing > 50:
        warnings.warn(
            f"The coordinates (CHR/POS) of many SNPs could not be found. Make sure the SNP column contains rsIDs matching the IDs of the reference panel. "
            f"Note that the coordinates are given in the build of the reference panel (GRCh37 (hg19) for the default one, use reference_panel = 38 for GRCh38 (hg38))."
        )
    return data


def fill_snpids_func(data, reference_panel_df, keep_indel):
    """
    Fill in the SNP column based on reference data.
    If EA and NEA are present, uses allele matching for more accurate SNP assignment.
    If some SNPids are still missing, they will be replaced by a standard name: CHR:POS:NEA:EA
    """
    for column in ["CHR", "POS"]:
        if not (column in data.columns):
            raise ValueError(
                f"The column {column} is not found in the data and is mandatory to fill snpID!"
            )
            
    data.drop(columns=["SNP"], inplace=True, errors="ignore")
    
    # If both EA and NEA are present, use allele matching
    if "EA" in data.columns and "NEA" in data.columns:
        # First preprocess the allele columns
        check_allele_column(data, "EA", keep_indel=keep_indel)
        check_allele_column(data, "NEA", keep_indel=keep_indel)
        
        # Create allele_min and allele_max columns (alleles in alphabetical order) for data and reference
        data['allele_min'] = np.minimum(data['EA'].fillna(""), data['NEA'].fillna(""))
        data['allele_max'] = np.maximum(data['EA'].fillna(""), data['NEA'].fillna(""))
        reference = _reference_at_positions(data, reference_panel_df, ["SNP", "A1", "A2"])
        reference['allele_min'] = np.minimum(reference['A1'].fillna(""), reference['A2'].fillna(""))
        reference['allele_max'] = np.maximum(reference['A1'].fillna(""), reference['A2'].fillna(""))
        # Reference records with the same position and alleles are the same variant: keep the first ID
        reference = reference.drop_duplicates(subset=["CHR", "POS", "allele_min", "allele_max"], keep="first")

        # Merge using position and alleles
        data = pd.merge(
            data,
            reference[["CHR", "POS", "SNP", "allele_min", "allele_max"]],
            on=["CHR", "POS", "allele_min", "allele_max"],
            how="left"
        )

        # Clean up temporary columns
        data.drop(columns=["allele_min", "allele_max"], inplace=True)
        n_ambiguous = 0

    else:
        # Merge only on position. Positions with several reference records (multi-allelic sites) are ambiguous.
        reference = _reference_at_positions(data, reference_panel_df, ["SNP"])
        ambiguous = reference.duplicated(subset=["CHR", "POS"], keep=False)
        n_ambiguous = _count_rows_in(data, ["CHR", "POS"], reference[ambiguous])
        data = pd.merge(
            data,
            reference[~ambiguous],
            on=["CHR", "POS"],
            how="left"
        )

    n_missing = data["SNP"].isna().sum()

    # Create standard names for missing SNPs if possible: CHR:POS:allele1:allele2 with the alleles in alphabetical order,
    # so that a variant gets the same ID whatever its effect allele.
    standard_name_condition = "EA" in data.columns and "NEA" in data.columns and n_missing > 0
    if standard_name_condition:
        missing_snp_condition = (
            data["SNP"].isna() & data["CHR"].notna() & data["POS"].notna() & data["EA"].notna() & data["NEA"].notna()
        )
        allele_1 = np.minimum(data.loc[missing_snp_condition, "EA"], data.loc[missing_snp_condition, "NEA"])
        allele_2 = np.maximum(data.loc[missing_snp_condition, "EA"], data.loc[missing_snp_condition, "NEA"])
        data.loc[missing_snp_condition, "SNP"] = (
            data.loc[missing_snp_condition, "CHR"].astype(str)
            + ":"
            + data.loc[missing_snp_condition, "POS"].astype(str)
            + ":"
            + allele_1
            + ":"
            + allele_2
        )
        n_named = missing_snp_condition.sum()
        print_statement = (
            f". {n_named} of them have their ID set to CHR:POS:allele1:allele2 (alleles in alphabetical order)"
            f"{', the others have missing CHR, POS, EA or NEA values and keep a nan ID' if n_named < n_missing else ''}"
        )

    perc_missing = n_missing / data.shape[0] * 100

    if n_missing == 0:
        print(
            f"The SNP column (rsID) has been created. All SNPs were found in the reference data."
        )
    else:
        print(
            f"The SNP column (rsID) has been created. {n_missing}({perc_missing:.3f}%) SNPs were not found in the reference data"
            f"{f' or are located at a position with several reference variants ({n_ambiguous} SNPs)' if n_ambiguous else ''}"
            f"{print_statement if standard_name_condition else ''}."
        )
        
    if perc_missing > 50:
        warnings.warn(
            f"The SNPid for many SNPs could not been found. Make sure you are using a reference panel in the same genome build as your data. The one used by default is GRCh37 (hg19). \n"
            f"You can use the GRCh38 (hg38) reference panel by setting reference_panel = 38."
        )

    return data

def parse_int_codes(values, chromosome=False):
    """
    Parse CHR values to non-negative integers or POS values to positive integers (Int64), invalid values become NA.

    Values are parsed numerically first (e.g. 1, "1", "1.0", "1e+05"). With chromosome=True, strings with a 'chr'
    prefix followed only by digits are accepted (e.g. "chr1"), and the chromosome names X, Y, XY, MT and M
    (with or without 'chr', case-insensitive) are converted to their numeric codes (constants.CHR_CODES).
    Anything else (e.g. "rs123", "chr6_qbl_hap6", "12,345", negative or non-integer numbers) becomes NA, as well as
    positions equal to 0 (positions are 1-based).

    Returns:
        tuple: (pd.Series of dtype Int64, number of chromosome names converted to codes)
    """
    minimum = 0 if chromosome else 1
    if pd.api.types.is_integer_dtype(values.dtype):
        parsed = values.astype("Int64")
        return parsed.where(parsed >= minimum).astype("Int64"), 0
    numeric = pd.to_numeric(values, errors="coerce").astype("float64")
    is_text = numeric.isna() & values.notna()
    n_codes = 0
    if is_text.any():
        text = values[is_text].astype(str).str.strip().str.upper()
        if chromosome:
            text = text.str.replace(r"^CHR", "", regex=True)
        parsed_text = pd.to_numeric(text.where(text.str.fullmatch(r"\d+")), errors="coerce")
        if chromosome:
            codes = text.map(CHR_CODES)
            n_codes = int(codes.notna().sum())
            parsed_text = parsed_text.fillna(codes)
        numeric[is_text] = parsed_text.to_numpy(dtype="float64")
    valid = np.isfinite(numeric) & (numeric >= minimum) & (numeric == np.round(numeric))
    return numeric.where(valid).astype("Int64"), n_codes


def check_int_column(data, int_col):
    """
    Set the type of the int_col column to Int64 and invalid values to NA. This function is used to check the validity of the CHR and POS columns.
    See :func:`parse_int_codes` for the accepted values (chromosome names X, Y, XY, MT/M are converted to 23, 24, 25, 26 in the CHR column).
    """
    nrows = data.shape[0]
    n_missing = data[int_col].isna().sum()
    n_zero = (pd.to_numeric(data[int_col], errors="coerce") == 0).sum() if int_col == "POS" else 0
    data[int_col], n_codes = parse_int_codes(data[int_col], chromosome=(int_col == "CHR"))
    if n_codes > 0:
        print(
            f"{n_codes}({n_codes/nrows*100:.3f}%) values of the CHR column are chromosome names and have been converted to their numeric codes (X=23, Y=24, XY=25, MT=26)."
        )
    n_nan = data[int_col].isna().sum()
    if n_nan > 0:
        print(
            f"The {int_col} column contains {n_nan}({n_nan/nrows*100:.3f}%) NaN values: {n_missing} missing and {n_nan - n_missing} set to NaN for not being valid "
            f"{'positive integers (including ' + str(n_zero) + ' positions equal to 0)' if int_col == 'POS' else 'non-negative integers'}."
        )
    return

def adjust_column_names(data, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns, FSTAT):
    """
    Rename columns to the standard names making sure that there are no duplicated names.
    Delete other columns if keep_columns=False, keep them if True.
    """

    # Check keep_columns argument
    if not isinstance(keep_columns, bool):
        raise TypeError(f"{keep_columns} only accepts values: True or False.")

    rename_dict = {
        CHR: "CHR",
        POS: "POS",
        SNP: "SNP",
        EA: "EA",
        NEA: "NEA",
        BETA: "BETA",
        SE: "SE",
        P: "P",
        EAF: "EAF",
        FSTAT: "FSTAT",
    }
    for key, value in rename_dict.items():
        if key != value and key not in data.columns:
            raise TypeError(f"Column {key} is not found in the dataframe.")

    if not keep_columns:
        cols_to_drop = [col for col in data.columns if col not in rename_dict.keys()]
        data.drop(columns=cols_to_drop, inplace=True)

    data.rename(columns=rename_dict, inplace=True)

    # Check duplicated column names
    column_counts = Counter(data.columns)
    duplicated_columns = [
        col
        for col, count in column_counts.items()
        if (count > 1) and (col in rename_dict.values())
    ]
    if duplicated_columns:
        raise ValueError(
            f"After adjusting the column names, the resulting dataframe has duplicated columns. Make sure your dataframe does not have a different column named {duplicated_columns}."
        )
    
    # Print columns found
    print(f"The following columns were found: {list(set(rename_dict.values()) & set(data.columns.to_list()))}")
    
    return data


def check_arguments(
    preprocessing,
    effect_column,
    fill_snpids,
    fill_coordinates,
    keep_indel,
    keep_dups,
):
    """
    Verify the arguments passed for the Geno initialization and apply logic based on the preprocessing value. See :class:`Geno` for more details.

    Returns:
        tuple: Tuple containing updated values for (keep_columns, keep_indel, keep_dups, fill_snpids, fill_coordinates)

    Raises:
        TypeError: For invalid data types or incompatible argument values.
    """

    # Validate preprocessing value
    if preprocessing not in ['None', 'Fill', 'Fill_delete']:
        raise ValueError(
            "preprocessing must be one of ['None', 'Fill', 'Fill_delete']. Refer to the Geno class docstring for details."
        )

    # Validate effect_column value
    if not ((effect_column is None) or (effect_column in ["OR", "BETA"])):
        raise ValueError("effect_column must be one of [None, 'OR', 'BETA'].")

    # Ensure all other arguments are either None or boolean type
    variables = {
        "fill_snpids": fill_snpids,
        "fill_coordinates": fill_coordinates,
        "keep_indel": keep_indel,
        "keep_dups": keep_dups,
    }
    for name, value in variables.items():
        if not (value is None or isinstance(value, bool)):
            raise TypeError(f"{name} only accepts values: None, True, or False.")

    # Helper functions for preprocessing logic
    def keeptype_column(arg):
        """Helper function to decide whether to keep indels/duplicates."""
        return True if arg is None and preprocessing in ['None', 'Fill'] else arg

    def filltype_column(arg):
        """Helper function to decide whether to fill snpids/coordinates."""
        return False if arg is None and preprocessing == 'None' else arg

    # Apply preprocessing logic
    keep_indel = keeptype_column(keep_indel)
    keep_dups = keeptype_column(keep_dups)
    fill_snpids = filltype_column(fill_snpids)
    fill_coordinates = filltype_column(fill_coordinates)

    return keep_indel, keep_dups, fill_snpids, fill_coordinates


def save_data(data, name, path="", fmt="parquet", sep="\t", header=True):
    """
    Save data to a specified file format.

    Supported formats: .parquet (default), .csv, .txt.
    Future supported formats: .vcf, .vcf.gz.

    Args:
        data (pd.DataFrame): DataFrame to be saved.
        name (str): A unique identifier for the data, used as the filename.
        path (str, optional): The directory where the file will be saved. Defaults to current directory.
        fmt (str, optional): The desired file format. Defaults to "parquet".
        sep (str, optional): Delimiter for text-based formats (.csv, .txt). Defaults to tab.
        header (bool, optional): Whether to include column names in text-based formats. Defaults to True.
    """
    path = os.path.join(path, name)
    if fmt == "parquet":
        data.to_parquet(
            f"{path}.parquet", engine="pyarrow", compression="brotli", index=True
        )
    elif fmt == "csv":
        data.to_csv(f"{path}.csv", sep=sep, header=header, index=False)
    elif fmt == "txt":
        data.to_csv(f"{path}.txt", sep=sep, header=header, index=False)
    else:
        raise ValueError(
            f"Format {fmt} is not supported. Supported formats are: parquet, csv, txt."
        )


def Combine_Geno(Gs):
    """
    Combine multiple Geno instances.

    Args:
        Gs (list): A list of Geno instances to combine.

    Returns:
        Geno: A new Geno instance containing the combined data.
    """
    from .Geno import Geno
    
    C = pd.DataFrame()

    for G in Gs:
        C = pd.concat([C, G.data])

    C = C.reset_index(drop=True)

    return Geno(C)


def update_eaf_func(data, reference_panel, object_name, ram=10000, fill=True):
    """
    Core logic to update or create the EAF (Effect Allele Frequency) column.

    This function calculates EAF from a reference panel. If CHR/POS are available,
    it uses a fast, coordinate-based extraction with PLINK. Otherwise, it falls
    back to SNP-ID-based extraction.

    A reference variant is assigned to a row only if its alleles match: {EA, NEA} = {REF, ALT} when the NEA
    column is present, EA in {REF, ALT} otherwise. Rows matching several reference variants with different
    frequencies are ambiguous and left unmatched.
    """
    ref_panel_path, ref_filetype = get_reference_panel_path(reference_panel)
    create_tmp()

    by_coordinate = "CHR" in data.columns and "POS" in data.columns
    freq_prefix = os.path.join("tmp_GENAL", f"{object_name}_eaf_freqs")
    plink_command = [
        get_plink_path(),
        "--pfile" if ref_filetype == "pgen" else "--bfile", ref_panel_path,
        "--memory", str(ram),
    ]

    # --- Match by CHR/POS or SNP ID ---
    if by_coordinate:
        print("CHR/POS columns present. SNPs searched based on genomic positions.")
        keys = ["CHR", "POS"]
        query = pd.DataFrame({
            "CHR": parse_int_codes(data["CHR"], chromosome=True)[0],
            "POS": parse_int_codes(data["POS"])[0],
        })

        # Write integer coordinates to a temp file for PLINK's --extract range, then run --freq adding POS to the output
        coord_path = os.path.join("tmp_GENAL", f"{object_name}_coord_list.txt")
        coords = query.dropna().drop_duplicates().astype("int64")
        coords[["CHR", "POS", "POS"]].to_csv(coord_path, sep='\t', index=False, header=False)
        plink_command += ["--extract", "range", coord_path, "--freq", "cols=+pos", "--out", freq_prefix]

    else:
        print("Using SNP IDs to extract frequencies.")
        if "SNP" not in data.columns:
            raise ValueError("SNP column is required when CHR/POS are not available.")
        keys = ["SNP"]
        query = pd.DataFrame({"SNP": data["SNP"].to_numpy()})

        snp_list_path = os.path.join("tmp_GENAL", f"{object_name}_snp_list.txt")
        data[["SNP"]].dropna().to_csv(snp_list_path, index=False, header=False)
        plink_command += ["--extract", snp_list_path, "--freq", "--out", freq_prefix]

    if query.dropna().empty:
        warnings.warn(f"No variant has a valid {'/'.join(keys)} value. Cannot update EAF.")
        return data
    run_plink_command(plink_command)

    freq_path = f"{freq_prefix}.afreq"
    if not os.path.exists(freq_path) or os.path.getsize(freq_path) == 0:
        warnings.warn("No variants from your data were found in the reference panel, or PLINK failed. Cannot update EAF.")
        return data

    freqs_df = read_afreq(freq_path)
    query["EA"] = data["EA"].astype(str).str.upper().to_numpy()
    if "NEA" in data.columns:
        query["NEA"] = data["NEA"].astype(str).str.upper().to_numpy()
    eaf_new, n_ambiguous, n_found, n_concordant = match_reference_eaf(query, freqs_df, keys)

    # Report the matching: coverage (variants absent from the panel) is informational, while a low allele concordance
    # among the variants found in the panel indicates a genome build mismatch (or a different allele coding).
    n_query = query[keys].notna().all(axis=1).sum()
    n_matched = int(np.isfinite(eaf_new).sum())
    print(
        f"The EAF of {n_matched}({n_matched/max(n_query, 1)*100:.3f}%) variants has been found in the reference panel "
        f"(matching {'/'.join(keys)} and {'EA/NEA alleles' if 'NEA' in query.columns else 'the EA allele'})."
    )
    if n_query - n_found > 0:
        absent = "have no reference variant at their position" if keys == ["CHR", "POS"] else "have a SNP ID not found in the reference panel"
        print(f"{n_query - n_found}({(n_query - n_found)/max(n_query, 1)*100:.3f}%) variants {absent}.")
    if n_found - n_concordant > 0:
        print(f"{n_found - n_concordant} variants found in the reference panel have alleles that do not match the reference variant(s).")
    if n_ambiguous > 0:
        print(f"{n_ambiguous} variants match several reference variants with different frequencies and are left unmatched.")
    if n_found > 0 and n_concordant / n_found < 0.9:
        found = "with a reference variant at their position" if keys == ["CHR", "POS"] else "with their SNP ID found in the reference panel"
        warnings.warn(
            f"Only {n_concordant/n_found*100:.1f}% of the {n_found} variants {found} match the alleles of the reference variant. "
            f"The data and the reference panel are likely on different genome builds"
            f"{' (or use a different allele coding)' if keys == ['SNP'] else ''}: "
            f"the EUR_37 default panel is in GRCh37 (hg19), use e.g. reference_panel='EUR_38' for GRCh38 (hg38)."
        )

    # Create updated EAF column
    if 'EAF' not in data.columns:
        data['EAF'] = np.nan
    if fill:
        data['EAF'] = np.where(np.isfinite(eaf_new), eaf_new, data['EAF'])
    else:
        data['EAF'] = eaf_new

    return data


def read_afreq(freq_path):
    """
    Read a PLINK2 .afreq file into a DataFrame with columns CHR, (POS), SNP, REF, ALT, ALT_FREQ, REF_FREQ.
    Multi-allelic records (ALT "A,T", ALT_FREQS "0.1,0.2") are split into one row per ALT allele,
    REF_FREQ being the frequency of the REF allele (1 - sum of the ALT frequencies).
    """
    freqs = pd.read_csv(freq_path, sep="\t", dtype={"#CHROM": str, "ID": str, "REF": str, "ALT": str, "ALT_FREQS": str})
    freqs = freqs.rename(columns={"#CHROM": "CHR", "ID": "SNP"})
    freqs["CHR"] = parse_int_codes(freqs["CHR"], chromosome=True)[0]
    freqs["ALT_FREQ"] = pd.to_numeric(freqs["ALT_FREQS"], errors="coerce")
    freqs["REF_FREQ"] = 1 - freqs["ALT_FREQ"]
    multi_allelic = freqs["ALT"].str.contains(",", na=False)
    if multi_allelic.any():
        split = freqs[multi_allelic].copy()
        split["ALT"] = split["ALT"].str.split(",")
        split["ALT_FREQ"] = split["ALT_FREQS"].str.split(",")
        split = split.explode(["ALT", "ALT_FREQ"])
        split["ALT_FREQ"] = pd.to_numeric(split["ALT_FREQ"], errors="coerce")
        split["REF_FREQ"] = 1 - split.groupby(level=0)["ALT_FREQ"].transform(lambda f: f.sum(skipna=False))
        freqs = pd.concat([freqs[~multi_allelic], split], ignore_index=True)
    return freqs.drop(columns=["ALT_FREQS"])


def match_reference_eaf(query, freqs, keys):
    """
    Assign to each row of query (keys, EA and optionally NEA columns) the frequency of its EA allele in freqs (see :func:`read_afreq`).
    A reference variant matches if {EA, NEA} = {REF, ALT} (or EA in {REF, ALT} without a NEA column).
    Rows matching several reference variants with different frequencies are left as nan.

    Returns:
        tuple: (np.ndarray of EAF values aligned with query rows, number of ambiguous rows, number of rows with at least
            one reference variant with the same keys, number of these rows whose alleles match at least one of them)
    """
    rows = query.reset_index(drop=True)
    rows["_row"] = np.arange(len(rows))
    candidates = rows.dropna(subset=keys).merge(
        freqs[keys + ["REF", "ALT", "ALT_FREQ", "REF_FREQ"]], on=keys, how="inner"
    )
    is_alt = candidates["EA"] == candidates["ALT"]
    is_ref = candidates["EA"] == candidates["REF"]
    if "NEA" in candidates.columns:
        is_alt &= candidates["NEA"] == candidates["REF"]
        is_ref &= candidates["NEA"] == candidates["ALT"]
    candidates["EAF"] = np.select([is_alt, is_ref], [candidates["ALT_FREQ"], candidates["REF_FREQ"]], default=np.nan)
    n_found = int(candidates["_row"].nunique())
    candidates = candidates[is_alt | is_ref]
    n_concordant = int(candidates["_row"].nunique())
    n_values = candidates.groupby("_row")["EAF"].nunique()
    unique_rows = n_values.index[n_values == 1]
    values = candidates.groupby("_row")["EAF"].first()
    eaf = np.full(len(rows), np.nan)
    eaf[unique_rows.to_numpy()] = values.loc[unique_rows].to_numpy()
    return eaf, int((n_values > 1).sum()), n_found, n_concordant
