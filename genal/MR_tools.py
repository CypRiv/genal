import pandas as pd
import numpy as np
import datetime
import os, subprocess
from pandas.api.types import is_numeric_dtype
from tqdm import tqdm

from .proxy import find_proxies, query_outcome_proxy
from .MR import *
from .MRpresso import mr_presso
from .constants import MR_METHODS_NAMES

REQUIRED_COLUMNS = ["SNP", "BETA", "SE", "EA", "NEA"]
RESULT_COLS = ["exposure", "outcome", "method", "nSNP", "b", "se", "pval", "Q", "Q_df", "Q_pval"]

def _ensure_mr_result_cols(result_df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure standard MR result columns exist and are type-stable.
    """
    required = ["method", "nSNP", "b", "se", "pval", "Q", "Q_df", "Q_pval"]
    for col in required:
        if col not in result_df.columns:
            result_df[col] = np.nan

    q_df = pd.to_numeric(result_df["Q_df"], errors="coerce")
    intlike = q_df.isna() | np.isclose(q_df, np.round(q_df))
    result_df["Q_df"] = np.round(q_df.where(intlike)).astype("Int64")
    return result_df


def mrpresso_func(
    data,
    action,
    eaf_threshold,
    n_iterations,
    outlier_test,
    distortion_test,
    significance_p,
    cpus,
):
    """
    Wrapper function corresponding to the :meth:`Geno.MRpresso` method.
    The MR-PRESSO algorithm is implemented here: :func:`MRpresso.mr_presso`
    Refer to them for more details regarding arguments and return values.

    Notes:
        - EAF column check if action is set to 2.
        - Data harmonization between exposure and outcome data based on action and eaf_threshold
        - NA check
        - MRpresso call and results return
    """

    # Check that action argument is a correct input
    if action not in [1, 2, 3]:
        raise ValueError("The action argument only takes 1,2 or 3 as value")

    # Unpack data (coming from MR_data attribute)
    df_exposure = data[0]
    df_outcome = data[1]
    name_outcome = data[2]

    # Check number of instruments
    if df_exposure.shape[0] < 3:
        print("Not enough instruments to run MRpresso. At least 3 are required.")
        return pd.DataFrame(), {}, pd.DataFrame(), {}, None
    
    # Check EAF columns if action = 2
    if action == 2:
        if "EAF" not in df_exposure.columns:
            print(
                "Warning: action = 2 but EAF column is missing from exposure data: palindromic SNPs will be deleted (action set to 3)."
            )
            action = 3
        elif "EAF" not in df_outcome.columns:
            print(
                "Warning: action = 2 but EAF column is missing from outcome data: palindromic SNPs will be deleted (action set to 3)."
            )
            action = 3

    # Harmonize exposure and outcome data
    df_mr = harmonize_MR(
        df_exposure, df_outcome, action=action, eaf_threshold=eaf_threshold
    )
    
    df_mr = df_mr_formatting(df_mr)
    if "SNP" not in df_mr.columns:
        raise ValueError("MR dataframe must contain a 'SNP' column after formatting.")

    df_mr = df_mr.set_index("SNP", drop=False)
    if df_mr.index.duplicated().any():
        duplicated_snps = df_mr.index[df_mr.index.duplicated()].unique().tolist()
        print(
            f"Warning: duplicated SNP IDs found in MR dataframe (keeping first occurrence): {duplicated_snps}"
        )
        df_mr = df_mr.loc[~df_mr.index.duplicated(keep="first")]

    # Call and return the results of MR-PRESSO
    return mr_presso(
        df_mr,
        ["BETA_e"],
        n_iterations,
        outlier_test,
        distortion_test,
        significance_p,
        cpus,
    )


def MR_func(
    data,
    methods,
    action,
    eaf_threshold,
    nboot,
    penk,
    phi,
    name_exposure,
    cpus,
    subset_data
):
    """
    Wrapper function corresponding to the :meth:`Geno.MR` method. Refer to them for more details regarding arguments and return values.
    The MR algorithms are implemented here: :func:`MR.mr_ivw`, :func:`MR.mr_weighted_median`, :func:`MR.mr_egger_regression`, :func:`MR.mr_simple_median`...

    Notes:
        - Validation of the action and methods arguments
        - EAF column check if action is set to 2.
        - Data harmonization between exposure and outcome data based on action and eaf_threshold
        - NA check
        - MR methods execution
        - Compiles results and return a pd.DataFrame
    """
    # Check the methods argument (contains either MR method names or "all")
    valid_methods = list(MR_METHODS_NAMES.keys())
    valid_methods.append("all")
    methods = methods if isinstance(methods, list) else [methods]
    if not all(m in valid_methods for m in methods):
        raise ValueError(
            f"The list of methods can only contain strings in {valid_methods}"
        )

    # If subset_data is provided, skip harmonization and use it directly
    if subset_data is not None:
        df_mr = subset_data
        name_outcome = data[2]
    else:
        # Check that action argument is a correct input
        if action not in [1, 2, 3]:
            raise ValueError("The action argument only takes 1,2 or 3 as value")

        # Unpack data (coming from MR_data attribute)
        df_exposure = data[0]
        df_outcome = data[1]
        name_outcome = data[2]

        # Check number of instruments
        if df_exposure.shape[0] < 2:
            print("Not enough instruments to run MR. At least 2 are required.")
            return pd.DataFrame(columns=RESULT_COLS), pd.DataFrame()

        # Check EAF columns if action = 2
        if action == 2:
            if "EAF" not in df_exposure.columns:
                print(
                    "Warning: action = 2 but EAF column is missing from exposure data: palindromic SNPs will be deleted (action set to 3)."
                )
                action = 3
            elif "EAF" not in df_outcome.columns:
                print(
                    "Warning: action = 2 but EAF column is missing from outcome data: palindromic SNPs will be deleted (action set to 3)."
                )
                action = 3

        # Harmonize exposure and outcome data
        df_mr = harmonize_MR(
            df_exposure, df_outcome, action=action, eaf_threshold=eaf_threshold
        )

        df_mr = df_mr_formatting(df_mr)

        # Check number of remaining instruments
        n_snps = df_mr.shape[0]
        if n_snps < 2:
            print(f"{n_snps} SNPs remaining after harmonization step but at least 2 are required to run MR.")
            return pd.DataFrame(columns=RESULT_COLS), df_mr

    # Prepare values for MR methods
    BETA_e, BETA_o, SE_e, SE_o = (
        df_mr["BETA_e"],
        df_mr["BETA_o"],
        df_mr["SE_e"],
        df_mr["SE_o"],
    )  
    
    print(
        f"Running Mendelian Randomization with {name_exposure} as exposure and {name_outcome} as outcome."
    )

    # Mapping the methods passed as argument to the corresponding functions and freeze arguments
    FUNCTION_MAP = {
        "Egger": partial(mr_egger_regression, BETA_e, SE_e, BETA_o, SE_o),
        "Egger-boot": partial(
            mr_egger_regression_bootstrap, BETA_e, SE_e, BETA_o, SE_o, nboot, cpus
        ),
        "WM": partial(mr_weighted_median, BETA_e, SE_e, BETA_o, SE_o, nboot),
        "WM-pen": partial(mr_pen_wm, BETA_e, SE_e, BETA_o, SE_o, nboot, penk),
        "Simple-median": partial(mr_simple_median, BETA_e, SE_e, BETA_o, SE_o, nboot),
        "IVW": partial(mr_ivw, BETA_e, SE_e, BETA_o, SE_o),
        "IVW-RE": partial(mr_ivw_re, BETA_e, SE_e, BETA_o, SE_o),
        "IVW-FE": partial(mr_ivw_fe, BETA_e, SE_e, BETA_o, SE_o),
        "UWR": partial(mr_uwr, BETA_e, SE_e, BETA_o, SE_o),
        "Sign": partial(mr_sign, BETA_e, BETA_o),
        "Simple-mode": partial(mr_simple_mode, BETA_e, SE_e, BETA_o, SE_o, phi, nboot, cpus),
        "Weighted-mode": partial(mr_weighted_mode, BETA_e, SE_e, BETA_o, SE_o, phi, nboot, cpus),
    }

    # Compute required MR methods and gather results
    results = []
    if "all" in methods:
        methods = list(MR_METHODS_NAMES.keys())
    for method in methods:
        func = FUNCTION_MAP.get(method, None)
        result = func()
        results.extend(result)

    res = pd.DataFrame(results)
    res = _ensure_mr_result_cols(res)
    res["exposure"], res["outcome"] = name_exposure, name_outcome
    res = res.reindex(columns=RESULT_COLS)

    return res, df_mr


def MR_loo_func(
    data,
    method,
    action,
    eaf_threshold,
    heterogeneity,
    nboot,
    penk,
    phi,
    name_exposure,
    cpus,
    subset_data,
):
    """
    Wrapper function corresponding to the :meth:`Geno.MR_loo` method.
    """
    valid_methods = list(MR_METHODS_NAMES.keys())
    if not isinstance(method, str) or method not in valid_methods or method == "all":
        raise ValueError(
            f"The method argument must be a single string in {valid_methods} (and cannot be 'all')."
        )

    name_outcome = data[2]

    if subset_data is not None:
        df_mr = subset_data.copy()
        if "SNP" not in df_mr.columns:
            df_mr["SNP"] = df_mr.index
    else:
        if action not in [1, 2, 3]:
            raise ValueError("The action argument only takes 1,2 or 3 as value")

        df_exposure = data[0]
        df_outcome = data[1]
        name_outcome = data[2]

        if action == 2:
            if "EAF" not in df_exposure.columns:
                print(
                    "Warning: action = 2 but EAF column is missing from exposure data: palindromic SNPs will be deleted (action set to 3)."
                )
                action = 3
            elif "EAF" not in df_outcome.columns:
                print(
                    "Warning: action = 2 but EAF column is missing from outcome data: palindromic SNPs will be deleted (action set to 3)."
                )
                action = 3

        df_mr = harmonize_MR(
            df_exposure, df_outcome, action=action, eaf_threshold=eaf_threshold
        )
        df_mr = df_mr_formatting(df_mr)

    if "SNP" not in df_mr.columns:
        raise ValueError("MR dataframe must contain a 'SNP' column for leave-one-out.")

    df_mr = df_mr.set_index("SNP", drop=False)
    if df_mr.index.duplicated().any():
        duplicated_snps = df_mr.index[df_mr.index.duplicated()].unique().tolist()
        print(
            f"Warning: duplicated SNP IDs found in MR dataframe (keeping first occurrence): {duplicated_snps}"
        )
        df_mr = df_mr.loc[~df_mr.index.duplicated(keep="first")]

    method_display_name = (
        MR_METHODS_NAMES[method][0]
        if isinstance(MR_METHODS_NAMES[method], (list, tuple))
        else MR_METHODS_NAMES[method]
    )

    if len(df_mr) < 3:
        print(
            "Not enough instruments to run leave-one-out MR. At least 3 are required."
        )
        all_row = {
            "method": method_display_name,
            "nSNP": np.nan,
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "Q": np.nan,
            "Q_df": np.nan,
            "Q_pval": np.nan,
        }
        return pd.DataFrame(), all_row, method_display_name

    def _run_method(df_sub: pd.DataFrame):
        BETA_e, BETA_o, SE_e, SE_o = (
            df_sub["BETA_e"],
            df_sub["BETA_o"],
            df_sub["SE_e"],
            df_sub["SE_o"],
        )

        if method == "Egger":
            return mr_egger_regression(BETA_e, SE_e, BETA_o, SE_o)
        if method == "Egger-boot":
            return mr_egger_regression_bootstrap(
                BETA_e, SE_e, BETA_o, SE_o, nboot, cpus, show_progress=False
            )
        if method == "WM":
            return mr_weighted_median(BETA_e, SE_e, BETA_o, SE_o, nboot)
        if method == "WM-pen":
            return mr_pen_wm(BETA_e, SE_e, BETA_o, SE_o, nboot, penk)
        if method == "Simple-median":
            return mr_simple_median(BETA_e, SE_e, BETA_o, SE_o, nboot)
        if method == "IVW":
            return mr_ivw(BETA_e, SE_e, BETA_o, SE_o)
        if method == "IVW-RE":
            return mr_ivw_re(BETA_e, SE_e, BETA_o, SE_o)
        if method == "IVW-FE":
            return mr_ivw_fe(BETA_e, SE_e, BETA_o, SE_o)
        if method == "UWR":
            return mr_uwr(BETA_e, SE_e, BETA_o, SE_o)
        if method == "Sign":
            return mr_sign(BETA_e, BETA_o)
        if method == "Simple-mode":
            return mr_simple_mode(
                BETA_e,
                SE_e,
                BETA_o,
                SE_o,
                phi,
                nboot,
                cpus,
                show_progress=False,
            )
        if method == "Weighted-mode":
            return mr_weighted_mode(
                BETA_e,
                SE_e,
                BETA_o,
                SE_o,
                phi,
                nboot,
                cpus,
                show_progress=False,
            )

        raise ValueError(f"Method '{method}' is not supported for leave-one-out.")

    def _select_causal_row(method_results):
        if not method_results:
            return {}
        if isinstance(MR_METHODS_NAMES[method], (list, tuple)):
            causal = next(
                (r for r in method_results if r.get("method") == method_display_name),
                None,
            )
            return causal if causal is not None else {}
        return method_results[0]

    all_res = _select_causal_row(_run_method(df_mr))
    all_row = {
        "method": method_display_name,
        "nSNP": all_res.get("nSNP", len(df_mr)),
        "b": all_res.get("b", np.nan),
        "se": all_res.get("se", np.nan),
        "pval": all_res.get("pval", np.nan),
        "Q": all_res.get("Q", np.nan),
        "Q_df": all_res.get("Q_df", np.nan),
        "Q_pval": all_res.get("Q_pval", np.nan),
    }

    loo_rows = []
    for snp in tqdm(df_mr.index, desc="Leave-one-out MR", total=len(df_mr), ncols=100):
        df_sub = df_mr.drop(index=snp)
        res = _select_causal_row(_run_method(df_sub))
        loo_rows.append(
            {
                "SNP": snp,
                "nSNP": res.get("nSNP", len(df_sub)),
                "b": res.get("b", np.nan),
                "se": res.get("se", np.nan),
                "pval": res.get("pval", np.nan),
                "Q": res.get("Q", np.nan),
                "Q_df": res.get("Q_df", np.nan),
                "Q_pval": res.get("Q_pval", np.nan),
            }
        )

    loo_df = pd.DataFrame(loo_rows)
    loo_df["exposure"] = name_exposure
    loo_df["outcome"] = name_outcome
    loo_df["method"] = method_display_name
    loo_df = _ensure_mr_result_cols(loo_df)
    loo_df = loo_df.reindex(
        columns=[
            "exposure",
            "outcome",
            "method",
            "SNP",
            "nSNP",
            "b",
            "se",
            "pval",
            "Q",
            "Q_df",
            "Q_pval",
        ]
    )

    if not heterogeneity:
        loo_df = loo_df.loc[
            :, ["exposure", "outcome", "method", "SNP", "nSNP", "b", "se", "pval"]
        ]

    return loo_df, all_row, method_display_name

def df_mr_formatting(df_mr):
    """
    Function to delete invalid values from the MR dataframe (after the harmonization step)
    """
    # Delete NAs or infinite values (or null values in SE columns, null values in BETA are accepted) and print the SNP names and if the invalid value came from exposure or outcome data.
    df_mr = df_mr[["SNP", "BETA_e", "SE_e", "BETA_o", "SE_o"]].copy()
    df_mr.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_mr.loc[:, ["SE_e", "SE_o"]] = df_mr.loc[:, ["SE_e", "SE_o"]].replace(0, np.nan)
    mask_exposure = df_mr[["BETA_e", "SE_e"]].isna().any(axis=1)
    mask_outcome = df_mr[["BETA_o", "SE_o"]].isna().any(axis=1)
    rows_to_delete_exposure = df_mr[mask_exposure]
    rows_to_delete_outcome = df_mr[mask_outcome]
    n_deleted_exposure = len(rows_to_delete_exposure)
    n_deleted_outcome = len(rows_to_delete_outcome)
    if n_deleted_exposure > 0:
        print(
            f"Deleting {n_deleted_exposure} SNPs with NA or infinite values in BETA/SE columns, or null values in SE column (exposure data): {rows_to_delete_exposure['SNP'].tolist()}"
        )
    if n_deleted_outcome > 0:
        print(
            f"Deleting {n_deleted_outcome} SNPs with NA or infinite values in BETA/SE columns, or null values in SE column (outcome data): {rows_to_delete_outcome['SNP'].tolist()}"
        )

    df_mr = df_mr.dropna(subset=["BETA_e", "SE_e", "BETA_o", "SE_o"]).reset_index(
        drop=True
    )

    return df_mr.loc[:, ["SNP", "BETA_e", "SE_e", "BETA_o", "SE_o"]]



def harmonize_MR(df_exposure, df_outcome, action=2, eaf_threshold=0.42):
    """
    Harmonize exposure and outcome for MR analyses.

    Parameters:
        - df_exposure (pd.DataFrame): Exposure data with "SNP","BETA","SE","EA","NEA" and "EAF" if action=2
        - df_outcome (pd.DataFrame): Outcome data with "SNP","BETA","SE","EA","NEA" and "EAF" if action=2
        - action (int, optional): Determines how to treat palindromes. Defaults to 2.
            1: Doesn't attempt to flip them (= Assume all alleles are coded on the forward strand)
            2: Use allele frequencies (EAF) to attempt to flip them (conservative, default)
            3: Remove all palindromic SNPs (very conservative).
        - eaf_threshold (float, optional): Maximal effect allele frequency accepted when attempting to flip palindromic SNPs (only applied if action = 2). Defaults to 0.42.

    Returns:
        - pd.DataFrame: Harmonized data.

    Notes:
        - Verify the presence of required columns in both dataframes and rename them
        - Merge exposure and outcome data
        - Identify palindromes
        - Classify SNPs into aligned / inverted / need to be flipped
        - Flip the ones that require flipping
        - Switch those that are inverted to align them
        - Remove those that are still not aligned
        - Treat palindromes based on action parameter
    """

    # Check required columns in both dataframes
    check_required_columns(df_exposure, REQUIRED_COLUMNS)
    check_required_columns(df_outcome, REQUIRED_COLUMNS)

    # Rename columns
    df_exposure = df_exposure.rename(
        columns={
            "EA": "EA_e",
            "NEA": "NEA_e",
            "EAF": "EAF_e",
            "BETA": "BETA_e",
            "SE": "SE_e",
            "P": "P_e",
        },
        errors="ignore",
    )
    df_outcome = df_outcome.rename(
        columns={
            "EA": "EA_o",
            "NEA": "NEA_o",
            "EAF": "EAF_o",
            "BETA": "BETA_o",
            "SE": "SE_o",
            "P": "P_o",
        },
        errors="ignore",
    )

    df_outcome = df_outcome[
        df_outcome.columns.intersection(
            ["SNP", "EA_o", "NEA_o", "EAF_o", "BETA_o", "SE_o", "P_o"]
        )
    ]

    # Merge the dataframes on SNP
    df = df_exposure.merge(df_outcome, on="SNP", how="left")

    # Default EAF columns if they do not exist
    df["EAF_e"] = df.get("EAF_e", 0.5)
    df["EAF_o"] = df.get("EAF_o", 0.5)

    # Identify palindromes
    condition1 = ((df["EA_e"] == "A") & (df["NEA_e"] == "T")) | (
        (df["EA_e"] == "T") & (df["NEA_e"] == "A")
    )
    condition2 = ((df["EA_e"] == "C") & (df["NEA_e"] == "G")) | (
        (df["EA_e"] == "G") & (df["NEA_e"] == "C")
    )
    df["palindrome"] = condition1 | condition2

    # Align effect alleles between exposure and outcome
    # Classify SNPs into aligned / inverted / need to be flipped
    df["aligned"] = (df.EA_e == df.EA_o) & (df.NEA_e == df.NEA_o)  # Already aligned
    df["inverted"] = (df.EA_e == df.NEA_o) & (df.NEA_e == df.EA_o)  # Inverted
    df["to_flip"] = (
        ~df["aligned"] & ~df["inverted"] & ~df["palindrome"]
    )  # Neither aligned nor inverted nor palindromic

    # Flip the SNPs to be flipped
    if df["to_flip"].sum() > 0:
        to_flip_idx = df[df["to_flip"]].index  # Get indices of SNPs to be flipped
        df.loc[to_flip_idx, "EA_o"] = flip_alleles(df.loc[to_flip_idx, "EA_o"])
        df.loc[to_flip_idx, "NEA_o"] = flip_alleles(df.loc[to_flip_idx, "NEA_o"])

    # Recheck inverted SNPS to flag those that are inverted after being flipped
    df["inverted"] = np.where(
        (df["EA_e"] == df["NEA_o"]) & (df["NEA_e"] == df["EA_o"]), True, False
    )

    # Switch the inverted SNPs to align them
    if df["inverted"].sum() > 0:
        inverted_idx = df[df["inverted"]].index  # Get indices of inverted SNPs
        df.loc[inverted_idx, ["EA_o", "NEA_o"]] = df.loc[
            inverted_idx, ["NEA_o", "EA_o"]
        ].values  # Swap outcome EA and NEA values
        df.loc[inverted_idx, "BETA_o"] *= -1  # Invert outcome BETA
        df.loc[inverted_idx, "EAF_o"] = (
            1 - df.loc[inverted_idx, "EAF_o"]
        )  # Invert outcome EAF

    # All the SNPs should be aligned at this point. If not, they have an allele mismatch and need to be removed
    df["aligned"] = (df["EA_e"] == df["EA_o"]) & (df["NEA_e"] == df["NEA_o"])  # Recheck aligned
    df["allele_mismatch"] = ~df[
        "aligned"
    ]  # If still not aligned: requires exclusion due to allele mismatch
    mismatched_snps = df[df["allele_mismatch"]].shape[0]
    if mismatched_snps > 0:
        print(
            f"{mismatched_snps} SNPs have been excluded due to a mismatch between the exposure and outcome alleles data."
        )
    df = df[~df["allele_mismatch"]]
    df.reset_index(inplace=True, drop=True)

    # Treat palindromes based on the action parameter
    if action == 3:  # Simply delete them
        snps_deleted = df[df.palindrome].SNP.values
        df = df[~df.palindrome]
        df.reset_index(drop=True, inplace=True)
        print(
            f"Action = 3: excluding {len(snps_deleted)} palindromic SNPs: {', '.join(snps_deleted)} \n"
        )
    elif action == 2:
        df = apply_action_2(df, eaf_threshold)
    elif action == 1:
        print(
            "Action = 1: Keeping all palindromic SNPs without attempting to flip them."
        )

    return df


def flip_alleles(x):
    """Flip the alleles."""
    x = x.str.upper()
    x = x.replace("C", "g").replace("G", "c").replace("A", "t").replace("T", "a")
    x = x.str.upper()
    return x


def check_required_columns(df, columns):
    """Check if the required columns are present in the dataframe."""
    missing_columns = [col for col in columns if col not in df.columns]
    if missing_columns:
        raise ValueError(
            f"The columns {', '.join(missing_columns)} are not found in the data and are necessary."
        )


def apply_action_2(df, eaf_threshold):
    """
    Use EAF_e and EAF_o to align palindromes if both EAFs are outside the intermediate allele frequency range.
        - Replace NA values in EAF columns by 0.5 (will be flagged and removed in step 3)
        - Set boundaries for intermediate allele frequencies
        - Identify palindromes that have an intermediate allele frequency and delete them
        - Among the remaining palindromes, identify the ones that need to be flipped and flip them
    """
    # If EAF is nan for a SNP, it will be removed
    df["EAF_e"] = np.where(df.EAF_e.isna(), 0.5, df.EAF_e)
    df["EAF_o"] = np.where(df.EAF_o.isna(), 0.5, df.EAF_o)

    # Set the boundaries for intermediate frequencies
    minf = np.minimum(eaf_threshold, 1 - eaf_threshold)
    maxf = 1 - minf

    # Identify palindromes that have an intermediate allele frequency and delete them
    df["ambiguous"] = df["palindrome"] & (
        ((minf <= df["EAF_e"]) & (df["EAF_e"] <= maxf))
        | ((minf <= df["EAF_o"]) & (df["EAF_o"] <= maxf))
    )
    snps_deleted = df[df.ambiguous].SNP.values
    df = df[~df.ambiguous]
    diff = len(snps_deleted)
    if diff > 0:
        print(
            f"Action = 2: {diff} SNPs excluded for being palindromic with intermediate allele frequencies: {', '.join(snps_deleted)} \n"
        )
    else:
        print(
            f"Action = 2: None of the SNPs are palindromic with intermediate allele frequency, keeping all of them."
        )

    # Identify palindromes that need to be flipped and flip them
    df.loc[:, "to_flip"] = df["palindrome"] & ((df.EAF_e - 0.5) * (df.EAF_o - 0.5) < 0)
    if df["to_flip"].sum() > 0:
        to_flip_idx = df[df["to_flip"]].index  # Get indices of SNPs to be flipped
        df.loc[to_flip_idx, "BETA_o"] *= -1  # Invert outcome BETA
        df.loc[to_flip_idx, "EAF_o"] = (
            1 - df.loc[to_flip_idx, "EAF_o"]
        )  # Invert outcome EAF
        print(f"Action = 2: {df.to_flip.sum()} palindromic SNPs have been flipped.")

    df.reset_index(drop=True, inplace=True)
    return df

### ___________________________ 
### Query outcome functions
### ___________________________

def query_outcome_func(
    data, outcome, name, proxy, reference_panel, kb, r2, window_snps, cpus
):
    """
    Wrapper function corresponding to the :meth:`Geno.query_outcome` method.
    Refer to it for more details on the arguments and return values.

    Notes:
        - Validation of the required columns
        - Load outcome data from Geno or path.
        - Identify SNPs present in the outcome data
        - Find proxies for the absent SNPs if needed
        - Return exposure dataframe, outcome dataframe, outcome name
    """
    # Check required columns in the exposure data
    for column in REQUIRED_COLUMNS:
        if column not in data.columns:
            raise ValueError(
                f"The column {column} is not found in the data and is necessary."
            )

    # Load the outcome dataframe (to be queried)
    import genal

    if isinstance(outcome, genal.Geno):
        df_outcome, name = load_outcome_from_geno_object(outcome)
    elif isinstance(outcome, str):
        df_outcome, name = load_outcome_from_filepath(outcome)
    else:
        raise ValueError(
            "You need to provide either a Geno object or filepath string to the outcome variable."
        )

    # Check necessary columns from outcome
    for column in REQUIRED_COLUMNS:
        if column not in df_outcome.columns:
            raise ValueError(
                f"The column {column} is not found in the outcome data and is necessary."
            )

    # Identify the exposure SNPs present in the outcome data
    print("Identifying the exposure SNPs present in the outcome data...")
    outcome_snps = set(df_outcome.SNP.values)
    exposure_snps = set(data.SNP.values)
    snps_present = exposure_snps & outcome_snps
    print(
        f"{len(snps_present)} SNPs out of {len(exposure_snps)} are present in the outcome data."
    )

    # Find proxies for absent SNPs if needed
    if proxy and (len(exposure_snps) - len(snps_present) > 0):
        snps_absent = exposure_snps - snps_present
        print(f"Searching proxies for {len(snps_absent)} SNPs...")
        ld = find_proxies(
            snps_absent,
            reference_panel=reference_panel,
            kb=kb,
            r2=r2,
            window_snps=window_snps,
            threads=cpus,
        )
        if isinstance(ld, pd.DataFrame) and not ld.empty:
            outcome = query_outcome_proxy(df_outcome, ld, snps_present, outcome_snps)
            exposure = data[data.SNP.isin(outcome.SNP)]
        else:
            print("No proxies found.")
            exposure = data[data.SNP.isin(snps_present)]
            outcome = df_outcome[df_outcome.SNP.isin(snps_present)]
    else:
        exposure = data[data.SNP.isin(snps_present)]
        outcome = df_outcome[df_outcome.SNP.isin(snps_present)]

    exposure.reset_index(drop=True, inplace=True)
    outcome.reset_index(drop=True, inplace=True)

    print(
        f"(Exposure data, Outcome data, Outcome name) stored in the .MR_data attribute."
    )

    return exposure, outcome, name


def load_outcome_from_geno_object(outcome):
    """Load outcome data from a Geno object."""
    df_outcome = outcome.data
    name = outcome.name
    print(f"Outcome data successfully loaded from '{name}' Geno instance.")
    return df_outcome, name


def load_outcome_from_filepath(outcome):
    """Load outcome data from a file path."""
    if not os.path.isfile(outcome):
        raise ValueError("The path provided doesn't lead to a file.")
    if not (outcome.endswith(".h5") or outcome.endswith(".hdf5")):
        raise ValueError("The file provided needs to be in .h5 or .hdf5 format.")
    df_outcome = pd.read_hdf(outcome, key="data")
    name = os.path.splitext(os.path.basename(outcome))[0]
    print(f"Outcome data successfully loaded from path provided.")
    return df_outcome, name
