import pandas as pd
import numpy as np
import os
import subprocess

from io import StringIO

from .tools import get_reference_panel_path, get_plink19_path

## TO DO: accept lists of CHR/POS instead of SNP names for these functions


def query_outcome_proxy(df, ld, snps_to_extract, snps_df=[]):
    """
    Extract the best proxies from a dataframe, as well as specific SNPs.

    Given a dataframe `df` (originating from GENO.data) and a dataframe of potential proxies
    (output from `find_proxies`), this function extracts the best proxies from `df` as well as
    the SNPs specified in `snps_to_extract`.
    This is suited for querying outcome data.

    Args:
        df (pd.DataFrame): Dataframe of SNP information with the usual GENO columns
                           (SNP, BETA, SE, EAF, EA, NEA). EAF is not necessary.
        ld (pd.DataFrame): Dataframe of proxies (output from `find_proxies`).
        snps_to_extract (list): List of SNPs to extract in addition to the proxies.
        snps_df (list, optional): List of SNPs to choose the proxy from. Should be the list of
                                  SNPs in df. Can be provided to avoid recomputing it. Defaults to an empty list.

    Returns:
        pd.DataFrame: Dataframe with queried SNPs and their proxies.
    """

    # If snps_df is empty, populate it with SNPs from df
    if not snps_df:
        snps_df = df.SNP.values

    # Filter proxies that are present in df
    ld = ld[ld.SNP_B.isin(snps_df)]

    # Remove original SNPs
    ld = ld[ld["SNP_A"] != ld["SNP_B"]]

    # Sort by r2 and select the best proxy for each SNP
    ld = ld.reindex(ld["R"].abs().sort_values(ascending=False).index)
    ld = ld.groupby("SNP_A").first().reset_index(drop=False)

    # Determine SNPs to query
    snps_to_query = set(snps_to_extract) | set(ld.SNP_B.values)
    df_queried = df[df.SNP.isin(snps_to_query)]

    # Merge dataframes and identify proxies
    output = df_queried.merge(ld, how="left", left_on="SNP", right_on="SNP_B")
    output["proxy"] = output["SNP_B"].notnull()

    # Flip BETA if the proxied SNP alleles are switched in the reference panel
    conditions = [
        (output["EA"] == output["B2"]),
        (output["EA"] == output["B1"]),
        (~output["proxy"]),
        (
            (output["EA"].isin([output["B1"], output["B2"]]) == False)
            & (output["proxy"])
        ),
    ]
    choices = [
        -output["BETA"],  # if EA == B2, flip the sign of BETA
        output["BETA"],  # if EA == B1, BETA does not change
        output["BETA"],  # if the original SNP was not proxied, BETA does not change
        np.nan,  # if the original SNP was proxied but EA is neither "B1" nor "B2", BETA is NaN
    ]
    output["BETA"] = np.select(conditions, choices)

    # Drop rows with NaN BETA values
    nrow = output.shape[0]
    output = output.dropna(subset=["BETA"])
    if output.shape[0] < nrow:
        print(
            f"Deleted {nrow-output.shape[0]} base SNPs that did not have matching alleles in reference data."
        )
    print(f"Found proxies for {output['proxy'].sum()} SNPs.")

    # Replace original SNPs with their proxy (if proxied)
    output["SNP"] = np.where(output["proxy"], output["SNP_A"], output["SNP"])
    output["POS"] = np.where(output["proxy"], output["BP_A"], output["POS"])
    output["CHR"] = np.where(output["proxy"], output["CHR_A"], output["CHR"])
    output["EA"] = np.where(output["proxy"], output["A1"], output["EA"])
    output["NEA"] = np.where(output["proxy"], output["A2"], output["NEA"])
    if "EAF" in output.columns:
        output["EAF"] = np.where(output["proxy"], output["MAF_A"], output["EAF"])

    # Drop columns related to ld
    output = output.drop(columns=ld.columns)

    return output


def apply_proxies(df, ld, searchspace=None):
    """
    Given a dataframe (coming from GENO.data attribute) and a dataframe of proxies
    (output from find_proxies), replace the SNPs in df with their best proxies, if they exist.
    This function is suited for exposure data (before running a PRS for instance).

    Args:
        df (DataFrame): Dataframe of SNP information with the usual GENO columns (SNP, BETA, SE, EAF, EA, NEA). EAF is not necessary.
        ld (DataFrame): Dataframe of proxies (output from find_proxies).
        searchspace (list, optional): List of SNPs to restrict the list of potential proxies. By default, includes all the proxies found. Using a searchspace can be done either at the find_proxies step or at this step, but it is much faster to use it at this step.

    Returns:
        DataFrame: A DataFrame with SNPs replaced by their best proxies, if they exist.
    """
    # Check mandatory columns
    mandatory_cols = ["EA", "SNP", "BETA"]
    for col in mandatory_cols:
        if col not in df.columns:
            raise ValueError(f"The column {col} is not found in the data!")
    
    # Filter by searchspace if provided
    if searchspace:
        print("Filtering the potential proxies with the searchspace provided.")
        ld = ld[ld.SNP_B.isin(searchspace)]

    # Remove original SNPs and sort by r2
    ld = ld[ld["SNP_A"] != ld["SNP_B"]]
    ld = ld.reindex(ld["R"].abs().sort_values(ascending=False).index)

    # Select the best proxy for each SNP
    ld = ld.groupby("SNP_A").first().reset_index(drop=False)

    # Merge the dataframes
    output = df.merge(ld, how="left", left_on="SNP", right_on="SNP_A")
    output["proxy"] = pd.notnull(output["SNP_B"])

    # Flip BETA if the original SNP alleles are switched in the reference panel
    conditions = [
        output["EA"] == output["A2"],
        output["EA"] == output["A1"],
        ~output["proxy"],
        ~output["EA"].isin([output["A1"], output["A2"]]) & output["proxy"],
    ]
    choices = [
        -output["BETA"],  # if EA == A2, flip the sign of BETA
        output["BETA"],  # if EA == A1, BETA does not change
        output["BETA"],  # if SNP_A is NaN (The original SNP was not proxied), BETA does not change
        np.nan,  # if the original SNP was proxied but EA is neither "A1" nor "A2", BETA is NaN
    ]
    output["BETA"] = np.select(conditions, choices)

    # Delete SNPs with mismatched alleles
    nrow = output.shape[0]
    output = output.dropna(subset=["BETA"])
    if output.shape[0] < nrow:
        print(
            f"Deleted {nrow-output.shape[0]} base SNPs that did not have matching alleles in reference data."
        )
    print(f"Found proxies for {output['proxy'].sum()} missing SNPs.")

    # Replace the original SNPs with their proxy (if proxied)
    output["SNP"] = np.where(output["proxy"], output["SNP_B"], output["SNP"])
    output["EA"] = np.where(output["proxy"], output["B1"], output["EA"])
    if "POS" in output.columns:
        output["POS"] = np.where(output["proxy"], output["BP_B"], output["POS"])
    if "CHR" in output.columns:
        output["CHR"] = np.where(output["proxy"], output["CHR_B"], output["CHR"])
    if "NEA" in output.columns:
        output["NEA"] = np.where(output["proxy"], output["B2"], output["NEA"])
    if "EAF" in output.columns:
        output["EAF"] = np.where(output["proxy"], output["MAF_B"], output["EAF"])

    # Drop ld columns
    output = output.drop(columns=ld.columns)

    return output


def find_proxies(
    snp_list,
    searchspace=None,
    reference_panel="eur",
    kb=5000,
    r2=0.6,
    window_snps=5000,
    threads=1,
):
    """
    Given a list of SNPs, return a table of proxies.

    Args:
        snp_list (list): List of rsids.
        searchspace (list, optional): List of SNPs to include in the search. By default, includes the whole reference panel.
        reference_panel (str, optional): The reference population to get linkage disequilibrium values and find proxies.
                                         Accepts values: "EUR", "SAS", "AFR", "EAS", "AMR".
                                         Alternatively, provide a path leading to a specific bed/bim/fam reference panel.
        kb (int, optional): Width of the genomic window to look for proxies. Defaults to 5000.
        r2 (float, optional): Minimum linkage disequilibrium value with the main SNP for a proxy to be included. Defaults to 0.6.
        window_snps (int, optional): Compute the LD value for SNPs that are not more than x SNPs apart from the main SNP. Defaults to 5000.
        threads (int, optional): Number of threads to use. Defaults to 1.

    Returns:
        DataFrame: A DataFrame containing the proxies. Only biallelic SNPs are returned.

    """
    # Ensure tmp_GENAL directory exists
    os.makedirs(f"tmp_GENAL/", exist_ok=True)

    # Convert snp_list to numpy array
    snp_list = np.array(list(snp_list))

    # Check if searchspace is provided
    if searchspace is None:
        extract_arg = ""
    else:
        print("Searching proxies in the provided searchspace.")
        with open("tmp_GENAL/searchspace.txt", "w") as file:
            for s in searchspace + snp_list:
                file.write(str(s) + "\n")
        extract_arg = "--extract tmp_GENAL/searchspace.txt"

    # Save snp_list to a file
    np.savetxt("tmp_GENAL/snps_to_proxy.txt", snp_list, fmt="%s", delimiter=" ")

    # Construct and execute the plink command
    command = f"{get_plink19_path()} --bfile {get_reference_panel_path(reference_panel)} {extract_arg} --keep-allele-order --r in-phase with-freqs gz --ld-snp-list tmp_GENAL/snps_to_proxy.txt --ld-window-kb {kb} --ld-window-r2 {r2} --ld-window {window_snps} --out tmp_GENAL/proxy.targets --threads {threads}"
    subprocess.run(command, shell=True, capture_output=True, text=True, check=True)

    # Read and process the output
    cmd = f"gunzip -c tmp_GENAL/proxy.targets.ld.gz"
    unzipped_content = subprocess.check_output(cmd, shell=True).decode("utf-8")
    ld = pd.read_csv(StringIO(unzipped_content), sep="\s+")

    # Filter out multiallelic SNPs
    ld["PHASE"] = ld["PHASE"].str.replace("/", "")
    ld = ld[ld["PHASE"].apply(len) == 4]
    ld = ld.reset_index(drop=True)

    # Split the "PHASE" column into separate characters
    temp = pd.DataFrame(
        ld["PHASE"].apply(list).to_list(), columns=["A1", "B1", "A2", "B2"]
    )
    ld = pd.concat([ld, temp], axis=1)

    # Convert integer columns to Int64 type
    for int_col in ["CHR_A", "CHR_B", "BP_A", "BP_B"]:
        ld[int_col] = ld[int_col].astype("Int64")

    return ld
