import pandas as pd
import numpy as np
import os
import re
import uuid

from .tools import get_reference_panel_path, get_plink_path, load_reference_panel, run_plink_command
from .geno_tools import missing_snp_ids

## TO DO: accept lists of CHR/POS instead of SNP names for these functions

_COORD_TOKEN_RE = re.compile(
    r"^(chr)?(?P<chr>[0-9]{1,2}|X|Y|XY|MT|M)[:_](?P<pos>[0-9]+)([:_](?P<a1>[A-Za-z]+)[:_](?P<a2>[A-Za-z]+))?$",
    re.IGNORECASE,
)

_CHR_ALIASES = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}


def numeric_effect_columns(df):
    """Return a copy of df with the BETA, SE and EAF columns (when present) converted to numbers; values that cannot be converted become NaN."""
    df = df.copy()
    for col in ["BETA", "SE", "EAF"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _drop_palindromic_proxies(ld):
    """Remove candidate proxies with palindromic alleles (A/T or C/G): their strand cannot be resolved from the alleles."""
    pair = ld["MAJ_B"].astype(str).str.upper() + ld["NONMAJ_B"].astype(str).str.upper()
    palindromic = pair.isin(["AT", "TA", "CG", "GC"])
    if palindromic.any():
        print(f"Excluded {int(palindromic.sum())} palindromic (A/T or C/G) candidate proxies.")
    return ld[~palindromic]


def _drop_missing_proxy_ids(ld):
    """Remove the LD pairs whose target (SNP_A) or proxy (SNP_B) has a missing ID (NaN or placeholder such as '.')."""
    missing = missing_snp_ids(ld["SNP_A"]) | missing_snp_ids(ld["SNP_B"])
    if missing.any():
        print(f"Excluded {int(missing.sum())} candidate proxy pairs with a missing SNP ID.")
    return ld[~missing]


def _best_proxies(candidates):
    """
    Assign proxies to target SNPs (SNP_A) greedily by decreasing |R|: each target gets its best available candidate and a
    proxy serves at most one target (a target whose best proxy is taken by a target in higher LD falls back to its next best).
    """
    order = candidates["R"].abs().sort_values(ascending=False, kind="mergesort").index
    used_targets, used_proxies, keep = set(), set(), []
    for idx, target, proxy in zip(order, candidates.loc[order, "SNP_A"], candidates.loc[order, "SNP_B"]):
        if target in used_targets or proxy in used_proxies:
            continue
        used_targets.add(target)
        used_proxies.add(proxy)
        keep.append(idx)
    return candidates.loc[keep]


def query_outcome_proxy(df, ld, snps_to_extract, snps_df=[]):
    """
    Extract the best proxies from a dataframe, as well as specific SNPs.

    Given a dataframe `df` (originating from Geno.data) and a dataframe of potential proxies
    (output from `find_proxies`), this function extracts the best proxies from `df` as well as
    the SNPs specified in `snps_to_extract`.
    This is suited for querying outcome data.

    Candidate proxies must be present in `df`, differ from their target SNP, not be part of
    `snps_to_extract` (those SNPs are used as themselves), have non-missing IDs, be non-palindromic, have alleles
    in `df` matching the reference panel alleles, and have valid BETA/SE values (finite, SE > 0). Proxies are then
    assigned greedily by decreasing |r|: each target gets its best remaining candidate and a proxy serves at most
    one target. The proxy association is reported for the target SNP, with the target's
    minor allele in the reference panel as effect allele. If `df` has an EAF column, the EAF of the proxied
    SNP is the reference panel frequency of this effect allele (the target's own frequency, used by the
    palindrome checks of the harmonization).

    Args:
        df (pd.DataFrame): Dataframe of SNP information with the usual Geno columns
                           (SNP, BETA, SE, EAF, EA, NEA). EAF is not necessary.
        ld (pd.DataFrame): Dataframe of proxies (output from `find_proxies`).
        snps_to_extract (list): List of SNPs to extract in addition to the proxies.
        snps_df (list, optional): List of SNPs to choose the proxy from. Should be the list of
                                  SNPs in df. Can be provided to avoid recomputing it. Defaults to an empty list.

    Returns:
        pd.DataFrame: Dataframe with queried SNPs and their proxies (the "proxy" column indicates proxied SNPs).
    """
    # If ld is None
    if not isinstance(ld, pd.DataFrame):
        raise ValueError("ld is None (The SNPs to be proxied were not found in the reference panel)")
    df = numeric_effect_columns(df)

    # If snps_df is empty, populate it with SNPs from df
    if len(snps_df) == 0:
        snps_df = df.SNP.values

    # Candidate proxies: present in df, different from the target and not extracted as themselves
    ld = ld[
        ld["SNP_B"].isin(snps_df)
        & (ld["SNP_A"] != ld["SNP_B"])
        & ~ld["SNP_B"].isin(snps_to_extract)
    ]
    ld = _drop_palindromic_proxies(_drop_missing_proxy_ids(ld))

    # Pair each candidate with its association in df and keep the allele-compatible pairs
    # (alleles in df equal to the reference panel alleles of the proxy, in either order)
    candidates = df.merge(ld, how="inner", left_on="SNP", right_on="SNP_B")
    ea = candidates["EA"].astype(str).str.upper()
    nea = candidates["NEA"].astype(str).str.upper()
    compatible = ((ea == candidates["NONMAJ_B"]) & (nea == candidates["MAJ_B"])) | (
        (ea == candidates["MAJ_B"]) & (nea == candidates["NONMAJ_B"])
    )
    if (~compatible).any():
        print(
            f"Excluded {int((~compatible).sum())} candidate proxies whose alleles do not match the reference panel alleles."
        )
    # Candidates with a missing or invalid association (non-finite BETA/SE or SE <= 0) cannot be used
    beta = pd.to_numeric(candidates["BETA"], errors="coerce")
    se = pd.to_numeric(candidates["SE"], errors="coerce")
    valid_effect = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if (compatible & ~valid_effect).any():
        print(
            f"Excluded {int((compatible & ~valid_effect).sum())} candidate proxies with missing or invalid BETA/SE values in the outcome data."
        )
    proxies = _best_proxies(candidates[compatible & valid_effect]).copy()

    # In the plink output, the alleles taken as reference for proxying are "MAJ_A" and "MAJ_B" (major alleles in the reference panel)
    # We want to use as effect allele for the original SNP its minor allele in the reference panel (NONMAJ_A)
    # So, we flip BETA if the proxy's effect allele is the major allele in the reference panel (BETA of NONMAJ_B)...
    ea_is_maj = proxies["EA"].astype(str).str.upper() == proxies["MAJ_B"]
    beta_nonmaj_b = np.where(ea_is_maj, -proxies["BETA"], proxies["BETA"])
    # ...and flip it again if R is negative: NONMAJ_A is then in phase with MAJ_B
    negative_r = np.sign(proxies["R"]) == -1
    proxies["BETA"] = np.where(negative_r, -beta_nonmaj_b, beta_nonmaj_b)

    # EAF of the proxied SNP: reference panel frequency of its effect allele NONMAJ_A
    if "EAF" in proxies.columns:
        proxies["EAF"] = proxies["NONMAJ_FREQ_A"]

    # Replace the proxied SNPs with the position and alleles of the original SNPs
    proxies["SNP"] = proxies["SNP_A"]
    proxies["EA"] = proxies["NONMAJ_A"]
    proxies["NEA"] = proxies["MAJ_A"]
    if "POS" in proxies.columns:
        proxies["POS"] = proxies["BP_A"]
    if "CHR" in proxies.columns:
        proxies["CHR"] = proxies["CHR_A"]
    print(f"Found proxies for {len(proxies)} SNPs.")

    output = pd.concat(
        [
            df[df.SNP.isin(snps_to_extract)].assign(proxy=False),
            proxies[df.columns].assign(proxy=True),
        ],
        ignore_index=True,
    )

    return output


def apply_proxies(df, ld, searchspace=None):
    """
    Given a dataframe (coming from GENO.data attribute) and a dataframe of proxies
    (output from find_proxies), replace the SNPs in df with their best proxies, if they exist.
    This function is suited for exposure data (before running a PRS for instance).

    Candidate proxies must differ from their target SNP, not already be a SNP of df, have non-missing IDs and be
    non-palindromic. The target SNP alleles in df (EA, and NEA if present) must match the reference panel alleles; SNPs
    with candidate proxies but mismatching alleles are deleted. Proxies are assigned greedily by decreasing |r|: each
    target gets its best remaining candidate and a proxy serves at most one target.

    Args:
        df (DataFrame): Dataframe of SNP information with the usual GENO columns (SNP, BETA, SE, EAF, EA, NEA). EAF is not necessary.
        ld (DataFrame): Dataframe of proxies (output from find_proxies).
        searchspace (list, optional): List of SNPs to restrict the list of potential proxies. By default, includes all the proxies found. Using a searchspace can be done either at the find_proxies step or at this step, but it is much faster to use it at this step.

    Returns:
        DataFrame: A DataFrame with SNPs replaced by their best proxies, if they exist.
    """
    # If ld is None
    if not isinstance(ld, pd.DataFrame):
        raise ValueError("ld is None (The SNPs to be proxied were not found in the reference panel)")
    
    # Check mandatory columns
    mandatory_cols = ["EA", "SNP", "BETA"]
    for col in mandatory_cols:
        if col not in df.columns:
            raise ValueError(f"The column {col} is not found in the data!")
    df = numeric_effect_columns(df)
    
    # Filter by searchspace if provided
    if searchspace:
        print("Filtering the potential proxies with the searchspace provided.")
        ld = ld[ld.SNP_B.isin(searchspace)]

    # Candidate proxies: different from the target and not already in df (those SNPs are used as themselves)
    ld = ld[(ld["SNP_A"] != ld["SNP_B"]) & ~ld["SNP_B"].isin(df["SNP"])]
    ld = _drop_palindromic_proxies(_drop_missing_proxy_ids(ld))

    # Keep the candidates of target SNPs whose alleles match the reference panel alleles of the target
    targets = df[["SNP", "EA"] + (["NEA"] if "NEA" in df.columns else [])].drop_duplicates(subset="SNP")
    candidates = targets.merge(ld, how="inner", left_on="SNP", right_on="SNP_A")
    ea = candidates["EA"].astype(str).str.upper()
    compatible = (ea == candidates["MAJ_A"]) | (ea == candidates["NONMAJ_A"])
    if "NEA" in candidates.columns:
        nea = candidates["NEA"].astype(str).str.upper()
        compatible &= ((ea == candidates["MAJ_A"]) & (nea == candidates["NONMAJ_A"])) | (
            (ea == candidates["NONMAJ_A"]) & (nea == candidates["MAJ_A"])
        )
    mismatched = set(candidates.loc[~compatible, "SNP_A"]) - set(candidates.loc[compatible, "SNP_A"])
    ld = _best_proxies(candidates.loc[compatible, ld.columns])

    # Merge the dataframes
    output = df[~df["SNP"].isin(mismatched)].merge(ld, how="left", left_on="SNP", right_on="SNP_A")
    output["proxy"] = pd.notnull(output["SNP_B"])
    if mismatched:
        print(
            f"Deleted {len(mismatched)} base SNPs that did not have matching alleles in reference data."
        )

    # In the plink output, the alleles taken as reference for proxying are "MAJ_A" and "MAJ_B" (major alleles in the reference panel)
    # We want to use as effect allele for the proxy SNP its minor allele in the reference panel
    # So, we flip BETA if the original SNP's effect allele is the major allele in the reference panel
    ea_is_maj = output["proxy"] & (output["EA"].astype(str).str.upper() == output["MAJ_A"])
    output["BETA"] = np.where(ea_is_maj, -output["BETA"], output["BETA"])

    # Flip BETA if the sign of R is negative: indicates that the positive correlation corresponds to MAJ_A with NONMAJ_B
    sign_r = np.sign(output["R"]) # Sign of R
    output["BETA"] = np.where(sign_r == -1, -output["BETA"], output["BETA"])
    print(f"Found proxies for {output['proxy'].sum()} SNPs.")

    # Replace the original SNPs with their proxy (if proxied)
    # As said above, we use as effect allele the minor allele in the reference panel
    output["SNP"] = np.where(output["proxy"], output["SNP_B"], output["SNP"])
    output["EA"] = np.where(output["proxy"], output["NONMAJ_B"], output["EA"])
    if "POS" in output.columns:
        output["POS"] = np.where(output["proxy"], output["BP_B"], output["POS"])
    if "CHR" in output.columns:
        output["CHR"] = np.where(output["proxy"], output["CHR_B"], output["CHR"])
    if "NEA" in output.columns:
        output["NEA"] = np.where(output["proxy"], output["MAJ_B"], output["NEA"])
    if "EAF" in output.columns:
        output["EAF"] = np.where(output["proxy"], output["NONMAJ_FREQ_B"], output["EAF"])

    # Drop ld columns
    output.drop(columns=ld.columns, inplace=True)

    return output


def _match_reference_ids(coord, ref, keys):
    """
    Map coordinate tokens to reference panel IDs on the columns `keys` (CHR/POS, optionally with the allele pair).
    Tokens matching several reference variants (different IDs) are left unresolved, so the result does not depend on
    the order of the reference panel. Returns the tokens with a SNP column (NaN if unresolved) and the set of the
    _order values of the ambiguous tokens.
    """
    ref = ref[keys + ["SNP"]].drop_duplicates()
    n_matches = ref.groupby(keys).size().rename("n_matches").reset_index()
    unique_ref = ref.merge(n_matches[n_matches["n_matches"] == 1][keys], on=keys)
    mapped = coord.merge(unique_ref, on=keys, how="left")
    ambiguous = coord.merge(n_matches[n_matches["n_matches"] > 1][keys], on=keys)["_order"]
    return mapped, set(ambiguous)


def _normalize_chr_token(chr_token):
    chr_token = str(chr_token).upper()
    if chr_token.isdigit():
        chr_int = int(chr_token)
        if 1 <= chr_int <= 26:
            return chr_int
        return None
    return _CHR_ALIASES.get(chr_token)


def _parse_coord_token(token):
    token = str(token)
    match = _COORD_TOKEN_RE.match(token)
    if not match:
        return None

    chr_int = _normalize_chr_token(match.group("chr"))
    if chr_int is None:
        return None

    pos_int = int(match.group("pos"))
    a1 = match.group("a1")
    a2 = match.group("a2")

    allele_min = None
    allele_max = None
    if a1 and a2:
        a1 = str(a1).upper()
        a2 = str(a2).upper()
        allele_min = min(a1, a2)
        allele_max = max(a1, a2)

    return {
        "raw": token,
        "CHR": chr_int,
        "POS": pos_int,
        "allele_min": allele_min,
        "allele_max": allele_max,
    }


def find_proxies(
    snp_list,
    searchspace=None,
    reference_panel="EUR_37",
    kb=5000,
    r2=0.8,
    window_snps=1000000,
    threads=1,
    name=None
):
    """
    Given a list of SNPs, return a table of proxies using PLINK 2.0.

    Args:
        snp_list (list): List of variant identifiers. Can be rsids (e.g. "rs123") or
            coordinate tokens like "1:12345", "chr1:12345", optionally with alleles
            ("1:12345:A:G"). Coordinate tokens are resolved to the reference panel ID
            namespace using CHR/POS (and alleles when provided). A token matching several reference
            variants (e.g. a multi-allelic position given without alleles) is not resolved.
        searchspace (list, optional): List of SNPs to include in the search. By default, includes the whole reference panel.
        reference_panel (str, optional): The reference population to get linkage disequilibrium values and find proxies.
            Acceptable populations are "EUR", "SAS", "AFR", "EAS", "AMR" and available builds are 37 and 38 ("EUR_38" or "AFR_37" etc...)
            Also accepts or a path to a specific bed/bim/fam or pgen/pvar/psam panel.
            Default is "EUR_37".
        kb (int, optional): Width of the genomic window to look for proxies. Defaults to 5000.
        r2 (float, optional): Minimum linkage disequilibrium value with the main SNP for a proxy to be included. Defaults to 0.8.
        window_snps (int, optional): Compute the LD value for SNPs that are not more than x SNPs apart from the main SNP. Defaults to 1000000 (equivalent to infinity).
        threads (int, optional): Number of threads to use. Defaults to 1.

    Returns:
        DataFrame: A DataFrame containing the proxies. Only biallelic SNPs are returned.
    """
    # Ensure tmp_GENAL directory exists
    os.makedirs(f"tmp_GENAL/", exist_ok=True)

    # Generate a default name if none is provided
    if name is None:
        name = str(uuid.uuid4())[:8]

    # Normalize input tokens to strings and parse coordinate tokens (CHR:POS[:A1:A2]).
    tokens = [str(s) for s in snp_list]
    id_tokens = []
    coord_records = []
    for order, token in enumerate(tokens):
        parsed = _parse_coord_token(token)
        if parsed is None:
            id_tokens.append(token)
        else:
            parsed["_order"] = order
            coord_records.append(parsed)

    # Get reference panel path and type (needed for both mapping and LD computation).
    ref_path, filetype = get_reference_panel_path(reference_panel)

    # Resolve coordinate tokens to the reference panel ID namespace using existing helpers.
    id_to_raw = {}
    resolved_ids = []
    if coord_records:
        ref_df = load_reference_panel(ref_path)
        coord_df = pd.DataFrame(coord_records)

        mapped_frames = []
        ambiguous_orders = set()
        needs_alleles = coord_df["allele_min"].notna().any()
        needs_pos = coord_df["allele_min"].isna().any()
        if needs_alleles:
            a1 = ref_df["A1"].fillna("").astype(str)
            a2 = ref_df["A2"].fillna("").astype(str)
            ref_allele = pd.DataFrame(
                {
                    "CHR": ref_df["CHR"],
                    "POS": ref_df["POS"],
                    "allele_min": np.minimum(a1, a2),
                    "allele_max": np.maximum(a1, a2),
                    "SNP": ref_df["SNP"],
                }
            )

            coord_with = coord_df[coord_df["allele_min"].notna()]
            if not coord_with.empty:
                mapped, ambiguous = _match_reference_ids(coord_with, ref_allele, ["CHR", "POS", "allele_min", "allele_max"])
                mapped_frames.append(mapped)
                ambiguous_orders |= ambiguous
            coord_without = coord_df[coord_df["allele_min"].isna()]
        else:
            coord_without = coord_df

        if needs_pos and not coord_without.empty:
            mapped, ambiguous = _match_reference_ids(coord_without, ref_df, ["CHR", "POS"])
            mapped_frames.append(mapped)
            ambiguous_orders |= ambiguous

        mapped_all = (
            pd.concat(mapped_frames, ignore_index=True)
            .sort_values("_order", kind="mergesort")
            .reset_index(drop=True)
        )
        if ambiguous_orders:
            print(
                f"{len(ambiguous_orders)} coordinate-based variants match several reference panel variants and were not resolved: "
                "give them as CHR:POS:A1:A2 or as rsIDs."
            )
        n_dropped = int(mapped_all["SNP"].isna().sum()) - len(ambiguous_orders)
        if n_dropped > 0:
            print(
                f"Dropped {n_dropped} coordinate-based variants not found in the reference panel."
            )

        mapped_resolved = mapped_all.dropna(subset=["SNP"])
        if not mapped_resolved.empty:
            id_to_raw = (
                mapped_resolved.drop_duplicates(subset=["SNP"], keep="first")
                .set_index("SNP")["raw"]
                .to_dict()
            )
            resolved_ids = mapped_resolved["SNP"].astype(str).tolist()

    # Final targets for PLINK LD computation (rsids + resolved coordinate-based IDs).
    targets_for_plink = id_tokens + resolved_ids
    if not targets_for_plink:
        return None

    # Check if searchspace is provided
    if searchspace is None:
        extract_arg = []
    else:
        print("Searching proxies in the provided searchspace.")
        with open(f"tmp_GENAL/{name}_searchspace.txt", "w") as file:
            for s in searchspace:
                file.write(str(s) + "\n")
            for s in targets_for_plink:
                file.write(str(s) + "\n")
        extract_arg = ["--extract", os.path.join("tmp_GENAL", f"{name}_searchspace.txt")]

    # Save targets to a file
    np.savetxt(
        f"tmp_GENAL/{name}_snps_to_proxy.txt",
        np.array(targets_for_plink, dtype=str),
        fmt="%s",
        delimiter=" ",
    )

    # Construct and execute the plink2 command (input flag based on filetype)
    command = [
        get_plink_path(),
        "--bfile" if filetype == "bed" else "--pfile", ref_path,
        *extract_arg,
        "--r-unphased", "cols=chrom,pos,id,maj,nonmaj,freq",
        "--ld-snp-list", os.path.join("tmp_GENAL", f"{name}_snps_to_proxy.txt"),
        "--ld-window-kb", str(kb),
        "--ld-window-r2", str(r2),
        "--ld-window", str(window_snps),
        "--threads", str(threads),
        "--out", os.path.join("tmp_GENAL", f"{name}_proxy.targets"),
    ]

    run_plink_command(command)

    # Read log file to return amount of SNPs to be proxied present in the ref panel
    log_path = os.path.join("tmp_GENAL", f"{name}_proxy.targets.log")
    log_content = open(log_path).read()
    match = re.search(r"(\d+) variants? remaining", log_content)
    if match:
        n_present = int(match.group(1))
        if n_present == 0:
            print("None of the SNPs to be proxied are present in the reference panel.")
            return None
        else:
            print(f"{n_present} SNPs to be proxied are present in the reference panel.")

    # Read and process the output
    try:
        ld = pd.read_csv(f"tmp_GENAL/{name}_proxy.targets.vcor", sep=r"\s+")
    except FileNotFoundError:
        print("No proxies found that meet the specified criteria.")
        return None
      
    # Rename columns to match the expected format
    ld.rename(columns={
        'ID_A': 'SNP_A',
        'ID_B': 'SNP_B',
        '#CHROM_A': 'CHR_A',
        'CHROM_B': 'CHR_B',
        'POS_A': 'BP_A',
        'POS_B': 'BP_B',
        'UNPHASED_R': 'R',
    }, inplace=True)

    # Create PHASE column for compatibility
    #ld['PHASE'] = ld['A1'] + ld['B1'] + ld['A2'] + ld['B2']

    # Filter out multiallelic SNPs
    #ld = ld[ld["PHASE"].str.len() == 4]
    #ld = ld.reset_index(drop=True)

    # Convert integer columns to Int64 type
    for int_col in ["CHR_A", "CHR_B", "BP_A", "BP_B"]:
        ld[int_col] = ld[int_col].astype("Int64")

    # Rewrite coordinate-based targets back to the caller namespace for downstream merges.
    if id_to_raw:
        ld["SNP_A"] = ld["SNP_A"].map(id_to_raw).fillna(ld["SNP_A"])

    return ld
