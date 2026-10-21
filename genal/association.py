import pandas as pd
import numpy as np
from pandas.api.types import is_bool_dtype, is_numeric_dtype
import scipy.stats as st
import os, glob, re

from .extract_prs import check_pfiles
from .tools import get_plink_path, run_plink_command, COMPLEMENT

PSAM_SAMPLE_COLUMNS = ["#FID", "FID", "#IID", "IID", "SID", "PAT", "MAT", "SEX"]


def association_test_func_plink2(data, covar_list, standardize, name, data_pheno, pheno_type):
    """
    Conduct single-SNP association tests against a phenotype.

    This function performs a series of operations:
        1. Checks for necessary preliminary steps.
        2. Updates the PSAM file with the phenotype data.
        3. Creates a covariate file if required.
        4. Runs a PLINK association test.
        5. Processes the results and returns them.

    Args:
        data (pd.DataFrame): Genetic data with the standard Geno columns.
        covar_list (list): List of column names in the data_pheno DataFrame to use as covariates.
        standardize (bool): Flag indicating if the phenotype needs standardization.
        name (str): Prefix for the filenames used during the process.
        data_pheno (pd.DataFrame): Phenotype data with at least an IID and PHENO columns.
        pheno_type (str): Type of phenotype ('binary' or 'quant').

    Returns:
        pd.DataFrame: Processed results of the association test.

    This function corresponds to the following Geno method: :meth:`Geno.association_test`.
    """

    # Check necessary files are available
    genetic_path = os.path.join("tmp_GENAL", f"{name}_allchr")
    if not check_pfiles(genetic_path):
        raise FileNotFoundError(
            "Run the extract_snps() method before performing association tests."
        )
    if data.shape[0] == 0:
        raise ValueError(
            "No SNPs for the association tests. Check the .data or .data_clumped dataframes."
        )

    # Each individual is used once
    duplicated = data_pheno.duplicated(subset=["FID", "IID"], keep="first") & data_pheno[["FID", "IID"]].notna().all(axis=1)
    if duplicated.sum() > 0:
        print(f"{duplicated.sum()} rows of the phenotype data have a duplicated FID/IID and are removed (the first row of each individual is kept).")
        data_pheno = data_pheno[~duplicated]

    # Update phenotype in the PSAM file
    psam = _prepare_psam_file(genetic_path, data_pheno, pheno_type, standardize)

    # Prepare the covariate file and the file of the analysed samples if covariates are provided
    covar_list, covar_filename, keep_filename = _handle_covariates(covar_list, data_pheno, name, psam)

    # Execute PLINK association test
    output = _run_plink2_assoc_test(
        genetic_path, name, covar_filename, covar_list, pheno_type, keep_filename
    )

    # Process and return results
    return _process_results_plink2(output, data, pheno_type)

def _run_plink2_assoc_test(
    genetic_path, name, covar_filename, covar_list, pheno_type, keep_filename=None
):
    """Helper function to execute the PLINK 2.0 association test (on the samples of keep_filename when given)."""
    
    print(
        f"Running {'linear' if pheno_type == 'quant' else 'logistic'} association tests on {genetic_path} data "
        f"{f'with adjustment for: {covar_list}' if len(covar_list) > 0 else 'without covariates. This is not recommended'}."
    )
    
    output = os.path.join("tmp_GENAL", name)
    
    # Build the PLINK 2.0 command (the extracted genetic data is a pgen/pvar/psam fileset)
    command = [
        get_plink_path(),
        "--pfile", genetic_path,
        "--glm", 
        *(["allow-no-covars"] if len(covar_list) == 0 else []),
        "no-x-sex",
        "cols=+omitted",
        "--no-input-missing-phenotype",
        "--pheno-name", "PHENO1"
    ]
    
    # Covariates are variance-standardized for numerical stability (the SNP effect estimates are unchanged)
    if len(covar_list) > 0:
        command.extend([
            "--covar", covar_filename,
            "--covar-name", ",".join(covar_list),
            "--covar-variance-standardize",
        ])

    # Complete-case sample of the requested covariates (also when no covariate column remains)
    if keep_filename is not None:
        command.extend(["--keep", keep_filename])

    command.extend(["--out", output])
    
    run_plink_command(command)
    
    return output

def _process_results_plink2(output, data, pheno_type):
    """Helper function to process results after the PLINK association test."""
    # Path to PLINK results: <output>.PHENO1.glm.<regression type>
    results_files = glob.glob(output + ".PHENO1.glm.*")
    if len(results_files) != 1:
        raise FileNotFoundError(
            f"Expected one PLINK association results file '{output}.PHENO1.glm.*', found {len(results_files)}. Check the {output}.log file."
        )
    results_path = results_files[0]
    regression = results_path.split(".glm.")[-1]
    if (pheno_type == "quant") != (regression == "linear"):
        raise ValueError(
            f"PLINK ran a '{regression}' regression for a phenotype of type '{pheno_type}'. Check the {output}.log file."
        )
    assoc = pd.read_csv(results_path, sep="\t")

    # Filter to keep only coefficients corresponding to our phenotype
    assoc = assoc[assoc["TEST"] == "ADD"].copy()

    # If logistic regression, log-transform the odds ratio
    assoc["BETA"] = np.log(assoc.OR) if pheno_type == "binary" else assoc.BETA

    # Rename columns
    assoc.rename(
        columns={"ID": "SNP", "LOG(OR)_SE": "SE", "A1": "A1_tested", "OMITTED": "OMITTED_tested", "REF": "REF_tested", "ALT": "ALT_tested"},
        errors="ignore", inplace=True,
    )
    
    # Merge results with the data based on the variant IDs used for the extraction
    data = data.drop(columns=["BETA", "SE", "P"], errors="ignore").merge(
        assoc[["SNP", "BETA", "SE", "A1_tested", "OMITTED_tested", "REF_tested", "ALT_tested", "P"]], how="inner", on="SNP"
    )

    # Align beta values on EA: each PLINK result is the effect of the allele A1 versus the OMITTED allele
    nea = data["NEA"] if "NEA" in data.columns else pd.Series(np.nan, index=data.index)
    matches = [
        _effect_sign(ea, n, a1, omitted, ref, alt)
        for ea, n, a1, omitted, ref, alt in zip(
            data["EA"], nea, data["A1_tested"], data["OMITTED_tested"], data["REF_tested"], data["ALT_tested"]
        )
    ]
    tested = data["BETA"].notna()  # results that are not NA in PLINK
    status = pd.Series([match for _, match in matches], index=data.index)[tested]
    data["BETA"] = data["BETA"] * np.array([sign for sign, _ in matches], dtype=float)
    n_na = int((~tested).sum())
    n_opposite = int((status == "opposite_strand").sum())
    n_other_contrast = int((status == "other_contrast").sum())
    n_undefined_contrast = int((status == "undefined_contrast").sum())
    if n_opposite > 0:
        print(f"{n_opposite} SNPs match the alleles of the tested variant on the opposite strand and are aligned on their EA.")
    if n_other_contrast > 0:
        print(
            f"{n_other_contrast} results of multiallelic variants are removed: the contrast tested by PLINK (A1 versus the OMITTED allele) is not EA versus NEA."
        )
    if n_undefined_contrast > 0:
        print(
            f"{n_undefined_contrast} results of multiallelic variants are removed: without NEA, EA does not define the contrast tested by PLINK."
        )

    # Drop unnecessary columns
    data = data.drop(
        columns=["A1_tested", "OMITTED_tested", "REF_tested", "ALT_tested"], errors="ignore"
    )

    # Remove rows with mismatches in allele columns and notify the user
    nrow_previous = data.shape[0]
    data = data.dropna(subset="BETA")
    delta_nrow = int((status == "mismatch").sum())
    if (delta_nrow > 0) or (n_na > 0):
        print(
            f"{f'{n_na}({n_na/max(1, nrow_previous)*100:.3f}%) SNP-trait tests returned NA value and ' if n_na>0 else ''}{delta_nrow}({delta_nrow/max(1, nrow_previous)*100:.3f}%) SNPs removed due to allele discrepancies between the main data and the genetic data."
        )
    return data


def _effect_sign(ea, nea, a1, omitted, ref, alt):
    """
    Sign aligning on the effect allele EA the effect estimated by PLINK for the allele A1 versus the OMITTED allele.

    With NEA, {EA, NEA} must be {A1, OMITTED}, also accepted on the opposite strand for SNVs that are not palindromic:
    EA = A1 keeps the sign and EA = OMITTED flips it. At a multiallelic variant (alleles REF and ALT), PLINK fits one
    term per allele and each result is the contrast of one allele versus the omitted allele: only the result of the
    EA versus NEA contrast is kept. Without NEA, the variant must be biallelic, with EA = A1 or EA = OMITTED.

    Returns:
        tuple: Sign (1, -1, or NaN when the result is not the effect of EA versus NEA) and the match type:
            "direct", "opposite_strand", "other_contrast" (multiallelic variant, the result is another contrast),
            "undefined_contrast" (multiallelic variant without NEA) or "mismatch" (alleles do not match).
    """
    alleles = {str(ref).upper(), *str(alt).upper().split(",")}
    multiallelic = len(alleles) > 2
    a1, omitted, ea = str(a1).upper(), str(omitted).upper(), str(ea).upper()
    if pd.isna(nea):
        if multiallelic:
            return np.nan, "undefined_contrast"
        if ea == a1:
            return 1.0, "direct"
        return (-1.0, "direct") if ea == omitted else (np.nan, "mismatch")
    nea = str(nea).upper()

    def sign(effect, other):
        if (effect, other) == (a1, omitted):
            return 1.0
        return -1.0 if (effect, other) == (omitted, a1) else np.nan

    direct = sign(ea, nea)
    if not np.isnan(direct):
        return direct, "direct"
    if ea in COMPLEMENT and nea in COMPLEMENT and COMPLEMENT[ea] != nea:
        opposite = sign(COMPLEMENT[ea], COMPLEMENT[nea])
        if not np.isnan(opposite):
            return opposite, "opposite_strand"
    if multiallelic and ea in alleles and nea in alleles and ea != nea:
        return np.nan, "other_contrast"
    return np.nan, "mismatch"


def _merge_on_psam_ids(psam, data_pheno, columns):
    """
    Match the rows of the phenotype data to the samples of the PSAM file.

    IDs are compared as strings, on FID and IID when the PSAM file has a #FID column and on IID otherwise. The requested
    columns are merged under internal names, so that their names never collide with the merge bookkeeping.

    Returns:
        tuple: pd.DataFrame with the requested columns of data_pheno, one row per PSAM sample in the PSAM order (NA for
            samples without phenotype data), and np.ndarray of booleans (PSAM samples found in data_pheno).
    """
    psam_keys, pheno_keys = (["#FID", "IID"], ["FID", "IID"]) if "#FID" in psam.columns else (["#IID"], ["IID"])
    internal = [f"_column_{i}" for i in range(len(columns))]
    pheno = data_pheno.dropna(subset=pheno_keys)
    pheno = pd.concat(
        [pheno[pheno_keys].astype(str).set_axis(psam_keys, axis=1), pheno[columns].set_axis(internal, axis=1)], axis=1
    )
    merged = psam[psam_keys].astype(str).merge(
        pheno.drop_duplicates(subset=psam_keys), how="left", on=psam_keys, indicator=True
    )
    matched = (merged["_merge"] == "both").to_numpy()
    return merged[internal].set_axis(columns, axis=1), matched


def _prepare_psam_file(genetic_path, data_pheno, pheno_type, standardize):
    """Helper function to write the phenotype data in the PHENO1 column of the PSAM file."""
    # Read the PSAM file, keeping the IDs as strings
    psam = pd.read_csv(genetic_path + ".psam", sep="\t", dtype=str, keep_default_na=False)

    # The phenotype columns already present in the genetic data are not used
    psam_phenotypes = [col for col in psam.columns if col not in PSAM_SAMPLE_COLUMNS]
    if len(psam_phenotypes) > 0:
        print(f"The phenotype columns of the genetic data ({', '.join(psam_phenotypes)}) are ignored.")
        psam = psam.drop(columns=psam_phenotypes)

    # Merge phenotype data with the PSAM dataframe
    values, matched = _merge_on_psam_ids(psam, data_pheno, ["PHENO"])

    # Verify that the merge was successful
    if matched.sum() == 0:
        raise ValueError(
            "The IDs in the phenotype dataframe are inconsistent with those in the genetic dataset. Call set_phenotype() method again, specifying the correct column names for the genetic IDs (IID and FID)."
        )
    psam["PHENO1"] = values["PHENO"].to_numpy()
    
    # Count the number of individuals with a valid phenotype trait
    n_non_na = psam["PHENO1"].notna().sum()
    print(
        f"{n_non_na} individuals are present in the genetic data and have a valid phenotype trait."
    )
    
    # Update phenotype values based on its type
    if pheno_type == "binary":
        # PLINK case/control coding: 1 = control, 2 = case
        psam["PHENO1"] = (psam["PHENO1"] + 1).astype("Int64")
    elif standardize:
        # Standardizing for quantitative phenotypes
        print(
            "Standardizing the phenotype to approximate a normal distribution. Use standardize = False if you do not want to standardize."
        )
        psam["PHENO1"] = (psam["PHENO1"] - psam["PHENO1"].mean(skipna=True)) / psam["PHENO1"].std(skipna=True)

    # PLINK reads a phenotype taking only the values 0 (missing), 1 and 2 as case/control
    if pheno_type == "quant":
        values = psam["PHENO1"].dropna()
        if len(values) > 0 and values.isin([0, 1, 2]).all():
            psam["PHENO1"] = psam["PHENO1"] + 3
            print(
                "The phenotype only takes the values 0, 1 or 2, which PLINK reads as case/control codes: "
                "3 is added to all its values so that it is analysed as a quantitative trait (this does not change the SNP effect estimates)."
            )

    # Make sure the SEX column is not empty without modifying existing values
    if "SEX" in psam.columns:
        psam["SEX"] = psam["SEX"].replace("", "NA")

    # PLINK2 expects missing values to be represented as 'NA' tokens (not empty fields).
    psam.to_csv(genetic_path + ".psam", sep="\t", index=False, na_rep="NA")
    return psam


def _handle_covariates(covar_list, data_pheno, name, psam):
    """
    Helper function to prepare the covariate file and the file of the analysed samples (sample IDs of the PSAM file).

    As in a complete-case regression, the analysed individuals are the samples of the genetic data with a phenotype and
    non-missing values for all the covariates; they are written to a keep file (PLINK --keep), so that the exclusion
    holds even when no covariate column remains. A covariate with a single value among them is removed. Non-numeric
    covariates (text or categories) are encoded as k - 1 indicator (0/1) columns named '<covariate>=<level>', the first
    level in sorted order being the reference. Numeric covariates are used as continuous variables (numeric codes, e.g.
    of assessment centres, must be converted to str or category). Column names are made PLINK-compatible and unique.

    Returns:
        tuple: Names of the covariate columns of the covariate file, covariate filename (None without covariate column)
            and keep filename (([], None, None) without covariates).
    """
    if len(covar_list) == 0:
        return [], None, None

    # Ensure all covariates are present in phenotype data (PHENO, IID and FID hold the phenotype and the sample IDs)
    covar_list = list(dict.fromkeys(covar_list))
    for col in covar_list:
        if col in ["PHENO", "IID", "FID"]:
            raise ValueError(
                f"'{col}' can't be used as a covariate: after set_phenotype(), the PHENO, IID and FID columns hold the "
                "phenotype and the sample IDs. Rename the covariate column and call set_phenotype() again."
            )
        if col not in data_pheno.columns:
            raise TypeError(f"The {col} column is not found in the .phenotype dataframe.")

    # Covariate values of the samples of the genetic data, and analysed individuals (phenotype and all covariates)
    # Sample IDs and covariate values are kept in separate tables, so that covariate names never collide with the IDs
    covariates, _ = _merge_on_psam_ids(psam, data_pheno, covar_list)
    id_columns = ["#FID", "IID"] if "#FID" in psam.columns else ["#IID"]
    ids = psam[id_columns].astype(str).reset_index(drop=True)
    with_pheno = psam["PHENO1"].notna().to_numpy()
    complete = covariates.notna().all(axis=1).to_numpy()
    n_excluded = (with_pheno & ~complete).sum()
    if n_excluded > 0:
        print(
            f"{n_excluded}({n_excluded/max(1, with_pheno.sum())*100:.3f}%) individuals with a valid phenotype have NA values in the covariates columns and will be excluded from the association tests."
        )
    if not (with_pheno & complete).any():
        raise ValueError("No individual of the genetic data with a phenotype has non-missing values for all the covariates.")
    ids, covariates = ids[complete], covariates[complete]
    analysed = with_pheno[complete]

    # Covariate columns: numeric covariates, and k - 1 indicator columns for the k levels of the non-numeric ones
    covar_data = ids.copy()
    used_names = set(id_columns)
    for col in covar_list:
        values = covariates[col]
        if is_numeric_dtype(values):
            if values[analysed].nunique() < 2:
                print(f"The {col} covariate contains only one value and is removed from the tests.")
                continue
            covar_data[_covariate_name(col, used_names)] = values.astype(int) if is_bool_dtype(values) else values
            continue
        values = values.astype(str).str.strip()
        levels = sorted(values[analysed].unique())
        if len(levels) < 2:
            print(f"The {col} covariate contains only one value and is removed from the tests.")
            continue
        for level in levels[1:]:
            covar_data[_covariate_name(f"{col}={level}", used_names)] = (values == level).astype(int)
        shown = ", ".join(levels[:10]) + (", ..." if len(levels) > 10 else "")
        print(
            f"The {col} covariate is categorical: its {len(levels)} levels ({shown}) are encoded as {len(levels) - 1} "
            f"indicator columns, with '{levels[0]}' as the reference level."
        )
    covar_columns = [col for col in covar_data.columns if col not in id_columns]

    # Write the analysed samples and, if any covariate column remains, the covariate file
    keep_filename = os.path.join("tmp_GENAL", f"{name}_covar.keep")
    ids.to_csv(keep_filename, sep="\t", header=True, index=False)
    covar_filename = None
    if len(covar_columns) > 0:
        covar_filename = os.path.join("tmp_GENAL", f"{name}_covar.cov")
        covar_data.to_csv(covar_filename, sep="\t", header=True, index=False)
    return covar_columns, covar_filename, keep_filename


def _covariate_name(name, used_names):
    """Covariate column name accepted by PLINK (letters, digits, '_', '.', '='), made unique among used_names (which is updated)."""
    base = re.sub(r"[^0-9A-Za-z_.=]+", "_", str(name))
    unique, i = base, 2
    while unique in used_names:
        unique, i = f"{base}_{i}", i + 1
    used_names.add(unique)
    return unique


### __________________________
### Set phenotype functions
### __________________________

def set_phenotype_func(data_original, PHENO, PHENO_type, IID, FID=None, alternate_control=False):
    """
    Set a phenotype dataframe containing individual IDs and phenotype columns formatted for single-SNP association testing.

    Args:
        data (pd.DataFrame): Contains at least an individual IDs column and one phenotype column.
        IID (str): Name of the individual IDs column in data.
        PHENO (str): Name of the phenotype column in data.
        PHENO_type (str, optional): Type of the phenotype column. Either "quant" for quantitative (continuous) or "binary".
            The function tries to infer the type if not provided.
        FID (str, optional): Name of the family ID column in data. If not provided, FID will be set to IID values.
        alternate_control (bool): For a binary trait, values 0/1 are read as 0 = control, 1 = case, values 1/2 as
            1 = control, 2 = case, and for any other pair of values the most frequent value is the control.
            Use True to invert the control and case values.

    Returns:
        pd.DataFrame: The modified data.
        str: The inferred or provided PHENO_type.
    """
    data = data_original.copy()
    _validate_columns_existence(data, PHENO, IID, FID)
    
    data = _standardize_column_names(data, PHENO, IID, FID)
    PHENO_type = _determine_phenotype_type(data, PHENO_type)
    data = _validate_and_process_phenotype(data, PHENO, PHENO_type, alternate_control)
    _report_na_values(data)

    print("The phenotype data is stored in the .phenotype attribute.")
    return data, PHENO_type


def _validate_columns_existence(data, PHENO, IID, FID):
    """Checks if columns exist and raises errors if not."""
    # Check if PHENO is a string
    if not isinstance(PHENO, str):
        raise ValueError("The PHENO argument must be a string containing the name of the phenotype column.")
    # Check if IID is a string
    if not isinstance(IID, str):
        raise ValueError("The IID argument must be a string containing the name of the individual IDs column.")
    
    for column in [PHENO, IID]:
        # Raise an error if the column name is not provided
        if column is None:
            raise ValueError(f"Please provide a name for the {column} variable.")
        # Raise an error if the column does not exist in the data
        if column not in data.columns:
            raise ValueError(
                f"The column '{column}' is not present in the dataset. This column is required!"
            )
        
    # Handle FID column
    if FID is not None and FID not in data.columns:
        raise ValueError(f"The column '{FID}' is not present in the provided dataset.")

    # If IID or FID is numerical, convert to integer
    if is_numeric_dtype(data[IID]):
        data[IID] = data[IID].astype("Int64")
    if FID is not None and is_numeric_dtype(data[FID]):
        data[FID] = data[FID].astype("Int64")
    
    if data.shape[0] == 0:
        raise ValueError("The phenotype dataframe is empty.")


def _standardize_column_names(data, PHENO, IID, FID):
    """Standardizes the column names to 'IID' and 'PHENO'."""
    # Drop redundant columns if they exist and rename the target columns to standard names
    if PHENO != "PHENO":
        data.drop(columns=["PHENO"], errors="ignore", inplace=True)
    if IID != "IID":
        data.drop(columns=["IID"], errors="ignore", inplace=True)
    data.rename(columns={IID: "IID", PHENO: "PHENO"}, inplace=True)
    
    if FID is not None:
        if FID != "FID":
            data.drop(columns=["FID"], errors="ignore", inplace=True)
        data.rename(columns={FID: "FID"}, inplace=True)
    else:
        data["FID"] = data["IID"]
        print(
            "The FID column was not provided. The FIDs are assumed to be the same as the IIDs."
        )
    
    return data


def _determine_phenotype_type(data, PHENO_type):
    """Guesses or validates the phenotype type."""
    # If phenotype type is not given, deduce it based on the unique values in the column
    if PHENO_type is None:
        if len(np.unique(data.PHENO.dropna())) == 2:
            print(
                "Detected a binary phenotype in the 'PHENO' column. Specify 'PHENO_type=\"quant\"' if this is incorrect."
            )
            return "binary"
        else:
            print(
                "Detected a quantitative phenotype in the 'PHENO' column. Specify 'PHENO_type=\"binary\"' if this is incorrect."
            )
            return "quant"
    return PHENO_type


def _validate_and_process_phenotype(data, PHENO, PHENO_type, alternate_control):
    """Validates the phenotype and processes it accordingly."""
    # Process the phenotype based on its type
    if PHENO_type == "binary":
        _process_binary_phenotype(data, PHENO, alternate_control)
    elif PHENO_type == "quant":
        _validate_quantitative_phenotype(data, PHENO)
    else:
        raise ValueError("Accepted values for 'PHENO_type' are 'binary' or 'quant'.")
    return data


def _process_binary_phenotype(data, PHENO, alternate_control):
    """
    Processes a binary phenotype: recode controls as 0 and cases as 1.

    Values 0/1 are read as 0 = control, 1 = case; values 1/2 as 1 = control, 2 = case (PLINK coding);
    for any other pair of values, the most frequent value is the control. alternate_control=True inverts the control and case values.
    """
    # Ensure that the phenotype is binary
    values = np.unique(data.PHENO.dropna())
    if len(values) != 2:
        raise ValueError(
            f"The '{PHENO}' column is not binary as it does not contain exactly two distinct values."
        )

    if set(values) == {0, 1}:
        code_control, code_case = 0, 1
        reason = "0/1 coding: 0 = control, 1 = case"
    elif set(values) == {1, 2}:
        code_control, code_case = 1, 2
        reason = "1/2 coding: 1 = control, 2 = case"
    else:
        code_control, code_case = data.PHENO.value_counts().index[:2]
        reason = "the most frequent value is the control"
    if alternate_control:
        code_control, code_case = code_case, code_control
        reason = f"{reason}, inverted by alternate_control=True"
    
    print(
        f"Identified {code_control} as the control code and {code_case} as the case code in '{PHENO}' ({reason}). "
        f"{'Set alternate_control=True to inverse this interpretation.' if not alternate_control else ''}"
    )

    # Update the control and case codings
    data["PHENO"] = data["PHENO"].map({code_control: 0, code_case: 1})

    # Print number and percentage of cases
    n_case = int(data.PHENO.sum())
    print(
        f"There are {n_case}({n_case/data.shape[0]*100:.3f}%) cases in the 'PHENO' column."
    )


def _validate_quantitative_phenotype(data, PHENO):
    """Validates a quantitative phenotype."""
    # Ensure that the phenotype is numeric
    if not is_numeric_dtype(data.PHENO):
        raise ValueError(
            f"The '{PHENO}' column must contain numeric values for a quantitative phenotype."
        )


def _report_na_values(data):
    """Reports the number of NA values in 'IID' and 'PHENO' columns."""
    nrows = data.shape[0]
    n_nan_id = data.IID.isna().sum()
    n_nan_pheno = data.PHENO.isna().sum()

    # Report NA values in ID and PHENO columns, if they exist
    if n_nan_id > 0:
        print(
            f"Detected {n_nan_id} NA values in the 'ID' column, accounting for {n_nan_id/nrows*100:.3f}% of entries. These will be omitted during analyses."
        )
    if n_nan_pheno > 0:
        print(
            f"Detected {n_nan_pheno} NA values in the 'PHENO' column, accounting for {n_nan_pheno/nrows*100:.3f}% of entries. These will be omitted during analyses."
        )
