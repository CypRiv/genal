import pandas as pd
import numpy as np
from pandas.api.types import is_numeric_dtype
import scipy.stats as st
import os, subprocess

from .extract_prs import check_bfiles
from .tools import get_plink19_path


def association_test_func(data, covar_list, standardize, name, data_pheno, pheno_type):
    """
    Conduct single-SNP association tests against a phenotype.

    This function performs a series of operations:
        1. Checks for necessary preliminary steps.
        2. Updates the FAM file with the phenotype data.
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
    if not check_bfiles(genetic_path):
        raise FileNotFoundError(
            "Run the extract_snps() method before performing association tests."
        )
    if data.shape[0] == 0:
        raise ValueError(
            "No SNPs for the association tests. Check the .data or .data_clumped dataframes."
        )

    # Update phenotype in the FAM file
    fam = _prepare_fam_file(genetic_path, data_pheno, pheno_type, standardize)

    # Prepare covariate file if covariates are provided
    covar, covar_filename = _handle_covariates(covar_list, data_pheno, name)

    # Execute PLINK association test
    output, method = _run_plink_assoc_test(
        genetic_path, name, pheno_type, covar, covar_filename, covar_list
    )

    # Process and return results
    return _process_results(output, method, data, pheno_type)


def _prepare_fam_file(genetic_path, data_pheno, pheno_type, standardize):
    """Helper function to prepare the FAM file with phenotype data."""
    # Read the FAM file
    fam = pd.read_csv(genetic_path + ".fam", header=None, delimiter=" ")
    # Extract relevant phenotype data
    data_pheno_trait = data_pheno[["IID", "PHENO"]].rename(columns={"IID": 0}).copy()
    # Merge phenotype data with the FAM dataframe
    fam = fam.merge(data_pheno_trait, how="left", on=0, indicator=True)
    fam[5] = fam.PHENO
    # Verify that the merge was successful
    if (fam["_merge"] == "both").sum() == 0:
        raise ValueError(
            "The IDs in the phenotype dataframe are inconsistent with those in the genetic dataset. Call set_phenotype() method again, specifying the correct column name for the genetic IDs."
        )
    fam.drop(axis=1, columns=["PHENO", "_merge"], inplace=True)
    # Count the number of individuals with a valid phenotype trait
    n_non_na = fam.shape[0] - fam[5].isna().sum()
    print(
        f"{n_non_na} individuals are present in the genetic data and have a valid phenotype trait."
    )
    # Update phenotype values based on its type
    if pheno_type == "binary":
        fam[5] = fam[5] + 1
        fam[5] = fam[5].astype("Int64")
    if (pheno_type == "quant") & (standardize == True):
        # Standardizing for quantitative phenotypes
        print(
            "Standardizing the phenotype to approximate a normal distribution. Use standardize = False if you do not want to standardize."
        )
        fam[5] = (fam[5] - fam[5].mean()) / fam[5].std()
    fam[5] = fam[5].fillna(-9)
    fam.to_csv(genetic_path + ".fam", header=None, index=False, sep=" ")
    return fam


def _handle_covariates(covar_list, data_pheno, name):
    """Helper function to prepare the covariate file."""
    if len(covar_list) > 0:
        # Ensure all covariates are present in phenotype data
        for col in covar_list:
            if col not in data_pheno.columns:
                raise TypeError(
                    f"The {col} column is not found in the .phenotype dataframe."
                )
        # Select required columns and rename columns
        data_cov = data_pheno[["IID", "IID"] + covar_list].copy()
        data_cov.columns = ["FID"] + list(data_cov.columns[1:])
        # Remove rows with NA values and print their number
        nrows = data_cov.shape[0]
        data_cov.dropna(inplace=True)
        removed_rows = nrows - data_cov.shape[0]
        if removed_rows > 0:
            print(
                f"{removed_rows}({removed_rows/nrows*100:.3f}%) individuals have NA values in the covariates columns and will be excluded from the association tests."
            )
        # Ensure the covariates are numeric and not trivial (lead to association fail)
        for col in covar_list:
            if data_pheno[col].nunique() == 1:
                print(
                    f"The {col} covariate contains only one value and is removed from the tests."
                )
                data_cov.drop(axis=1, columns=[col], inplace=True)
            if not pd.api.types.is_numeric_dtype(data_pheno[col]):
                print(
                    f"The {col} covariate is not numeric and is removed from the tests."
                )
                data_cov.drop(axis=1, columns=[col], inplace=True, errors="ignore")
        # Define the covariate filename
        covar_filename = os.path.join("tmp_GENAL", f"{name}_covar.cov")
        data_cov.to_csv(covar_filename, sep=" ", header=True, index=False)
        covar = True
    else:
        covar = False
        covar_filename = None
    return covar, covar_filename


def _run_plink_assoc_test(
    genetic_path, name, pheno_type, covar, covar_filename, covar_list
):
    """Helper function to execute the PLINK association test."""
    method = "logistic" if pheno_type == "binary" else "linear"
    print(
        f"Running single-SNP {method} regression tests on {genetic_path} data {f'with adjustment for: {covar_list}' if covar else 'without covariate adjustments'}."
    )
    # Formulate the covariate argument for the PLINK command
    covar_argument = f"--covar {covar_filename} --hide-covar" if covar == True else ""
    output = os.path.join("tmp_GENAL", name)
    # Construct and run the PLINK command
    command = f"{get_plink19_path()} --bfile {genetic_path} --{method} {covar_argument} --out {output}"
    subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
    return output, method


def _process_results(output, method, data, pheno_type):
    """Helper function to process results after the PLINK association test."""
    # Path to PLINK results
    results_path = output + f".assoc." + method
    assoc = pd.read_csv(results_path, delimiter="\s+")

    # If logistic regression, log-transform the odds ratio
    assoc["BETA"] = np.log(assoc.OR) if pheno_type == "binary" else assoc.BETA
    
    n_na = assoc["BETA"].isna().sum()
    
    # Merge results with the clumped data
    data = data.drop(axis=1, columns=["BETA", "CHR", "P"], errors="ignore").merge(
        assoc, how="inner", on="SNP"
    )

    # Adjust beta values based on allele match
    data["BETA"] = np.where(
        data.EA == data.A1, data.BETA, np.where(data.NEA == data.A1, -data.BETA, np.nan)
    )

    # Calculate and set standard error values
    data["SE"] = np.abs(data.BETA / st.norm.ppf(data.P / 2))
    # Drop unnecessary columns
    data = data.drop(
        axis=1, columns=["A1", "TEST", "NMISS", "OR", "STAT", "BP"], errors="ignore"
    )

    # Remove rows with mismatches in allele columns and notify the user
    nrow_previous = data.shape[0]
    data = data.dropna(subset="BETA")
    delta_nrow = nrow_previous - data.shape[0] - n_na
    if (delta_nrow > 0) or (n_na > 0):
        print(
            f"{f'{n_na}({n_na/nrow_previous*100:.3f}%) SNP-trait tests returned NA value and ' if n_na>0 else ''}{delta_nrow}({delta_nrow/nrow_previous*100:.3f}%) SNPs removed due to allele discrepancies between the main data and the genetic data."
        )
    return data


def set_phenotype_func(data_original, PHENO, PHENO_type, IID, alternate_control):
    """
    Set a phenotype dataframe containing individual IDs and phenotype columns formatted for single-SNP association testing.

    Args:
        data (pd.DataFrame): Contains at least an individual IDs column and one phenotype column.
        IID (str): Name of the individual IDs column in data.
        PHENO (str): Name of the phenotype column in data.
        PHENO_type (str, optional): Type of the phenotype column. Either "quant" for quantitative (continuous) or "binary".The function tries to infer the type if not provided.
        alternate_control (bool): Assumes that for a binary trait, the controls are coded with the most frequent value. Use True to reverse the assumption.

    Returns:
        pd.DataFrame: The modified data.
        str: The inferred or provided PHENO_type.

    Raises:
        ValueError: For inconsistencies in the provided data or arguments.

    This function corresponds to the following GENO method: :meth:`GENO.set_phenotype`.
    """
    data = data_original.copy()
    _validate_columns_existence(data, PHENO, IID)
    data = _standardize_column_names(data, PHENO, IID)
    PHENO_type = _determine_phenotype_type(data, PHENO_type)
    data = _validate_and_process_phenotype(data, PHENO, PHENO_type, alternate_control)
    _report_na_values(data)

    print("The phenotype data is stored in the .phenotype attribute.")
    return data, PHENO_type


def _validate_columns_existence(data, PHENO, IID):
    """Checks if columns exist and raises errors if not."""
    for column in [PHENO, IID]:
        # Raise an error if the column name is not provided
        if column is None:
            raise ValueError(f"Please provide a name for the {column} variable.")
        # Raise an error if the column does not exist in the data
        if column not in data.columns:
            raise ValueError(
                f"The column '{column}' is not present in the dataset. This column is required!"
            )
    if data.shape[0] == 0:
        raise ValueError("The phenotype dataframe is empty.")


def _standardize_column_names(data, PHENO, IID):
    """Standardizes the column names to 'IID' and 'PHENO'."""
    # Drop redundant columns if they exist and rename the target columns to standard names
    if PHENO != "PHENO":
        data.drop(axis=1, columns=["PHENO"], errors="ignore", inplace=True)
    if IID != "IID":
        data.drop(axis=1, columns=["IID"], errors="ignore", inplace=True)
    data.rename(columns={IID: "IID", PHENO: "PHENO"}, inplace=True)
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
    else:
        if not PHENO_type in ["binary", "quant"]:
            raise ValueError(
                f"The only possible values for the PHENO_type argument are 'binary' or 'quant'"
            )
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
    """Processes a binary phenotype."""
    # Ensure that the phenotype is binary
    if len(np.unique(data.PHENO.dropna())) != 2:
        raise ValueError(
            f"The '{PHENO}' column is not binary as it contains more than two distinct values."
        )

    if alternate_control:
        code_control = data.PHENO.value_counts().index[1]
        code_case = data.PHENO.value_counts().index[0]
    else:
        code_control = data.PHENO.value_counts().index[0]
        code_case = data.PHENO.value_counts().index[1]
    
    print(
        f"Identified {code_control} as the control code in 'PHENO'. {'Set alternate_control=True to inverse this interpretation.' if not alternate_control else ''}"
    )

    # Update the control and case codings
    data.replace({"PHENO": {code_control: 0, code_case: 1}}, inplace=True)


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
