import pandas as pd
import numpy as np
import scipy.stats as st
import os
import warnings
from collections import Counter

from .constants import STANDARD_COLUMNS
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


def check_snp_column(data):
    """Remove duplicates in the SNP column."""
    duplicate_indices = data[data.duplicated(subset=["SNP"], keep="first")].index
    n_del = len(duplicate_indices)
    n_initial = data.shape[0]
    if n_del > 0:
        data.drop(index=duplicate_indices, inplace=True)
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
    # If SE is missing
    if ("P" in data.columns) & ("BETA" in data.columns) & ("SE" not in data.columns):
        data["SE"] = np.where(
            data["P"] < 1, np.abs(data.BETA / st.norm.ppf(data.P / 2)), 0
        )
        print("The SE (Standard Error) column has been created.")
    # If P is missing
    if ("SE" in data.columns) & ("BETA" in data.columns) & ("P" not in data.columns):
        data["P"] = np.where(
            data["SE"] > 0, 2 * (1 - st.norm.cdf(np.abs(data.BETA) / data.SE)), 1
        )
        print("The P (P-value) column has been created.")
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


def check_beta_column(data, effect_column, preprocessing):
    """
    If the BETA column is a column of odds ratios, log-transform it.
    If no effect_column argument is specified, determine if the BETA column are beta estimates or odds ratios.
    """
    if effect_column is None:
        if preprocessing == 'None':
            return data
        median = np.median(data.BETA)
        has_negative = (data.BETA < 0).any()
        
        # Odds Ratios cannot be negative. If they are, it's a Beta.
        if 0.5 < median < 1.5 and not has_negative:
            effect_column = "OR"
            print(
                "The BETA column looks like Odds Ratios. Use effect_column='BETA' if it is a column of Beta estimates."
            )
        else:
            effect_column = "BETA"
            print(
                "The BETA column looks like Beta estimates. Use effect_column='OR' if it is a column of Odds Ratios."
            )

    ## Log transform the effect column if appropriate
    if effect_column not in ["BETA", "OR"]:
        raise ValueError(
            "The argument effect_column accepts only 'BETA' or 'OR' as values."
        )
    if effect_column == "OR":
        data["BETA"] = np.log(data["BETA"].clip(lower=0.01))
        data.drop(columns="SE", errors="ignore", inplace=True)
        print("The BETA column has been log-transformed to obtain Beta estimates.")
    return


def fill_ea_nea(data, reference_panel_df):
    """Fill in the EA and NEA columns based on reference data."""
    if "BETA" in data.columns:
        warnings.warn(
            "Warning: You have specified an effect (BETA) column but no effect allele (EA) column. An effect estimate is only meaningful if paired with the corresponding effect allele."
        )
    data = data.merge(
        reference_panel_df[["CHR", "POS", "A1", "A2"]], on=["CHR", "POS"], how="left"
    )
    n_missing = data["A1"].isna().sum()
    data.rename(columns={"A1": "EA", "A2": "NEA"}, inplace=True)

    perc_missing = n_missing / data.shape[0] * 100
    print(
        f"Alleles columns created: effect (EA) and non-effect allele (NEA). {n_missing}({perc_missing:.3f}%) values are set to nan because SNPs were not found in the reference data."
    )
    if perc_missing > 50:
        warnings.warn(
            f"The EA (Effect Allele) and NEA (Non-Effect Allele) for many SNPs could not been found. Make sure the CHR/POS coordinates are in the same build as the reference panel (GRCh37 (hg19) for the default one). If not, you can first use the .lift() method to lift them. For instance: .lift(start='hg38', end='hg19', replace=True) if they are in build GRCh38 (hg38)."
        )
    return data


def fill_nea(data, reference_panel_df):
    """Fill in the NEA column based on reference data."""
    data = data.merge(
        reference_panel_df[["CHR", "POS", "A1", "A2"]], on=["CHR", "POS"], how="left"
    )
    conditions = [data["EA"] == data["A1"], data["EA"] == data["A2"]]
    choices = [data["A2"], data["A1"]]
    data["NEA"] = np.select(conditions, choices, default=np.nan)
    n_missing = data["NEA"].isna().sum()
    data.drop(columns=["A1", "A2"], inplace=True)

    perc_missing = n_missing / data.shape[0] * 100
    print(
        f"The NEA (Non Effect Allele) column has been created. {n_missing}({perc_missing:.3f}%) values are set to nan because SNPs were not found in the reference data."
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
    data = data.merge(reference_panel_df[["CHR", "POS", "SNP"]], on="SNP", how="left")
    n_missing = data["CHR"].isna().sum()
    data["CHR"] = data["CHR"].astype("Int64")
    data["POS"] = data["POS"].astype("Int64")
    print(
        f"The coordinates columns (CHR for chromosome and POS for position) have been created. {n_missing}({n_missing/data.shape[0]*100:.3f}%) SNPs were not found in the reference data and their values set to nan."
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
        
        # Create allele_min and allele_max columns for data
        data['allele_min'] = np.minimum(data['EA'].fillna(""), data['NEA'].fillna(""))
        data['allele_max'] = np.maximum(data['EA'].fillna(""), data['NEA'].fillna(""))
        
        # Create or use allele_min and allele_max columns for reference panel
        if not ('allele_min' in reference_panel_df.columns and 'allele_max' in reference_panel_df.columns):
            reference_panel_df['allele_min'] = np.minimum(reference_panel_df['A1'].fillna(""), reference_panel_df['A2'].fillna(""))
            reference_panel_df['allele_max'] = np.maximum(reference_panel_df['A1'].fillna(""), reference_panel_df['A2'].fillna(""))
        
        # Merge using position and alleles
        data = pd.merge(
            data,
            reference_panel_df[["CHR", "POS", "SNP", "allele_min", "allele_max"]],
            on=["CHR", "POS", "allele_min", "allele_max"],
            how="left"
        )
        
        # Clean up temporary columns
        data.drop(columns=["allele_min", "allele_max"], inplace=True)
        
    else:
        # Default behavior - merge only on position
        data = pd.merge(
            data,
            reference_panel_df[["CHR", "POS", "SNP"]],
            on=["CHR", "POS"],
            how="left"
        )

    n_missing = data["SNP"].isna().sum()
    
    # Create standard names for missing SNPs if possible
    standard_name_condition = "EA" in data.columns and "NEA" in data.columns and n_missing > 0
    if standard_name_condition:
        missing_snp_condition = data["SNP"].isna()
        data.loc[missing_snp_condition, "SNP"] = (
            data.loc[missing_snp_condition, "CHR"].astype(str)
            + ":"
            + data.loc[missing_snp_condition, "POS"].astype(str)
            + ":"
            + data.loc[missing_snp_condition, "NEA"].astype(str)
            + ":"
            + data.loc[missing_snp_condition, "EA"].astype(str)
        )
        print_statement = f" and their ID set to CHR:POS:NEA:EA"
    
    perc_missing = n_missing / data.shape[0] * 100
    
    if n_missing == 0:
        print(
            f"The SNP column (rsID) has been created. All SNPs were found in the reference data."
        )
    else:
        print(
            f"The SNP column (rsID) has been created. {n_missing}({perc_missing:.3f}%) SNPs were not found in the reference data{print_statement if standard_name_condition else ''}."
        )
        
    if perc_missing > 50:
        warnings.warn(
            f"The SNPid for many SNPs could not been found. Make sure you are using a reference panel in the same genome build as your data. The one used by default is GRCh37 (hg19). \n"
            f"You can use the GRCh38 (hg38) reference panel by setting reference_panel = 38."
        )

    return data

def check_int_column(data, int_col):
    """Set the type of the int_col column to Int64 and non-numeric values to NA. This function is used to check the validity of the CHR and POS columns."""
    nrows = data.shape[0]
    # Remove any non-digit characters, convert to numeric, setting non-numeric to NaN
    data[int_col] = pd.to_numeric(data[int_col].astype(str).str.extract('(\d+)', expand=False), errors='coerce')
    # Convert to Int64 which handles NaN values, using round() first to handle floats
    data[int_col] = data[int_col].round().astype('Int64')
    n_nan = data[int_col].isna().sum()
    if n_nan > 0:
        print(
            f"The {int_col} column contains {n_nan}({n_nan/nrows*100:.3f}%) values set to NaN (due to being missing or non-integer)."
        )
    return

def adjust_column_names(data, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns):
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


def save_data(data, name, path="", fmt="h5", sep="\t", header=True):
    """
    Save data to a specified file format.

    Supported formats: .h5 (default), .csv, .txt.
    Future supported formats: .vcf, .vcf.gz.

    Args:
        data (pd.DataFrame): DataFrame to be saved.
        name (str): A unique identifier for the data, used as the filename.
        path (str, optional): The directory where the file will be saved. Defaults to current directory.
        fmt (str, optional): The desired file format. Defaults to "h5".
        sep (str, optional): Delimiter for text-based formats (.csv, .txt). Defaults to tab.
        header (bool, optional): Whether to include column names in text-based formats. Defaults to True.
    """
    path = os.path.join(path, name)
    if fmt == "h5":
        data.to_hdf(f"{path}.h5", key="data", mode="w", format="table")
    elif fmt == "csv":
        data.to_csv(f"{path}.csv", sep=sep, header=header, index=False)
    elif fmt == "txt":
        data.to_csv(f"{path}.txt", sep=sep, header=header, index=False)
    else:
        print(f"Format {fmt} is not supported yet.")


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
    """
    ref_panel_path, ref_filetype = get_reference_panel_path(reference_panel)
    create_tmp()

    by_coordinate = "CHR" in data.columns and "POS" in data.columns

    # --- Match by CHR/POS or SNP ID ---
    if by_coordinate:
        print("CHR/POS columns present. SNPs searched based on genomic positions.")
        
        # 1. Write coordinates to a temp file for PLINK's --extract range
        coord_path = os.path.join("tmp_GENAL", f"{object_name}_coord_list.txt")
        data[['CHR', 'POS', 'POS']].dropna().to_csv(coord_path, sep='\t', index=False, header=False)

        # 2. Run --freq directly, extracting by range and adding POS to output
        freq_prefix = os.path.join("tmp_GENAL", f"{object_name}_eaf_freqs")
        plink_command = (
            f"{get_plink_path()} --{'pfile' if ref_filetype == 'pgen' else 'bfile'} {ref_panel_path} "
            f"--memory {ram} "
            f"--extract range {coord_path} "
            f"--freq cols=+pos "
            f"--out {freq_prefix}"
        )
        run_plink_command(plink_command)

        # 3. Load frequency results
        freq_path = f"{freq_prefix}.afreq"
        if not os.path.exists(freq_path) or os.path.getsize(freq_path) == 0:
            warnings.warn("No variants from your data were found in the reference panel by coordinate, or PLINK failed.")
            return data.copy()
        
        freqs_df = pd.read_csv(freq_path, sep='\t')
        freqs_df.rename(columns={'#CHROM': 'CHR', 'ALT': 'ALT_calc', 'ALT_FREQS': 'EAF_ref'}, inplace=True)
        
        # 4. Merge with original data to get EA and compute final EAF
        data = data.merge(freqs_df[["CHR", "POS", "ALT_calc", "EAF_ref"]], on=['CHR', 'POS'], how='left')
        
    else:
        print("Using SNP IDs to extract frequencies.")
        if "SNP" not in data.columns:
            raise ValueError("SNP column is required when CHR/POS are not available.")
                    
        snp_list_path = os.path.join("tmp_GENAL", f"{object_name}_snp_list.txt")
        data[["SNP"]].dropna().to_csv(snp_list_path, index=False, header=False)

        freq_prefix = os.path.join("tmp_GENAL", f"{object_name}_eaf_freqs")
        plink_command = (
            f"{get_plink_path()} --{'pfile' if ref_filetype == 'pgen' else 'bfile'} {ref_panel_path} "
            f"--memory {ram} "
            f"--extract {snp_list_path} "
            f"--freq "
            f"--out {freq_prefix}"
        )
        run_plink_command(plink_command)

        freq_path = f"{freq_prefix}.afreq"
        if not os.path.exists(freq_path):
            warnings.warn("PLINK did not generate a frequency file. Cannot update EAF.")
            return data.copy()

        freqs_df = pd.read_csv(freq_path, sep='\t')
        freqs_df.rename(columns={"#CHROM": "CHR", "ID": "SNP", "ALT": "ALT_calc", "ALT_FREQS": "EAF_ref"}, inplace=True)

        data = data.merge(freqs_df[["SNP", "ALT_calc", "EAF_ref"]], on="SNP", how="left")

    if data["EAF_ref"].isna().all():
        warnings.warn("No matching SNPs found in the reference panel.")
        return data.copy()

    # Handle allele direction to ensure correct EAF is returned
    conditions = [
        data["EA"] == data["ALT_calc"],
        data["NEA"] == data["ALT_calc"],
    ]
    choices = [
        data["EAF_ref"],
        1 - data["EAF_ref"],
    ]
    data["EAF_new"] = np.select(conditions, choices, default=np.nan)

    # Create updated EAF column
    if 'EAF' not in data.columns:
        data['EAF'] = np.nan
    if fill:
        data['EAF'] = np.where(pd.notna(data["EAF_new"]), data["EAF_new"], data['EAF'])
    else:
        data['EAF'] = np.where(pd.notna(data["EAF_new"]), data["EAF_new"], np.nan)
    data.drop(columns=["EAF_new", "EAF_ref", "ALT_calc"], inplace=True)

    return data