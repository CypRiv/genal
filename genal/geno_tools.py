import pandas as pd
import numpy as np
import scipy.stats as st
import os, subprocess
import shutil
import warnings
from collections import Counter
import wget

from .constants import STANDARD_COLUMNS, BUCKET_URL
from .tools import read_config



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
        if 0.5 < median < 1.5:
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
    """Set the type of the int_col column to Int64 and non-numeric values to NA."""
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
    Save a DataFrame to a file in a given format.

    Args:
    - data (pd.DataFrame): The data to be saved.
    - name (str): The name of the file without extension.
    - path (str, optional): Directory path for saving. Default is the current directory.
    - fmt (str, optional): Format for the file, e.g., "h5", "csv", "txt". Default is "h5".
    - sep (str, optional): Delimiter for csv or txt files. Default is tab.
    - header (bool, optional): Whether to include header in csv or txt files. Default is True.

    Returns:
    None. But saves the data to a file and prints the file path.

    Raises:
    - ValueError: If the provided format is not recognized.
    """
    if path:
        path_name = os.path.join(path, f"{name}.{fmt}")
    else:
        path_name = f"{name}.{fmt}"

    if fmt == "h5":
        df = data.copy()
        for col in df.select_dtypes(include="integer").columns:
            df[col] = df[col].astype("float64")
        df.to_hdf(path_name, mode="w", key="data")

    elif fmt in ["csv", "txt"]:
        data.to_csv(path_name, sep=sep, header=header, index=False)

    else:
        raise ValueError(
            "The fmt argument takes value in (h5 (default), csv, txt)."
        )

    print(f"Data saved to {path_name}")


def Combine_Geno(Gs):
    """
    Combine a list of GWAS objects into one.

    Args:
    - Gs (list): List of GWAS objects.

    Returns:
    Geno object: Combined Geno object.
    """
    from .Geno import Geno
    
    C = pd.DataFrame()

    for G in Gs:
        C = pd.concat([C, G.data])

    C = C.reset_index(drop=True)

    return Geno(C)

def filter_by_gene_func(data, gene_identifier, id_type="symbol", window_size=1000000, build="37"):
    """
    Filtering the data to include only variants that are within a specified distance of a specific gene.
    Corresponds to the :meth:`Geno.filter_by_gene` method.
    Args:
        data (pd.DataFrame): Input data with at least 'CHR' and 'POS' columns.
        gene_identifier (str): Identifier for the gene/protein to filter variants around.
        id_type (str, optional): Type of identifier provided. Options are:
            - "symbol": Gene symbol (e.g., "APOE")
            - "HGNC": HGNC ID (e.g., "HGNC:613")
            - "name": Full gene name (e.g., "apolipoprotein E")
            - "Ensembl": Ensembl gene ID (e.g., "ENSG00000130203")
            - "NCBI": NCBI gene ID (e.g., "348")
            - "UCSC": UCSC gene ID (e.g., "uc001hbu.2")
            - "Vega": Vega gene ID (e.g., "OTTHUMG00000019505")
            Default is "symbol".
        window_size (int, optional): Size of the window around the gene in base pairs. Default is 1,000,000 (1Mb).
        build (str, optional): Genome build of the data. Default is "37".
        
    Returns:
        pd.DataFrame: Filtered DataFrame containing only variants within the specified window 
            around the gene, with additional column 'Distance'.

    Notes:
        - Distance is calculated from the nearest gene boundary (start or end position)
        - Null distances indicate the variant is within the gene
    """
        
    # Validate id_type
    valid_id_types = ["symbol", "HGNC_id", "name", "gene_id", "NCBI_id", "UCSC_id", "Vega_id"]
    if id_type in ["HGNC", "NCBI", "UCSC", "Vega"]:
        id_type = id_type + "_id"
    if id_type == "Ensembl":
        id_type = "gene_id"
    if id_type not in valid_id_types:
        raise ValueError(f"Invalid id_type. Must be one of: {', '.join(valid_id_types)}")
    
    # Validate build
    if int(build) not in [37, 38]:
        raise ValueError(f"Invalid build. Must be one of: 37, 38")
    
    # Download the gene info file if not already present in the reference folder
    config = read_config()
    ref_path = config["paths"]["ref_path"]
    gene_info_file = os.path.join(ref_path, "gene_id_mapping_filtered.parquet")
    if not os.path.exists(gene_info_file):
        # Download parquet file
        print(f"Downloading gene info file to {gene_info_file}...")    
        url = BUCKET_URL + "gene_id_mapping_filtered.parquet"
        try:
            wget.download(url, gene_info_file)
            print("\nDownload complete.")
        except Exception as e:
            if os.path.exists(gene_info_file):
                os.remove(gene_info_file)
            raise RuntimeError(f"Failed to download gene info: {e}")

    df_gene_info = pd.read_parquet(gene_info_file, engine="pyarrow")
    
    # Find gene coordinates
    gene_data = df_gene_info[df_gene_info[id_type] == gene_identifier]
    
    if gene_data.empty:
        raise ValueError(f"Gene with {id_type}='{gene_identifier}' not found in gene info database.")
    
    if len(gene_data) > 1:
        print(f"Warning: Multiple entries found for {id_type}='{gene_identifier}'. Using the first entry.")
    gene_data = gene_data.iloc[0,:]

    print(f"Filtering variants within {window_size}bp window based on genome build {build} around gene: {', '.join(f'{col}: {gene_data[col]}' for col in valid_id_types)}")
    
    # Extract gene location information
    chrom = gene_data['CHR']
    # Convert to integer if possible
    if str(chrom).isdigit():
        chrom = int(chrom)
    elif chrom=="X":
        chrom=23
    else:
        raise ValueError(f"Gene {gene_identifier} is located on chromosome {chrom}, which is not supported.")
    
    gene_start = int(gene_data[f'gene_start_{build}'])
    gene_end = int(gene_data[f'gene_end_{build}'])

    # Define the window boundaries
    window_start = max(0, gene_start - window_size/2)
    window_end = gene_end + window_size/2
    
    # Filter variants within the window
    filtered = data[
        (data['CHR'] == chrom) & 
        (data['POS'] >= window_start) & 
        (data['POS'] <= window_end)
    ].copy()

    if not filtered.empty:
        # Calculate distance from gene: if inside the gene, distance is 0, if before, distance is negative, if after, distance is positive
        filtered.loc[:, 'Distance'] = np.nan
        
        # Create boolean masks
        mask_inside = filtered['POS'].between(gene_start, gene_end)
        mask_before = filtered['POS'] < gene_start
        mask_after  = filtered['POS'] > gene_end

        filtered.loc[mask_inside, 'Distance'] = 0
        filtered.loc[mask_before, 'Distance'] = filtered['POS'] - gene_start
        filtered.loc[mask_after, 'Distance']  = filtered['POS'] - gene_end

        filtered["Distance"] = filtered["Distance"].astype("Int64")
        
        print(f"Found {len(filtered)} variants.")
    else:
        print(f"No variants found in a {window_size}bp window around {gene_identifier}")
    
    return filtered