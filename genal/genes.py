import pandas as pd
import numpy as np
import os

from .constants import GENE_INFO_URL, CHR_CODES
from .tools import download_file, reference_data_dir



def filter_by_gene_func(data, gene_identifier, id_type="symbol", window_size=1000000, *, build):
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
            - "NCBI": NCBI gene ID (e.g., "348" or 348)
            - "UCSC": UCSC gene ID (e.g., "uc001hbu.2")
            - "Vega": Vega gene ID (e.g., "OTTHUMG00000019505")
            Default is "symbol". Text identifiers are matched case-insensitively.
        window_size (int, optional): Size of the window around the gene in base pairs. Default is 1,000,000 (1Mb).
        build (str): Genome build of the CHR/POS coordinates of the data, "37" or "38" (required, keyword only).

    Returns:
        pd.DataFrame: Filtered DataFrame containing only variants within the specified window
            around the gene, with additional column 'Distance'.

    Raises:
        ValueError: If the gene is not found, or if the gene information file has no chromosome or no coordinates
            in the requested build for this gene.

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
    build = str(build).strip()
    if build not in ["37", "38"]:
        raise ValueError(f"Invalid build: {build!r}. Must be one of: 37, 38")
    
    # Download the gene info file if not already present in the reference data folder
    data_dir = reference_data_dir()
    gene_info_file = os.path.join(data_dir, "gene_id_mapping_filtered.parquet")
    if not os.path.exists(gene_info_file):
        # Download parquet file
        print(f"Downloading gene info file to {gene_info_file}...")
        os.makedirs(data_dir, exist_ok=True)
        try:
            download_file(GENE_INFO_URL, gene_info_file, description="gene info parquet")
        except Exception as e:
            if os.path.exists(gene_info_file):
                os.remove(gene_info_file)
            raise RuntimeError(f"Failed to download gene info: {e}")

    df_gene_info = pd.read_parquet(gene_info_file, engine="pyarrow")
    
    # Find gene coordinates: NCBI ids are compared as integers, text identifiers case-insensitively
    if id_type == "NCBI_id":
        try:
            ncbi_id = int(str(gene_identifier).strip())
        except ValueError:
            raise ValueError(f"Invalid NCBI gene ID: {gene_identifier!r}. NCBI gene IDs are integers (e.g. 348).")
        matches = pd.to_numeric(df_gene_info[id_type], errors="coerce") == ncbi_id
    else:
        matches = (
            df_gene_info[id_type].astype("string").str.upper() == str(gene_identifier).strip().upper()
        ).fillna(False).astype(bool)
    gene_data = df_gene_info[matches]

    if gene_data.empty:
        raise ValueError(f"Gene with {id_type}='{gene_identifier}' not found in gene info database.")

    if len(gene_data) > 1:
        print(f"Warning: Multiple entries found for {id_type}='{gene_identifier}'. Using the first entry.")
    gene_data = gene_data.iloc[0,:]

    print(f"Filtering variants within {window_size}bp window (+/- {window_size/2}bp on each side) based on genome build {build} around gene: {', '.join(f'{col}: {gene_data[col]}' for col in valid_id_types)}")

    # Extract gene location information
    chrom = gene_data['CHR']
    if pd.isna(chrom):
        raise ValueError(
            f"Gene {gene_identifier} has no chromosome in the gene information file and cannot be used to filter variants."
        )
    # Convert to the numeric chromosome code (X=23, Y=24, XY=25, MT=26)
    chrom = str(chrom).strip().upper().removeprefix("CHR")
    if chrom.isdigit():
        chrom = int(chrom)
    elif chrom in CHR_CODES:
        chrom = CHR_CODES[chrom]
    else:
        raise ValueError(f"Gene {gene_identifier} is located on chromosome {chrom}, which is not supported.")

    if pd.isna(gene_data[f'gene_start_{build}']) or pd.isna(gene_data[f'gene_end_{build}']):
        raise ValueError(
            f"Gene {gene_identifier} has no coordinates in build {build} in the gene information file."
        )
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
        # Computed by position (not index labels, which may be repeated)
        pos = filtered['POS'].to_numpy(dtype="int64")
        distance = np.where(pos < gene_start, pos - gene_start, np.where(pos > gene_end, pos - gene_end, 0))
        filtered["Distance"] = pd.array(distance, dtype="Int64")
        
        print(f"Found {len(filtered)} variants.")
    else:
        print(f"No variants found in a {window_size}bp window around {gene_identifier}")
    
    return filtered
