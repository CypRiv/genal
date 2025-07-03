import pandas as pd
import numpy as np
import os
import wget

from .constants import BUCKET_URL
from .tools import read_config



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