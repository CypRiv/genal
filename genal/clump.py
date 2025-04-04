import os
import pandas as pd
import uuid
import re

from .tools import get_reference_panel_path, get_plink_path, run_plink_command

def clump_data_plink2(
    data,
    reference_panel="EUR_37",
    kb=10000,
    r2=0.01,
    p1=5e-8,
    p2=0.01,
    name=None,
    ram=10000,
):
    """
    Perform clumping on the given data using plink. Corresponds to the :meth:`Geno.clump` method.

    Args:
        data (pd.DataFrame): Input data with at least 'SNP' and 'P' columns.
        reference_panel (str, optional): The reference population to get linkage disequilibrium values and find proxies.
            Acceptable populations are "EUR", "SAS", "AFR", "EAS", "AMR" and available builds are 37 and 38 ("EUR_38" or "AFR_37" etc...)
            Also accepts or a path to a specific bed/bim/fam or pgen/pvar/psam panel.
            Default is "EUR_37".
        kb (int, optional): Clumping window in terms of thousands of SNPs. Default is 10000.
        r2 (float, optional): Linkage disequilibrium threshold, values between 0 and 1. Default is 0.01.
        p1 (float, optional): P-value threshold during clumping. SNPs above this value are not considered. Default is 5e-8.
        p2 (float, optional): P-value threshold post-clumping to further filter the clumped SNPs. If p2 < p1, it won't be considered. Default is 0.01.
        name (str, optional): Name used for the files created in the tmp_GENAL folder.
        ram (int, optional): Amount of RAM in MB to be used by plink.

    Returns:
        pd.DataFrame: Data after clumping, if any.
    """

    # Create unique ID for the name if none is passed
    if name is None:
        name = str(uuid.uuid4())[:8]

    # Save the relevant data columns to a temporary file
    to_clump_filename = os.path.join("tmp_GENAL", f"{name}_to_clump.txt")
    data[["SNP", "P"]].to_csv(to_clump_filename, index=False, sep="\t")

    # Get reference panel path and type
    ref_path, filetype = get_reference_panel_path(reference_panel)

    # Construct and execute the plink clumping command
    output_path = os.path.join("tmp_GENAL", name)
    
    # Base command differs based on filetype
    base_cmd = f"{get_plink_path()} --memory {ram}"
    if filetype == "bed":
        base_cmd += f" --bfile {ref_path}"
    else:  # pgen
        base_cmd += f" --pfile {ref_path}"
        
    plink_command = f"{base_cmd} --rm-dup force-first --clump {to_clump_filename} --clump-kb {kb} \
                     --clump-r2 {r2} --clump-p1 {p1} --clump-p2 {p2} --out {output_path}"
    
    run_plink_command(plink_command)
    
    # Read log file to get the number of missing top variant IDs
    log_content = open(os.path.join("tmp_GENAL", f"{name}.log")).read()
    match = re.search(r"(\d+)\s+top\s+variant\s+ID", log_content)
    if match:
        missing_variants = int(match.group(1))
        print(f"Warning: {missing_variants} top variant IDs missing")

    if "No significant --clump results." in log_content:
        print("No SNPs remaining after clumping.")
        return
    
    match = re.search(r"(\d+)\s+clump[s]?\s+formed\s+from\s+(\d+)\s+index", log_content)
    if match:
        print(f"{match.group(1)} clumps formed from {match.group(2)} top variants.")

    # Extract the list of clumped SNPs and get the relevant data subset
    clumped_filename = os.path.join("tmp_GENAL", f"{name}.clumps")
    if not os.path.exists(clumped_filename):
        raise FileNotFoundError(f"'{clumped_filename}' is missing.")
    plink_clumped = pd.read_csv(clumped_filename, sep="\s+", usecols=["ID"])
    clumped_data = data[data["SNP"].isin(plink_clumped["ID"])]
    clumped_data.reset_index(drop=True, inplace=True)
    return clumped_data