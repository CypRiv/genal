import os
import subprocess
import pandas as pd

from .geno_tools import *
from .tools import *


def clump_data(
    data,
    reference_panel="eur",
    plink19_path=get_plink19_path(),
    kb=250,
    r2=0.1,
    p1=5e-8,
    p2=0.01,
    name="noname",
    ram=10000,
    checks=[],
):
    """
    Perform clumping on the given data using plink. Corresponds to the :meth:`GENO.clump` method.
    
    Args:
        data (pd.DataFrame): Input data with at least 'SNP' and 'P' columns.
        reference_panel (str): The reference population for linkage disequilibrium values. Accepts values "eur", "sas", "afr", "eas", "amr". Alternatively, a path leading to a specific bed/bim/fam reference panel can be provided. Default is "eur".
        plink19_path (str): Path to the plink 1.9 software.
        kb (int, optional): Clumping window in terms of thousands of SNPs. Default is 250.
        r2 (float, optional): Linkage disequilibrium threshold, values between 0 and 1. Default is 0.1.
        p1 (float, optional): P-value threshold during clumping. SNPs above this value are not considered. Default is 5e-8.
        p2 (float, optional): P-value threshold post-clumping to further filter the clumped SNPs. If p2 < p1, it won't be considered. Default is 0.01.
        name (str): Name used for the files created in the tmp_GENAL folder.
        ram (int): Amount of RAM in MB to be used by plink.
        checks (list): List of column checks already performed on the data.
    
    Returns:
        pd.DataFrame: Data after clumping, if any.
        list: Updated checks.
    """

    # Ensure required columns exist in the data
    for column in ["SNP", "P"]:
        if column not in data.columns:
            raise ValueError(f"The column {column} is not found in the data")

    # Validate and process SNP and P columns, if not already done
    if "SNP" not in checks:
        data = check_snp_column(data)
        checks.append("SNP")
    if "P" not in checks:
        initial_rows = data.shape[0]
        data.dropna(subset=["SNP", "P"], inplace=True)
        deleted_rows = initial_rows - data.shape[0]
        if deleted_rows > 0:
            print(
                f"{deleted_rows} ({deleted_rows/initial_rows*100:.3f}%) rows with NA values in columns SNP or P have been deleted."
            )
        checks.append("P")

    # Create tmp directory if it doesn't exist
    if not os.path.exists("tmp_GENAL"):
        try:
            os.makedirs("tmp_GENAL")
        except OSError:
            raise OSError(
                "Unable to create the 'tmp_GENAL' directory. Check permissions."
            )

    # Save the relevant data columns to a temporary file
    to_clump_filename = os.path.join("tmp_GENAL", f"{name}_to_clump.txt")
    data[["SNP", "P"]].to_csv(to_clump_filename, index=False, sep="\t")

    # Construct and execute the plink clumping command
    output_path = os.path.join("tmp_GENAL", name)
    plink_command = f"{plink19_path} --memory {ram} --bfile {get_reference_panel_path(reference_panel)} \
                     --clump {to_clump_filename} --clump-kb {kb} --clump-r2 {r2} --clump-p1 {p1} \
                     --clump-p2 {p2} --out {output_path}"
    output = subprocess.run(
        plink_command, shell=True, capture_output=True, text=True, check=True
    )

    # Check and print the outputs for relevant information
    if output.returncode != 0:
        raise RuntimeError(
            f"PLINK execution failed with the following error: {output.stderr}"
        )
    if "more top variant IDs missing" in output.stderr:
        missing_variants = output.stderr.split("more top variant IDs missing")[0].split(
            "\n"
        )[-1]
        print(f"Warning: {missing_variants} top variant IDs missing")
    if "No significant --clump results." in output.stderr:
        print("No SNPs remaining after clumping.")
        return
    print(output.stdout.split("--clump: ")[1].split("\n")[0])

    # Extract the list of clumped SNPs and get the relevant data subset
    clumped_filename = os.path.join("tmp_GENAL", f"{name}.clumped")
    if not os.path.exists(clumped_filename):
        raise FileNotFoundError(f"'{clumped_filename}' is missing.")
    plink_clumped = pd.read_csv(clumped_filename, sep="\s+", usecols=["SNP"])
    clumped_data = data[data["SNP"].isin(plink_clumped["SNP"])]
    clumped_data.reset_index(drop=True, inplace=True)
    return clumped_data, checks
