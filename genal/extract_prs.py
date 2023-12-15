import pandas as pd
import numpy as np
import os, subprocess, re
from functools import partial
from concurrent.futures import ProcessPoolExecutor

from .tools import check_bfiles, read_config, write_config, setup_genetic_path, get_plink19_path


def prs_func(data, weighted=True, path=None, ram=10000, name=""):
    """
    Compute a PRS (Polygenic Risk Score) using provided SNP-level data. Corresponds to the :meth:`Geno.prs` method

    Args:
        data (pd.DataFrame): Dataframe containing SNP-level data.
        weighted (bool, optional): Perform a weighted PRS using BETA column estimates.
            If False, perform an unweighted PRS (equivalent to BETAs set to 1). Defaults to True.
        path (str, optional): Path to bed/bim/fam set of genetic files to use for PRS calculation.
            If not provided, it uses the genetic data extracted with the `.extract_snps` method or attempts to use the path in the config file. Defaults to None.
        ram (int, optional): RAM memory in MB to be used by plink. Defaults to 10000.
        name (str, optional): Name used for naming output and intermediate files. Defaults to "".

    Returns:
        pd.DataFrame: DataFrame containing PRS results.

    Raises:
        ValueError: If mandatory columns are missing in the data.
        TypeError: If valid bed/bim/fam files are not found.
        ValueError: If PRS computation was not successful.
    """
    
    # Call extract_snps
    path = extract_snps_func(data.SNP, name, path)

    # Set BETA values to 1 if unweighted PRS is required
    if not weighted:
        data["BETA"] = 1
        print(f"Computing an unweighted PRS using {path} data.")
    else:
        print(f"Computing a weighted PRS using {path} data.")

    # Check for empty dataframe
    n_snps = data.shape[0]
    if n_snps == 0:
        raise ValueError(
            "No SNPs remain for the polygenic risk score (PRS) calculation."
        )

    # Write processed data to file and run plink on it
    data = data[["SNP", "EA", "BETA"]]
    data_path = os.path.join("tmp_GENAL", f"{name}_to_prs.txt")
    output_path = os.path.join("tmp_GENAL", f"{name}_prs")

    data.to_csv(data_path, sep="\t", index=False, header=True)
    plink_command = f"{get_plink19_path()} --memory {ram} --bfile {path} \
                     --score {data_path} 1 2 3 header --out {output_path}"
    output = subprocess.run(
        plink_command, shell=True, capture_output=True, text=True, check=True
    )

    # Read and process PRS results
    prs_file = output_path + ".profile"
    log_file = output_path + ".log"
    if os.path.isfile(prs_file): #If the profile file exists: PRS was successful
        #Extracts the number of SNPs used for the PRS computation
        pattern = r'--score: (\d+) valid predictors loaded.'
        with open(log_file, 'r') as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    n_predictors = int(match.group(1))
        #Return scores
        print(
            f"The PRS computation was successful and used {n_predictors}/{n_snps} ({n_predictors/n_snps*100:.3f}%) SNPs."
        )
        df_score = pd.read_csv(prs_file, sep="\s+")
        return df_score[["FID", "IID", "SCORE"]]
    else:
        print(output.stdout)
        raise ValueError(
            f"The PRS computation was not successful. Check the {output_path + '.log'} file."
        )


def extract_snps_func(snp_list, name, path=None):
    """
    Extracts a list of SNPs from the given path. This function corresponds to the following Geno method: :meth:`Geno.extract_snps`.

    Args:
        snp_list (List[str]): List of SNPs to extract.
        name (str): Name prefix for the output files.
        path (str, optional): Path to the dataset. Defaults to the path from the configuration.

    Raises:
        TypeError: Raises an error when no valid path is saved or when there's an incorrect format in the provided path.
    """
    path = setup_genetic_path(path)

    snp_list, snp_list_path, nrow = prepare_snp_list(snp_list, name)

    filetype = "split" if "$" in path else "combined"

    output_path = os.path.join("tmp_GENAL", f"{name}_allchr")
    if filetype == "split":
        merge_command, bedlist_path = extract_snps_from_split_data(
            name, path, output_path, snp_list_path
        )
        handle_multiallelic_variants(name, merge_command, bedlist_path)
    else:
        extract_snps_from_combined_data(name, path, output_path, snp_list_path)
        
    print(f"Created bed/bim/fam fileset with extracted SNPs: {output_path}")

    # Report SNPs not found
    report_snps_not_found(nrow, name)

    return output_path


def prepare_snp_list(snp_list, name):
    """Prepare the SNP list for extraction."""
    snp_list = snp_list.dropna()
    snp_list_name = f"{name}_list.txt"
    snp_list_path = os.path.join("tmp_GENAL", snp_list_name)
    snp_list.to_csv(snp_list_path, sep=" ", index=False, header=None)
    nrow = len(snp_list)
    return snp_list, snp_list_path, nrow


def extract_command_parallel(task_id, name, path, snp_list_path):
    """
    Helper function to run SNP extraction in parallel for different chromosomes.
    Args:
        task_id (int): Identifier for the task/chromosome.
        name (str): Name prefix for the output files.
        path (str): Path to the data set.
        snp_list_path (str): Path to the list of SNPs to extract.
    Returns:
        int: Returns the task_id if no valid bed/bim/fam files are found.
    """
    bfile_path = path.replace("$", str(task_id))

    if not check_bfiles(bfile_path):
        return task_id

    output_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{task_id}")
    command = f"{get_plink19_path()} --bfile {bfile_path} --extract {snp_list_path} --make-bed --out {output_path}"
    subprocess.run(
        command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )


def create_bedlist(bedlist, output_name, not_found):
    """
    Creates a bedlist file for SNP extraction.
    Args:
        bedlist (str): Path to save the bedlist file.
        output_name (str): Base name for the output files.
        not_found (List[int]): List of chromosome numbers for which no bed/bim/fam files were found.
    """
    with open(bedlist, "w+") as bedlist_file:
        for i in range(1, 23):
            if i in not_found:
                print(f"bed/bim/fam files not found for chr{i}.")
            elif check_bfiles(f"{output_name}_chr{i}"):
                bedlist_file.write(f"{output_name}_chr{i}\n")
                print(f"SNPs extracted for chr{i}.")
            else:
                print(f"No SNPs extracted for chr{i}.")


def extract_snps_from_split_data(name, path, output_path, snp_list_path):
    """Extract SNPs from data split by chromosome."""
    print("Extracting SNPs for each chromosome...")
    num_tasks = 22
    partial_extract_command_parallel = partial(
        extract_command_parallel, name=name, path=path, snp_list_path=snp_list_path
    )  # Wrapper function
    with ProcessPoolExecutor() as executor:
        not_found = list(
            executor.map(partial_extract_command_parallel, range(1, num_tasks + 1))
        )

    # Merge extracted SNPs from each chromosome
    bedlist_name = f"{name}_bedlist.txt"
    bedlist_path = os.path.join("tmp_GENAL", bedlist_name)
    create_bedlist(
        bedlist_path, os.path.join("tmp_GENAL", f"{name}_extract"), not_found
    )
    print("Merging SNPs extracted from each chromosome...")
    merge_command = f"{get_plink19_path()} --merge-list {bedlist_path} --make-bed --out {output_path}"
    subprocess.run(
        merge_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )

    return merge_command, bedlist_path


def handle_multiallelic_variants(name, merge_command, bedlist_path):
    """Handle multiallelic variants detected during merging."""

    def remove_multiallelic():
        snps_to_exclude = pd.read_csv(
            os.path.join("tmp_GENAL", f"{name}_allchr-merge.missnp"), header=None
        )
        for i in range(1, 23):
            bim_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{i}.bim")
            if not os.path.isfile(bim_path):
                continue
            bim = pd.read_csv(bim_path, sep="\t", header=None)
            # If no SNPs would be left for this chr: remove corresponding bedlist line
            n_to_exclude = len(set(bim[1]).intersection(set(snps_to_exclude[0])))
            if n_to_exclude == len(set(bim[1])):
                print(f"No SNPs remaining for chromosome {i}.")
                tmp_filename = os.path.join("tmp_GENAL", "tmp_multiallelic")
                with open(bedlist_path, "r") as file, open(
                    tmp_filename, "w"
                ) as temp_file:
                    output_name = os.path.join("tmp_GENAL", f"{name}_extract")
                    line_to_exclude = f"{output_name}_chr{i}\n"
                    for current_line_number, line in enumerate(file, start=1):
                        if line != line_to_exclude:
                            temp_file.write(line)
                # Replace the original file with the temporary file
                os.replace(tmp_filename, bedlist_path)

            # If there is at least one multiallelic SNP for this chr
            elif n_to_exclude > 0:
                bfile_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{i}")
                missnp_path = os.path.join(
                    "tmp_GENAL", f"{name}_allchr-merge.missnp"
                )
                output_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{i}")
                command = f"{get_plink19_path()} --bfile {bfile_path} --exclude {missnp_path} --make-bed --out {output_path}"
                subprocess.run(
                    command,
                    shell=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
        return len(snps_to_exclude)

    log_content = open(os.path.join("tmp_GENAL", f"{name}_allchr.log")).read()
    if "with 3+ alleles present" in log_content:
        print("Multiallelic variants detected in the genetic files: removing them before merging.")
        n_multiallelic = remove_multiallelic()
        print(f"Reattempting the merge after exclusion of {n_multiallelic} multiallelic variants.")
        subprocess.run(
            merge_command,
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        ## In rare cases, plink fails to identify all the multiallelic SNPs with the first pass.
        log_content = open(os.path.join("tmp_GENAL", f"{name}_allchr.log")).read()
        if "with 3+ alleles present" in log_content:
            print("Multiallelic variants still detected: removing them before merging.")
            remove_multiallelic()
            print("Reattempting the merge after deletion of multiallelic variants.")
            subprocess.run(
                merge_command,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )


def extract_snps_from_combined_data(name, path, output_path, snp_list_path):
    """Extract SNPs from combined data."""
    print("Extracting SNPs...")
    extract_command = f"{get_plink19_path()} --bfile {path} --extract {snp_list_path} --make-bed --out {output_path}"
    subprocess.run(
        extract_command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    


def report_snps_not_found(nrow, name):
    """Report the number of SNPs not found in the data."""

    def count_lines(filepath):
        with open(filepath, "r") as file:
            return sum(1 for line in file)

    file_path = os.path.join("tmp_GENAL", f"{name}_allchr.bim")
    extracted_snps_count = count_lines(file_path)
    delta_nrow = nrow - extracted_snps_count
    if delta_nrow > 0:
        print(
            f"Extraction completed. {delta_nrow}({delta_nrow/nrow*100:.3f}%) SNPs were not extracted from the genetic data."
        )
