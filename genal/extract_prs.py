import pandas as pd
import os, subprocess, re, uuid
from functools import partial
from concurrent.futures import ProcessPoolExecutor

from .tools import check_bfiles, check_pfiles, setup_genetic_path, get_plink_path, create_tmp


MIN_RAM_PER_WORKER_MB = 3500 # Minimum RAM per PLINK process (conservative for large genotype files)

### ____________________
### PRS functions
### ____________________

def prs_func(data, weighted=True, path=None, ram=20000, cpus=4, name=None):
    """
    Compute a PRS (Polygenic Risk Score) using provided SNP-level data. Corresponds to the :meth:`Geno.prs` method
    """
    required = {"SNP", "EA", "BETA"}
    missing = required - set(data.columns)
    if missing: 
        raise ValueError(f"Missing required columns for PRS: {sorted(missing)}")
    if data.shape[0] == 0:
        raise ValueError("No SNPs were extracted from the genetic data and the PRS can't be computed.")

    # Get path and filetype
    path, filetype = setup_genetic_path(path)

    # Generate a default name if none is provided
    if name is None:
        name = str(uuid.uuid4())[:8]

    # Call extract_snps
    extracted_path = extract_snps_func(data.SNP, name, path, ram=ram, cpus=cpus)

    if extracted_path == "FAILED":
        raise ValueError("No SNPs were extracted from the genetic data and the PRS can't be computed.")

    # Additional check to ensure there are no duplicates in the data (need to think more about this, should be done upstream)
    data.drop_duplicates(subset=["SNP"], keep="first", inplace=True)
    if "CHR" in data.columns and "POS" in data.columns:
        data.drop_duplicates(subset=["CHR", "POS"], keep="first", inplace=True)
    
    # Write processed data to file and run plink on it
    data = data[["SNP", "EA", "BETA"]]
    data_path = os.path.join("tmp_GENAL", f"{name}_to_prs.txt")
    output_path = os.path.join("tmp_GENAL", f"{name}_prs")

    # Set BETA values to 1 if unweighted PRS is required
    if not weighted:
        data["BETA"] = 1
        print(f"Computing an unweighted PRS using {extracted_path} data.")
    else:
        print(f"Computing a weighted PRS using {extracted_path} data.")

    data.to_csv(data_path, sep="\t", index=False, header=True)
    
    # We can use --pfile since extract_snps now creates pgen files
    plink_command = f"{get_plink_path()} --memory {ram} --pfile {extracted_path} --threads {cpus} \
                     --score {data_path} 1 2 3 header --out {output_path} --allow-no-sex"

    # Check for empty dataframe
    n_snps = data.shape[0]
    if n_snps == 0:
        raise ValueError(
            "No SNPs remain for the polygenic risk score (PRS) calculation."
        )

    try:
        output = subprocess.run(
            plink_command, shell=True, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running PLINK command: {e}")
        print(f"PLINK stdout: {e.stdout}")
        print(f"PLINK stderr: {e.stderr}")
        raise ValueError("PLINK command failed. Check the error messages above for details.")

    # Read and process PRS results
    prs_file = output_path + ".sscore"
    log_file = output_path + ".log"
    if os.path.isfile(prs_file): #If the profile file exists: PRS was successful
        #Extracts the number of SNPs used for the PRS computation
        log_content = open(log_file).read()
        match = re.search(r'--score: (\d+) variant[s] processed', log_content)
        if match:
            n_predictors = int(match.group(1))
            print(
            f"The PRS computation was successful and used {n_predictors}/{n_snps} ({n_predictors/n_snps*100:.3f}%) SNPs."
            )
        else:
            print("Could not extract the number of SNPs used for the PRS computation.")
        #Return scores
        df_score = pd.read_csv(prs_file, sep="\s+")
        df_score.rename(columns={"#FID": "FID"}, inplace=True)
        return df_score
    else:
        print(output.stdout)
        raise ValueError(
            f"The PRS computation was not successful. Check the {output_path + '.log'} file."
        )


### _____________________
### Extract SNPs functions
### _____________________   

# We are currently excluding all multiallelic variants by forcing first on all duplicates. 
# Could be improved by keeping the relevant version of the multiallelic SNPs based on allele matching
def extract_snps_func(snp_list, name=None, path=None, ram=20000, cpus=4):
    """
    Extracts a list of SNPs from the given path. This function corresponds to the following Geno method: :meth:`Geno.extract_snps`.

    Args:
        snp_list (pd.Series): Series of SNPs to extract.
        name (str): Name prefix for the output files.
        path (str, optional): Path to the dataset. Defaults to the path from the configuration.

    Returns:
        str: path to the genetic files containing the extracted SNPs

    Raises:
        TypeError: Raises an error when no valid path is saved or when there's an incorrect format in the provided path.
    """
    # Check if snp_list is empty Series
    if snp_list.empty:
        print("The provided SNP list is empty.")
        return "FAILED"
    
    # Generate a default name if none is provided
    if name is None:
        name = str(uuid.uuid4())[:8]

    # Get path and filetype
    path, filetype = setup_genetic_path(path)
    create_tmp()  # Ensure tmp_GENAL exists

    # Prepare the SNP list
    snp_list = snp_list.dropna().drop_duplicates()
    snp_list_name = f"{name}_list.txt"
    snp_list_path = os.path.join("tmp_GENAL", snp_list_name)
    snp_list.to_csv(snp_list_path, sep=" ", index=False, header=None)
    nrow = len(snp_list)

    # Check if the data is split by chromosome
    filetype_split = "split" if "$" in path else "combined"

    output_path = os.path.join("tmp_GENAL", f"{name}_allchr")

    # Guard against empty SNP list (applies to both split and combined)
    if nrow == 0:
        print("The SNP list is empty after deduplication.")
        return "FAILED"

    if filetype_split == "split":
        # Calculate workers based on memory budget (not SNP count)
        max_workers_by_ram = max(1, int(ram // MIN_RAM_PER_WORKER_MB))
        workers = max(1, min(max_workers_by_ram, cpus, 22))  # Cap at 22 chromosomes, min 1

        # Allocate RAM per worker
        per_worker_ram = int(ram // workers)

        #print(f"Parallelizing extraction across {workers} workers with {per_worker_ram}MB RAM each")

        merge_command, bedlist_path = extract_snps_from_split_data(
            name, path, output_path, snp_list_path, filetype,
            workers=workers, per_worker_ram=per_worker_ram, ram=ram
        )
        handle_multiallelic_variants(name, merge_command, bedlist_path)
    else:
        extract_snps_from_combined_data(name, path, output_path, snp_list_path, filetype, ram=ram)

    #Check that at least 1 variant has been extracted. If not, return "FAILED" to warn downstream functions (prs, association_test)
    log_path = output_path + ".log"
    with open(log_path, 'r') as log_file:
        if " 0 variants remaining" in log_file.read():
            print("None of the provided SNPs were found in the genetic data.")
            return "FAILED"
        else:
            if check_pfiles(output_path):
                print(f"Created pgen/pvar/psam fileset with extracted SNPs: {output_path}")
            else:
                print(f"Could not extract the SNPs from the provided genetic data: check plink .log file")
            # Report SNPs not found
            report_snps_not_found(nrow, name)

    return output_path


def extract_command_parallel(task_id, name, path, snp_list_path, filetype, per_worker_ram=4000):
    """
    Helper function to run SNP extraction in parallel for different chromosomes.
    Args:
        task_id (int): Identifier for the task/chromosome.
        name (str): Name prefix for the output files.
        path (str): Path to the data set.
        snp_list_path (str): Path to the list of SNPs to extract.
        filetype (str): Type of genetic files ("bed" or "pgen")
        per_worker_ram (int): RAM limit in MB for this PLINK process.
    Returns:
        int: Returns the task_id if no valid files are found.
        dict: Returns error dict {'failed': True, 'chr': task_id, ...} if extraction fails.
        None: Returns None on success.
    """
    input_path = path.replace("$", str(task_id))

    # Check if files exist based on filetype
    if filetype == "bed" and not check_bfiles(input_path):
        return task_id
    elif filetype == "pgen" and not check_pfiles(input_path):
        return task_id

    output_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{task_id}")

    # Build command based on filetype
    base_cmd = f"{get_plink_path()}"
    if filetype == "bed":
        base_cmd += f" --bfile {input_path}"
    else:  # pgen
        base_cmd += f" --pfile {input_path}"

    command = f"{base_cmd} --extract {snp_list_path} --memory {per_worker_ram} --threads 1 --rm-dup force-first --make-pgen --out {output_path}"

    result = subprocess.run(
        command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )

    # Check for failures and return diagnostic info (diagnostics are in .log file)
    if result.returncode != 0:
        return {'failed': True, 'chr': task_id, 'log': f"{output_path}.log", 'returncode': result.returncode}

    # Also check if output files were created
    if not check_pfiles(output_path):
        return {'failed': True, 'chr': task_id, 'log': f"{output_path}.log", 'returncode': -1}

    return None  # Success


def extract_snps_from_split_data(name, path, output_path, snp_list_path, filetype, workers=4, per_worker_ram=4000, ram=20000):
    """Extract SNPs from data split by chromosome."""
    print("Extracting SNPs for each chromosome...")
    num_tasks = 22
    partial_extract_command_parallel = partial(
        extract_command_parallel,
        name=name,
        path=path,
        snp_list_path=snp_list_path,
        filetype=filetype,
        per_worker_ram=per_worker_ram
    )  # Wrapper function

    # First attempt with calculated workers
    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(
            executor.map(partial_extract_command_parallel, range(1, num_tasks + 1))
        )

    # Check for failures (non-None returns indicate errors)
    failed_chrs = [r for r in results if r is not None and isinstance(r, dict) and r.get('failed')]
    not_found = [r for r in results if r is not None and not isinstance(r, dict)]

    # Retry failed chromosomes with reduced workers if any failures occurred
    if failed_chrs and workers > 1:
        print(f"{len(failed_chrs)} chromosome(s) failed. Retrying with reduced parallelization...")
        retry_workers = max(1, workers // 2)
        # Recalculate RAM per worker based on original total budget
        total_ram_budget = per_worker_ram * workers
        per_worker_ram_retry = int(total_ram_budget // retry_workers)

        partial_retry = partial(
            extract_command_parallel,
            name=name,
            path=path,
            snp_list_path=snp_list_path,
            filetype=filetype,
            per_worker_ram=per_worker_ram_retry
        )

        failed_chr_ids = [r['chr'] for r in failed_chrs]
        with ProcessPoolExecutor(max_workers=retry_workers) as executor:
            retry_results = list(executor.map(partial_retry, failed_chr_ids))

        # Update results - surface errors for persistent failures
        for orig_id, retry_result in zip(failed_chr_ids, retry_results):
            if retry_result is not None and isinstance(retry_result, dict) and retry_result.get('failed'):
                # Still failed - surface the error
                log_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{orig_id}.log")
                if os.path.exists(log_path):
                    print(f"Chr{orig_id} failed after retry. Check log: {log_path}")
                    try:
                        with open(log_path, 'r') as f:
                            lines = f.readlines()
                            print(f"Last 10 lines of log:\n{''.join(lines[-10:])}")
                    except Exception:
                        pass

    # Merge extracted SNPs from each chromosome
    bedlist_name = f"{name}_bedlist.txt"
    bedlist_path = os.path.join("tmp_GENAL", bedlist_name)

    # Create the bedlist file
    output_name = os.path.join("tmp_GENAL", f"{name}_extract")
    with open(bedlist_path, "w+") as bedlist_file:
        found = []
        for i in range(1, 23):
            if i in not_found:
                print(f"bed/bim/fam or pgen/pvar/psam files not found for chr{i}.")
            elif check_pfiles(f"{output_name}_chr{i}"):
                bedlist_file.write(f"{output_name}_chr{i}\n")
                found.append(i)
                print(f"SNPs extracted for chr{i}.")
            else:
                print(f"No SNPs extracted for chr{i}.")

    if len(found) == 0:
        raise Warning("No SNPs were extracted from any chromosome.")
    
    # If only one chromosome was extracted, no need to merge, simply rename the files
    if len(found) == 1:
        chr_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{found[0]}")
        for ext in [".pgen", ".pvar", ".psam", ".log"]:
            os.rename(f"{chr_path}{ext}", f"{output_path}{ext}")
        return None, bedlist_path

    print("Merging SNPs extracted from each chromosome...")
    merge_command = f"{get_plink_path()} --memory {ram} --pmerge-list {bedlist_path} pfile --out {output_path}"
    try:
        subprocess.run(
            merge_command, shell=True, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running PLINK command: {e}")
        print(f"PLINK stdout: {e.stdout}")
        print(f"PLINK stderr: {e.stderr}")
        raise ValueError("PLINK command failed. Check the error messages above for details.")

    return merge_command, bedlist_path


def extract_snps_from_combined_data(name, path, output_path, snp_list_path, filetype, ram=20000):
    """Extract SNPs from combined data."""
    print("Extracting SNPs...")

    # Build command based on filetype
    base_cmd = f"{get_plink_path()}"
    if filetype == "bed":
        base_cmd += f" --bfile {path}"
    else:  # pgen
        base_cmd += f" --pfile {path}"

    extract_command = f"{base_cmd} --memory {ram} --extract {snp_list_path} --rm-dup force-first --make-pgen --out {output_path}"
    
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

    file_path = os.path.join("tmp_GENAL", f"{name}_allchr.pvar")
    extracted_snps_count = count_lines(file_path)-1 #pvar files include column names
    delta_nrow = nrow - extracted_snps_count
    if delta_nrow > 0:
        print(
            f"{delta_nrow}({delta_nrow/nrow*100:.3f}%) SNPs were not extracted from the genetic data."
        )

# TODO: Check if this function is still needed with plink2
def handle_multiallelic_variants(name, merge_command, bedlist_path):
    """Handle multiallelic variants detected during merging."""

    if merge_command is None:
        return
    
    def remove_multiallelic():
        missnp_path = os.path.join(
            "tmp_GENAL", f"{name}_allchr.vmiss"
            )
        if not os.path.exists(missnp_path):
            return 0
            
        snps_to_exclude = pd.read_csv(missnp_path, header=None)
        for i in range(1, 23):
            pvar_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{i}.pvar")
            if not os.path.isfile(pvar_path):
                continue
            pvar = pd.read_csv(pvar_path, sep="\t", header=None)
            # If no SNPs would be left for this chr: remove corresponding bedlist line
            n_to_exclude = len(set(pvar[2]).intersection(set(snps_to_exclude[0])))
            if n_to_exclude == len(set(pvar[2])):
                print(f"No SNPs remaining for chromosome {i}.")
                tmp_filename = os.path.join("tmp_GENAL", "tmp_multiallelic")
                with open(bedlist_path, "r") as file, open(
                    tmp_filename, "w"
                ) as temp_file:
                    output_name = os.path.join("tmp_GENAL", f"{name}_extract")
                    line_to_exclude = f"{output_name}_chr{i}\n"
                    for line in file:
                        if line != line_to_exclude:
                            temp_file.write(line)
                # Replace the original file with the temporary file
                os.replace(tmp_filename, bedlist_path)

            # If there is at least one multiallelic SNP for this chr
            elif n_to_exclude > 0:
                pfile_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{i}")
                command = f"{get_plink_path()} --pfile {pfile_path} --exclude {missnp_path} --make-pgen --out {pfile_path}"
                subprocess.run(
                    command,
                    shell=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
        return len(snps_to_exclude)

    log_content = open(os.path.join("tmp_GENAL", f"{name}_allchr.log")).read()
    if "Error: Multiple" in log_content:
        print("Multiallelic variants detected in the genetic files: removing them before merging.")
        n_multiallelic = remove_multiallelic()
        print(f"Reattempting the merge after exclusion of {n_multiallelic} multiallelic variants.")
        subprocess.run(
            merge_command,
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
