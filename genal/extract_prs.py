import numpy as np
import pandas as pd
import os, subprocess, re, uuid
from functools import partial
from concurrent.futures import ThreadPoolExecutor

from .tools import (
    check_pfiles,
    setup_genetic_path,
    get_plink_path,
    create_tmp,
    run_plink_command,
    split_chromosome_path,
    SPLIT_CHROMOSOMES,
)
from .geno_tools import missing_snp_ids


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
    data = data.copy()

    # Rows without SNP ID can't be extracted nor scored
    missing = missing_snp_ids(data["SNP"])
    if missing.any():
        print(f"{missing.sum()} rows with a missing SNP ID are removed from the PRS.")
        data = data[~missing.to_numpy()]
        if data.shape[0] == 0:
            raise ValueError("No SNPs were extracted from the genetic data and the PRS can't be computed.")

    # Get path and filetype
    path, filetype = setup_genetic_path(path)

    # Generate a default name if none is provided
    if name is None:
        name = str(uuid.uuid4())[:8]

    # Call extract_snps
    extracted_path = extract_snps_func(data.SNP, name, path, ram=ram, cpus=cpus, data=data)

    if extracted_path == "FAILED":
        raise ValueError("No SNPs were extracted from the genetic data and the PRS can't be computed.")

    # Each variant is scored once: remove rows with a duplicated SNP ID or with the same position and alleles
    n_rows = data.shape[0]
    data = data.drop_duplicates(subset=["SNP"], keep="first")
    if "CHR" in data.columns and "POS" in data.columns:
        variant_key = data[["CHR", "POS"]].copy()
        if "NEA" in data.columns:
            ea, nea = data["EA"].astype(str).to_numpy(), data["NEA"].astype(str).to_numpy()
            variant_key["A_min"], variant_key["A_max"] = np.where(ea < nea, ea, nea), np.where(ea < nea, nea, ea)
        else:
            variant_key["A_min"] = data["EA"].astype(str)
        duplicated = variant_key.duplicated(keep="first") & variant_key[["CHR", "POS"]].notna().all(axis=1)
        data = data[~duplicated.to_numpy()]
    n_duplicates = n_rows - data.shape[0]
    if n_duplicates > 0:
        print(f"{n_duplicates} rows with a duplicated SNP ID or variant (same position and alleles) have been removed.")

    # Write processed data to file and run plink on it
    data_path = os.path.join("tmp_GENAL", f"{name}_to_prs.txt")
    output_path = os.path.join("tmp_GENAL", f"{name}_prs")

    if weighted:
        weights = data[["SNP", "EA", "BETA"]].copy()
        print(f"Computing a weighted PRS using {extracted_path} data.")
    else:
        # Score the trait-increasing allele of each SNP with a weight of 1
        weights = data[["SNP", "EA"]].copy()
        negative = (data["BETA"] < 0).to_numpy()
        if "NEA" in data.columns:
            swap = negative & data["NEA"].notna().to_numpy()
            weights.loc[swap, "EA"] = data["NEA"].to_numpy()[swap]
        else:
            swap = np.zeros(data.shape[0], dtype=bool)
        weights["BETA"] = np.where(negative & ~swap, -1, 1)
        print(
            f"Computing an unweighted PRS using {extracted_path} data: each SNP contributes the dosage of its trait-increasing allele with a weight of 1."
        )
        if swap.sum() > 0:
            print(f"For {swap.sum()} SNPs with a negative BETA, the NEA is the trait-increasing allele and is the scored allele.")
        n_sign = (negative & ~swap).sum()
        if n_sign > 0:
            print(f"{n_sign} SNPs with a negative BETA have no NEA: their EA is scored with a weight of -1.")

    weights.to_csv(data_path, sep="\t", index=False, header=True)

    # Check for empty dataframe
    n_snps = weights.shape[0]
    if n_snps == 0:
        raise ValueError(
            "No SNPs remain for the polygenic risk score (PRS) calculation."
        )

    plink_command = [
        get_plink_path(), "--memory", str(ram), "--pfile", extracted_path, "--threads", str(cpus),
        "--score", data_path, "1", "2", "3", "header", "cols=+scoresums",
        "--out", output_path,
    ]
    output = run_plink_command(plink_command)

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
        print(
            "SCORE1_SUM is the sum over the SNPs used of the weight times the dosage of the scored allele. "
            "SCORE1_AVG is SCORE1_SUM divided by the number of alleles of these SNPs (2 per SNP)."
        )
        #Return scores
        df_score = pd.read_csv(prs_file, sep="\t")
        df_score.rename(columns={"#FID": "FID", "#IID": "IID"}, inplace=True)
        return df_score
    else:
        print(output.stdout)
        raise ValueError(
            f"The PRS computation was not successful. Check the {output_path + '.log'} file."
        )


### _____________________
### Extract SNPs functions
### _____________________

def extract_snps_func(snp_list, name=None, path=None, ram=20000, cpus=4, data=None):
    """
    Extracts a list of SNPs from the given path. This function corresponds to the following Geno method: :meth:`Geno.extract_snps`.

    Args:
        snp_list (pd.Series): Series of SNPs to extract.
        name (str): Name prefix for the output files.
        path (str, optional): Path to the dataset. Defaults to the path from the configuration.
        data (pd.DataFrame, optional): SNP-level data with the SNP and EA columns (NEA is used if present).
            When several records of the genetic data share a variant ID (e.g. the biallelic records of a split
            multiallelic variant), the record matching these alleles is kept (see :func:`resolve_duplicate_ids`).

    Returns:
        str: path to the genetic files containing the extracted SNPs ("FAILED" if no SNP was extracted)

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

    # Prepare the SNP list (missing or placeholder IDs can't be searched in the genetic data)
    missing = missing_snp_ids(snp_list)
    if missing.any():
        print(f"{missing.sum()} SNPs with a missing SNP ID can't be searched in the genetic data and are ignored.")
    snp_list = snp_list[~missing].drop_duplicates()
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

        extracted = extract_snps_from_split_data(
            name, path, output_path, snp_list_path, filetype,
            workers=workers, per_worker_ram=per_worker_ram, ram=ram, data=data
        )
        if not extracted:
            print("None of the provided SNPs were found in the genetic data.")
            return "FAILED"
    else:
        extract_snps_from_combined_data(name, path, output_path, snp_list_path, filetype, ram=ram)

    #Check that at least 1 variant has been extracted. If not, return "FAILED" to warn downstream functions (prs, association_test)
    log_path = output_path + ".log"
    with open(log_path, 'r') as log_file:
        if " 0 variants remaining" in log_file.read():
            print("None of the provided SNPs were found in the genetic data.")
            return "FAILED"
    if not check_pfiles(output_path):
        print(f"Could not extract the SNPs from the provided genetic data: check the {log_path} file.")
        return "FAILED"

    if filetype_split == "combined":
        counts = resolve_duplicate_ids(output_path, data, ram)
        report_duplicate_ids(counts)
        if not check_pfiles(output_path):
            print("No SNPs remain in the extracted genetic data.")
            return "FAILED"

    print(f"Created pgen/pvar/psam fileset with extracted SNPs: {output_path}")
    # Report SNPs not found
    report_snps_not_found(nrow, name)

    return output_path


def extract_command_parallel(task_id, name, path, snp_list_path, filetype, per_worker_ram=4000):
    """
    Helper function to run SNP extraction in parallel for different chromosomes.
    Args:
        task_id (int): Chromosome number (files of chromosomes 23-26 can also be named X, Y, XY, MT).
        name (str): Name prefix for the output files.
        path (str): Path to the data set.
        snp_list_path (str): Path to the list of SNPs to extract.
        filetype (str): Type of genetic files ("bed" or "pgen")
        per_worker_ram (int): RAM limit in MB for this PLINK process.
    Returns:
        int: Returns the task_id if no valid files are found.
        dict: Returns error dict {'failed': True, 'chr': task_id, ...} if extraction fails.
        None: Returns None on success or when none of the SNPs are on this chromosome.
    """
    input_path = split_chromosome_path(path, task_id, filetype)
    if input_path is None:
        return task_id

    output_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{task_id}")

    # Build command based on filetype
    command = [
        get_plink_path(), "--bfile" if filetype == "bed" else "--pfile", input_path,
        "--extract", snp_list_path, "--memory", str(per_worker_ram), "--threads", "1",
        "--make-pgen", "--out", output_path,
    ]

    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Check for failures and return diagnostic info (diagnostics are in .log file)
    log_path = f"{output_path}.log"
    if result.returncode != 0:
        # PLINK exits with an error when none of the SNPs are on this chromosome
        if os.path.exists(log_path) and " 0 variants remaining" in open(log_path).read():
            return None
        return {'failed': True, 'chr': task_id, 'log': log_path, 'returncode': result.returncode}

    # Also check if output files were created
    if not check_pfiles(output_path):
        return {'failed': True, 'chr': task_id, 'log': log_path, 'returncode': -1}

    return None  # Success


def extract_snps_from_split_data(name, path, output_path, snp_list_path, filetype, workers=4, per_worker_ram=4000, ram=20000, data=None):
    """Extract SNPs from data split by chromosome. Returns False if no SNPs were extracted from any chromosome."""
    print("Extracting SNPs for each chromosome...")
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
    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(
            executor.map(partial_extract_command_parallel, SPLIT_CHROMOSOMES)
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
        with ThreadPoolExecutor(max_workers=retry_workers) as executor:
            retry_results = list(executor.map(partial_retry, failed_chr_ids))
        failed_chrs = [r for r in retry_results if r is not None and isinstance(r, dict) and r.get('failed')]

    # Stop if the extraction still fails for some chromosomes (the results would be incomplete)
    if failed_chrs:
        for failure in failed_chrs:
            print(f"The extraction failed for chr{failure['chr']}. Check log: {failure['log']}")
            if os.path.exists(failure["log"]):
                with open(failure["log"], "r") as f:
                    print(f"Last 10 lines of log:\n{''.join(f.readlines()[-10:])}")
        raise RuntimeError(
            "The PLINK extraction failed for chromosome(s) "
            + ", ".join(f"{failure['chr']} (log: {failure['log']})" for failure in failed_chrs)
            + ". The analysis is stopped because its results would miss these chromosomes."
        )

    # Keep one record per variant ID in each chromosome (merging requires unique IDs)
    output_name = os.path.join("tmp_GENAL", f"{name}_extract")
    counts = {"matched": 0, "first": 0, "removed": 0}
    for i in SPLIT_CHROMOSOMES:
        if i not in not_found and check_pfiles(f"{output_name}_chr{i}"):
            counts_chr = resolve_duplicate_ids(f"{output_name}_chr{i}", data, per_worker_ram)
            counts = {key: counts[key] + counts_chr[key] for key in counts}

    # Merge extracted SNPs from each chromosome
    bedlist_name = f"{name}_bedlist.txt"
    bedlist_path = os.path.join("tmp_GENAL", bedlist_name)

    # Create the bedlist file
    with open(bedlist_path, "w+") as bedlist_file:
        found = []
        for i in SPLIT_CHROMOSOMES:
            if i in not_found:
                if i <= 22:
                    print(f"bed/bim/fam or pgen/pvar/psam files not found for chr{i}.")
            elif check_pfiles(f"{output_name}_chr{i}"):
                bedlist_file.write(f"{output_name}_chr{i}\n")
                found.append(i)
                print(f"SNPs extracted for chr{i}.")
            else:
                print(f"No SNPs extracted for chr{i}.")

    if len(found) == 0:
        report_duplicate_ids(counts)
        return False

    check_sample_sets({i: f"{output_name}_chr{i}" for i in found})

    # If only one chromosome was extracted, no need to merge, simply rename the files
    if len(found) == 1:
        report_duplicate_ids(counts)
        chr_path = os.path.join("tmp_GENAL", f"{name}_extract_chr{found[0]}")
        for ext in [".pgen", ".pvar", ".psam", ".log"]:
            os.rename(f"{chr_path}{ext}", f"{output_path}{ext}")
        return True

    print("Merging SNPs extracted from each chromosome...")
    merge_command = [get_plink_path(), "--memory", str(ram), "--pmerge-list", bedlist_path, "pfile", "--out", output_path]
    run_plink_command(merge_command)

    # IDs shared by records of different chromosomes
    counts_merged = resolve_duplicate_ids(output_path, data, ram)
    report_duplicate_ids({key: counts[key] + counts_merged[key] for key in counts})

    return True


def extract_snps_from_combined_data(name, path, output_path, snp_list_path, filetype, ram=20000):
    """Extract SNPs from combined data."""
    print("Extracting SNPs...")

    extract_command = [
        get_plink_path(), "--bfile" if filetype == "bed" else "--pfile", path,
        "--memory", str(ram), "--extract", snp_list_path, "--make-pgen", "--out", output_path,
    ]

    subprocess.run(extract_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def resolve_duplicate_ids(fileset, data, ram):
    """
    Keep a single record per variant ID in a pgen/pvar/psam fileset produced by an extraction.

    Records sharing a variant ID (e.g. the biallelic records of a split multiallelic variant) are resolved with the
    alleles of the SNP-level data: the first record whose alleles include both EA and NEA is kept.
    When the SNP-level data has no NEA for this SNP, or when no record matches,
    all the records with this ID are removed. Without SNP-level alleles (data is None or has no EA column),
    the first record is kept.

    Args:
        fileset (str): Path prefix of the pgen/pvar/psam fileset (modified in place).
        data (pd.DataFrame or None): SNP-level data with the SNP and EA columns (NEA optional).
        ram (int): RAM in MB for PLINK.

    Returns:
        dict: Number of duplicated IDs resolved to the record matching the alleles ("matched"),
            to their first record ("first"), and removed ("removed").
    """
    counts = {"matched": 0, "first": 0, "removed": 0}
    pvar_path = fileset + ".pvar"
    with open(pvar_path, "r") as f:
        lines = f.readlines()
    n_header = 0
    while n_header < len(lines) and lines[n_header].startswith("#"):
        n_header += 1
    columns = lines[n_header - 1].lstrip("#").rstrip("\n").split("\t")
    col = {c: columns.index(c) for c in ["ID", "REF", "ALT"]}
    fields = [line.rstrip("\n").split("\t") for line in lines[n_header:]]
    ids = pd.Series([f[col["ID"]] for f in fields])
    duplicated = ids.duplicated(keep=False) & (ids != ".")
    if not duplicated.any():
        return counts

    has_alleles = data is not None and "EA" in data.columns
    if has_alleles:
        info = data[~missing_snp_ids(data["SNP"]).to_numpy()].drop_duplicates(subset=["SNP"]).set_index("SNP")

    drop = set()
    for snp, rows in ids[duplicated].groupby(ids[duplicated]).groups.items():
        rows = list(rows)
        if not has_alleles or snp not in info.index:
            drop.update(rows[1:])
            counts["first"] += 1
            continue
        row = info.loc[snp]
        ea = str(row["EA"]).upper()
        nea = str(row["NEA"]).upper() if "NEA" in info.columns and pd.notna(row["NEA"]) else None
        keep = None
        if nea is not None and ea != nea:
            for r in rows:
                alleles = {fields[r][col["REF"]].upper(), *fields[r][col["ALT"]].upper().split(",")}
                if ea in alleles and nea in alleles:
                    keep = r
                    break
        if keep is None:
            drop.update(rows)
            counts["removed"] += 1
        else:
            drop.update(r for r in rows if r != keep)
            counts["matched"] += 1

    # Give the records to remove a missing ID, then keep only the records with an ID
    for r in drop:
        fields[r][col["ID"]] = "."
    with open(pvar_path, "w") as f:
        f.writelines(lines[:n_header])
        f.writelines("\t".join(field) + "\n" for field in fields)
    kept_ids = sorted(set(ids[~ids.index.isin(list(drop))]) - {"."})
    if len(kept_ids) == 0:
        for ext in [".pgen", ".pvar", ".psam"]:
            os.remove(fileset + ext)
        return counts
    kept_ids_path = fileset + "_unique_ids.txt"
    pd.Series(kept_ids).to_csv(kept_ids_path, index=False, header=False)
    unique_fileset = fileset + "_unique"
    run_plink_command([
        get_plink_path(), "--memory", str(ram), "--pfile", fileset,
        "--extract", kept_ids_path, "--make-pgen", "--out", unique_fileset,
    ])
    for ext in [".pgen", ".pvar", ".psam"]:
        os.replace(unique_fileset + ext, fileset + ext)
    return counts


def report_duplicate_ids(counts):
    """Print how variant IDs shared by several records of the genetic data were resolved."""
    n_duplicated = sum(counts.values())
    if n_duplicated == 0:
        return
    message = f"{n_duplicated} variant IDs are shared by several records of the genetic data (e.g. split multiallelic variants)."
    if counts["matched"] > 0:
        message += f" For {counts['matched']} of them, the record matching the EA/NEA alleles is used."
    if counts["removed"] > 0:
        message += f" {counts['removed']} of them are removed because no single record can be identified from the EA/NEA alleles."
    if counts["first"] > 0:
        message += f" For {counts['first']} of them, the first record is used (no allele information)."
    print(message)


def check_sample_sets(filesets):
    """Warn when the filesets extracted from each chromosome do not contain the same samples."""
    samples = {}
    for chrom, fileset in filesets.items():
        psam = pd.read_csv(fileset + ".psam", sep="\t", dtype=str, keep_default_na=False)
        id_cols = [c for c in psam.columns if c in ("#FID", "FID", "#IID", "IID", "SID")]
        samples[chrom] = set(zip(*[psam[c] for c in id_cols]))
    all_samples = set().union(*samples.values())
    incomplete = all_samples - set.intersection(*samples.values())
    if len(incomplete) > 0:
        details = ", ".join(
            f"chr{chrom}: {len(all_samples - s)}" for chrom, s in samples.items() if len(all_samples - s) > 0
        )
        print(
            f"Warning: the genetic files of the different chromosomes do not contain the same samples. "
            f"{len(incomplete)} of the {len(all_samples)} samples are absent from at least one chromosome file (absent samples per file: {details}). "
            f"The genotypes of a sample absent from a chromosome file are missing for all the SNPs of this chromosome: "
            f"they are mean-imputed in a PRS and the sample is excluded from the association tests of these SNPs."
        )


def report_snps_not_found(nrow, name):
    """Report the number of SNPs not found in the data."""

    file_path = os.path.join("tmp_GENAL", f"{name}_allchr.pvar")
    with open(file_path, "r") as file:
        extracted_snps_count = sum(1 for line in file if not line.startswith("#"))
    delta_nrow = nrow - extracted_snps_count
    if delta_nrow > 0:
        print(
            f"{delta_nrow}({delta_nrow/nrow*100:.3f}%) SNPs were not extracted from the genetic data."
        )
