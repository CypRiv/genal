import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile

import numpy as np
import pandas as pd
import requests

from .constants import REF_PANELS, REF_PANELS_URL, CONFIG_DIR, BUILDS, REF_PARQUET_URL, CHR_CODES, REF_DATA_VERSION

config_path = os.path.join(CONFIG_DIR, "config.json")
default_ref_path = os.path.join(CONFIG_DIR, "Reference_files")
PLINK_EXTENSIONS = (".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam")
SPLIT_CHROMOSOMES = range(1, 27)  # Chromosomes searched in paths split by chromosome ('$' placeholder)
CHR_NAMES = {23: "X", 24: "Y", 25: "XY", 26: "MT"}  # Names also searched for chromosomes 23-26
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}  # Complementary bases (strand flips)


def default_config():
    """Returns default config values"""
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
    default_config = {
        "paths": {
            "plink2_path": "",
            "liftover_path": "",
            "geno_path": "",
            "geno_filetype": "",
            "ref_path": default_ref_path,
            "ref_filetype": "bed"
        }
    }
    return default_config


def read_config():
    """Get config file data"""
    with open(config_path, "r") as f:
        config = json.load(f)
    return config


def reference_data_dir():
    """Folder of the reference data downloaded from the genal bucket: the version subfolder (REF_DATA_VERSION) of the reference folder."""
    return os.path.join(read_config()["paths"]["ref_path"], REF_DATA_VERSION)


def write_config(config):
    """Write data to config file (atomically: concurrent readers see the previous or the new file, never a partial one)"""
    os.makedirs(CONFIG_DIR, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=CONFIG_DIR, prefix=".config_", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(config, f, indent=4)
        os.replace(tmp_path, config_path)
    except BaseException:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    return


def download_file(url, destination, description="file", chunk_size=8192, timeout=120):
    """Download a file with cleanup on partial failures."""
    destination_dir = os.path.dirname(destination)
    if destination_dir:
        os.makedirs(destination_dir, exist_ok=True)

    print(f"Downloading {description}...")
    try:
        with requests.get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            with open(destination, "wb") as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
    except (requests.RequestException, OSError) as e:
        if os.path.exists(destination):
            os.remove(destination)
        raise RuntimeError(f"Failed to download {description}: {e}")

    print("Download complete.")
    return destination

def strip_plink_extension(path):
    """Remove a PLINK file extension (.bed, .bim, .fam, .pgen, .pvar, .psam) from a path, keeping any other dots."""
    root, ext = os.path.splitext(path)
    return root if ext in PLINK_EXTENSIONS else path

def setup_genetic_path(path):
    """Configure the genetic data path based on user input and saved configuration."""
    config = read_config()
    if path is None:
        if not config["paths"].get("geno_path"):
            raise TypeError("No path has been saved in the config file. Please provide one.")
        path = config["paths"]["geno_path"]
        filetype = config["paths"].get("geno_filetype", "bed")  # Default to bed if not specified
        print(f"Using path saved in config file: {path}")
        return path, filetype

    # Ensure correct path format
    if path.count("$") > 1:
        raise TypeError(
            "The path should contain at most 1 '$' sign. Use it to indicate the chromosome number if the data is split by chromosomes."
        )
    # Check if the path is valid
    path = strip_plink_extension(path)  # Remove file extension
    if path.count("$") == 0:
        is_bed = check_bfiles(path)
        is_pgen = check_pfiles(path)
        if not (is_bed or is_pgen):
            raise TypeError("The path does not lead to valid bed/bim/fam or pgen/pvar/psam files.")
        filetype = "bed" if is_bed else "pgen"
    else:
        check_split_bed = [split_chromosome_path(path, i, "bed") is not None for i in SPLIT_CHROMOSOMES]
        check_split_pgen = [split_chromosome_path(path, i, "pgen") is not None for i in SPLIT_CHROMOSOMES]
        if not (any(check_split_bed) or any(check_split_pgen)):
            raise TypeError("The path does not lead to valid bed/bim/fam or pgen/pvar/psam files.")
        filetype = "bed" if any(check_split_bed) else "pgen"

    # Save path and filetype to config.json file (used when path=None)
    if (config["paths"].get("geno_path"), config["paths"].get("geno_filetype")) != (path, filetype):
        config["paths"]["geno_path"] = path
        config["paths"]["geno_filetype"] = filetype
        write_config(config)
        print(f"Saved '{path}' as the default genotype path (used when path=None).")
    return path, filetype

def split_chromosome_path(path, chrom, filetype):
    """
    Get the path of the genetic files of one chromosome, for genetic data split by chromosome.

    The '$' sign of the path is replaced by the chromosome number and, for chromosomes 23, 24, 25 and 26,
    also by their name (X, Y, XY and MT): for instance, "data_chr$" leads to "data_chr23" or "data_chrX".

    Args:
        path (str): Path with a '$' sign in place of the chromosome.
        chrom (int): Chromosome number.
        filetype (str): "bed" for bed/bim/fam files or "pgen" for pgen/pvar/psam files.

    Returns:
        str or None: Path of the genetic files of this chromosome (without extension), None if there are none.

    Raises:
        ValueError: If genetic files exist under both the number and the name of the chromosome.
    """
    check_files = check_bfiles if filetype == "bed" else check_pfiles
    labels = [str(chrom)] + ([CHR_NAMES[chrom]] if chrom in CHR_NAMES else [])
    found = [path.replace("$", label) for label in labels if check_files(path.replace("$", label))]
    if len(found) > 1:
        raise ValueError(
            f"Genetic files of chromosome {chrom} exist under several names ({', '.join(found)}). Keep only one of them."
        )
    return found[0] if found else None

def check_bfiles(filepath):
    """Check if the path specified leads to a bed/bim/fam triple."""
    if (
        os.path.exists("{}.bed".format(filepath))
        and os.path.exists("{}.bim".format(filepath))
        and os.path.exists("{}.fam".format(filepath))
    ):
        return True
    return False

def check_pfiles(filepath):
    """Check if the path specified leads to a pgen/pvar/psam triple."""
    if (
        os.path.exists("{}.pgen".format(filepath))
        and os.path.exists("{}.pvar".format(filepath))
        and os.path.exists("{}.psam".format(filepath))
    ):
        return True
    return False

def create_tmp():
    """Create the temporary folder if not present"""
    if not os.path.exists("tmp_GENAL"):
        try:
            os.makedirs("tmp_GENAL")
        except OSError:
            raise OSError(
                "Unable to create the 'tmp_GENAL' directory. Check permissions."
            )
        
def delete_tmp():
    """Delete the tmp_GENAL folder of the current directory."""
    if os.path.isdir("tmp_GENAL"):
        def _onerror(func, path, exc_info):
            if isinstance(exc_info[1], FileNotFoundError):
                return
            raise exc_info[1]
        shutil.rmtree("tmp_GENAL", onerror=_onerror)
        print("The tmp_GENAL folder has been successfully deleted.")
    else:
        print("There is no tmp_GENAL folder to delete in the current directory.")
    return

def set_reference_folder(path=""):
    """
    Set a folder path to store reference data.

    This function allows users to specify a directory where reference data will be stored.
    If the directory doesn't exist, it will be created. If no path is provided, a default
    directory named 'Reference_files' in the .genal folder at root will be used.

    Parameters:
        path (str, optional): The desired directory path for storing reference data. Defaults to a temporary folder in the current working directory.

    Raises:
        OSError: If the directory cannot be created.

    Returns:
        None: The function prints messages to inform the user of the status and any errors.
    """

    # If no path is provided, set default path to root/.genal/Reference_files
    if not path:
        path = default_ref_path
        print(f"No path provided, defaulting to {default_ref_path}.")

    # If the directory doesn't exist, attempt to create it
    if not os.path.isdir(path):
        try:
            os.makedirs(path)
            print(f"Creating the '{path}' directory.")
        except OSError:
            raise OSError(
                f"Unable to create the '{path}' directory. Check permissions."
            )

    # Check if the directory is readable
    if not os.access(path, os.R_OK):
        print(f"Error: The directory '{path}' is not readable.")
        return

    # Check if the directory is writable
    if not os.access(path, os.W_OK):
        print(f"Error: The directory '{path}' is not writable.")
        return

    # Update the configuration with the new path
    config = read_config()
    config["paths"]["ref_path"] = path
    write_config(config)

    print(f"Reference files will be downloaded and stored in: '{path}'")


def get_reference_panel_path(reference_panel="EUR_38"):
    """
    Retrieve the path of the specified reference panel.

    This function checks if the provided reference panel is a valid path to bed/bim/fam or pgen/pvar/psam files.
    If not, it checks if the reference panel exists in the reference folder (`<ref_path>/v2/<PANEL>/`, ref_path being
    the reference folder of the configuration file, ~/.genal/Reference_files by default). If it doesn't exist,
    the function downloads it.

    Args:
        reference_panel (str): Reference panel identifier (e.g., "EUR_38", "AFR_37")
            or path to custom reference panel (bed/bim/fam or pgen/pvar/psam files). Default is "EUR_38".

    Returns:
        tuple: (path to reference panel, filetype ['bed' or 'pgen'])

    Raises:
        ValueError: If reference_panel is neither a valid path nor a valid panel identifier, or is a path split by
            chromosome ('$'), which must first be merged into one fileset.
        RuntimeError: If the download or extraction of the reference panel fails.
    """

    # Remove file extension and check if it's a path to a bed/bim/fam or pgen/pvar/psam triple
    reference_panel = strip_plink_extension(reference_panel)
    if check_bfiles(reference_panel):
        print(f"Using the provided path as the reference panel (bed format).")
        return reference_panel, "bed"
    elif check_pfiles(reference_panel):
        print(f"Using the provided path as the reference panel (pgen format).")
        return reference_panel, "pgen"
    else:
        if "$" in reference_panel:
            raise ValueError(
                f"Reference panels split by chromosome ('$' in {reference_panel!r}) are not supported. Merge the chromosome "
                "filesets into one with PLINK 2 (plink2 --pmerge-list <file listing the fileset prefixes> --make-pgen "
                "--out <merged prefix>; add 'bfile' after the list file for bed/bim/fam filesets) and use the merged prefix."
            )

        # Standardize panel name
        reference_panel = reference_panel.upper()
        
        # Validate reference panel name
        if reference_panel not in REF_PANELS:
            raise ValueError(
                f"Invalid reference panel: {reference_panel}. Must be one of: {', '.join(REF_PANELS)} "
                "or a path to a valid reference panel in bed/bim/fam or pgen/pvar/psam format."
            )
        
        # Set paths
        panel_dir = os.path.join(reference_data_dir(), reference_panel)
        panel_path = os.path.join(panel_dir, reference_panel)
        population, build = reference_panel.split("_")
        
        # Check if panel already exists
        if os.path.exists(panel_dir):
            if check_pfiles(panel_path):
                print(f"Using the {population} (build {build}) path as the reference panel (pgen format).")
                return panel_path, "pgen"
            elif check_bfiles(panel_path):
                print(f"Using the {population} (build {build}) path as the reference panel (bed format).")
                return panel_path, "bed"
        
        # Download and extract panel
        print(f"Downloading reference panel {reference_panel} in {panel_dir} ...")
        os.makedirs(panel_dir, exist_ok=True)
        
        # Download tar.gz file
        url = REF_PANELS_URL.format(panel=reference_panel)
        tar_path = os.path.join(panel_dir, f"{reference_panel}.tar.gz")
        try:
            download_file(url, tar_path, description=f"reference panel {reference_panel}")
            print(f"Extracting {reference_panel}...")
            with tarfile.open(tar_path) as tar:
                tar.extractall(panel_dir)
            os.remove(tar_path)
        except Exception as e:
            if os.path.exists(tar_path):
                os.remove(tar_path)
            raise RuntimeError(f"Failed to download/extract reference panel: {e}")
        
        # Default reference panels are in pgen format
        if check_pfiles(panel_path):
            return panel_path, "pgen"
        else:
            raise RuntimeError(f"Reference panel files not found after extraction")


def load_reference_panel(reference_panel="38"):
    """
    Load the reference panel variants.
    
    Args:
        reference_panel (str): Can be:
            - A path to a .bim or .pvar file
            - A build number ("37" or "38")
            - A population with build number (e.g. "EUR_37", "AFR_38"), where only the build number is considered
            If a build number is provided, loads the provided reference variants for that build (based on 1000G phase 3)
    
    Returns:
        pd.DataFrame: Reference panel DataFrame with standardized columns
    """
    # Make sure it's a string
    reference_panel = str(reference_panel)

    # Check if it's a path to a .bim or .pvar file
    reference_panel = strip_plink_extension(reference_panel)
    if os.path.exists(reference_panel + ".bim"):
        print(f"Using the provided bim file as the reference dataset.")
        reference_panel_df = check_reference_panel(read_variant_file(reference_panel + ".bim"))

    elif os.path.exists(reference_panel + ".pvar"):
        print(f"Using the provided pvar file as the reference dataset.")
        reference_panel_df = check_reference_panel(read_variant_file(reference_panel + ".pvar"))

    # Else, check if it points to a standard 37 or 38 reference panel
    else:
        # Extract build number if population is provided
        if reference_panel.upper() in REF_PANELS:
            _, build = reference_panel.upper().split("_")

        else:
            build = reference_panel
            if build not in BUILDS:
                raise ValueError(
                    f"Invalid reference panel: {reference_panel}. Must be one of: {', '.join(REF_PANELS)} "
                    "\n or a path to a valid reference panel in bed/bim/fam or pgen/pvar/psam format."
                )
        
        # Set up paths
        variants_dir = os.path.join(reference_data_dir(), "reference_variants")
        variants_file = os.path.join(variants_dir, f"reference_variants_{build}.parquet")
        
        # Download if not exists
        if not os.path.exists(variants_file):
            print(f"Downloading reference variants for build {build} in {variants_dir} ...")
            os.makedirs(variants_dir, exist_ok=True)
            
            # Download parquet file
            url = REF_PARQUET_URL.format(build=build)
            try:
                download_file(
                    url,
                    variants_file,
                    description=f"reference variants for build {build}",
                )
            except Exception as e:
                if os.path.exists(variants_file):
                    os.remove(variants_file)
                raise RuntimeError(f"Failed to download reference variants: {e}")
        
        # Load parquet file
        try:
            import pyarrow.parquet as pq

            reference_panel_df = arrow_to_pandas(pq.read_table(variants_file))
            print(f"Using reference variants from build {build}.")
        except Exception as e:
            raise RuntimeError(f"Failed to load reference variants: {e}")
        
    
    return reference_panel_df

def normalize_chr(chrom):
    """
    Convert chromosome codes to numbers: remove a 'chr' prefix and map X, Y, XY and MT (or M) to 23, 24, 25 and 26.

    Args:
        chrom (pd.Series): Chromosome codes (numbers or strings).

    Returns:
        pd.Series: Chromosome numbers (float), NaN for codes that are not chromosomes (e.g. unplaced contigs).
    """
    codes = pd.Series(pd.unique(chrom))
    names = codes.astype(str).str.strip().str.replace(r"^chr", "", regex=True, case=False).str.upper()
    numbers = pd.to_numeric(names.replace({code: str(number) for code, number in CHR_CODES.items()}), errors="coerce")
    return chrom.map(pd.Series(numbers.to_numpy(), index=codes)).astype(float)

def _read_variants_arrow(path, n_meta, names, tab_only):
    """
    Variant table of :func:`read_variant_file` read with pyarrow's multi-threaded CSV reader (tab delimiter), with the
    columns and values of the pandas reader: CHR as category, strings as read (no missing values), POS as read.
    Returns None when the file needs the pandas reader: a pyarrow parsing error, or, for files that pandas splits on
    any whitespace (tab_only), a field that is empty or contains whitespace or a missing position (pandas would merge
    consecutive separators or split the field).

    Args:
        path (str): Path to the .bim or .pvar file.
        n_meta (int): Number of '##' header lines to skip.
        names (list or None): Column names of a file without header line (.bim column order), None for a .pvar file
            with a '#CHROM' header line.
        tab_only (bool): Whether pandas reads the file with any whitespace as separator.
    """
    import pyarrow as pa
    import pyarrow.csv as pacsv
    import pyarrow.compute as pc

    if names is None:
        columns, string_columns, category = ["#CHROM", "POS", "ID", "REF", "ALT"], ["ID", "REF", "ALT"], "#CHROM"
    else:
        columns, string_columns, category = names, [c for c in names if c not in ("CHR", "POS")], "CHR"
    try:
        table = pacsv.read_csv(
            path,
            read_options=pacsv.ReadOptions(skip_rows=n_meta, column_names=names),
            parse_options=pacsv.ParseOptions(delimiter="\t"),
            convert_options=pacsv.ConvertOptions(
                include_columns=columns,
                # POS as integers (a value that is not an integer, e.g. a date, leaves the file to the pandas reader)
                column_types={category: pa.dictionary(pa.int32(), pa.string()), "POS": pa.int64(), **{c: pa.string() for c in string_columns}},
                strings_can_be_null=False,
            ),
        )
    except (pa.ArrowException, KeyError):  # Any pyarrow failure: the pandas reader reads the file (or raises its error)
        return None
    if tab_only:
        for column in [category, *string_columns]:
            values = table[column]
            if pa.types.is_dictionary(values.type):
                values = values.cast(pa.string())
            if pc.any(pc.or_(pc.equal(values, ""), pc.match_substring_regex(values, r"\s"))).as_py():
                return None
        del values
        if table["POS"].null_count > 0:
            return None
    table = table.select([c for c in columns if c != "CM"])
    return arrow_to_pandas(table)


def arrow_to_pandas(table):
    """
    Convert a pyarrow table to a DataFrame with the dtypes of pandas' readers (strings: object with pandas 2, str with
    pandas 3). The conversion is column by column, each Arrow column being freed once converted, and Python strings are
    not deduplicated: the deduplication is very slow on large tables and saves little, single-character alleles being
    shared by Python anyway.
    """
    import pyarrow as pa

    string_dtype = pd.Series([""], dtype=str).dtype
    types_mapper = None if string_dtype == object else {pa.string(): string_dtype, pa.large_string(): string_dtype}.get
    return table.to_pandas(types_mapper=types_mapper, deduplicate_objects=False, self_destruct=True, split_blocks=True)


def read_variant_file(path):
    """
    Read the variants of a .bim or .pvar file.

    .pvar files: the '##' header lines are skipped and the columns are selected by name from the '#CHROM' header line
    (QUAL, FILTER, INFO or CM columns are ignored), with A1 = REF and A2 = ALT. A .pvar file without '#CHROM' header line
    is read in .bim column order, as PLINK does.
    .bim files: whitespace-delimited columns CHR, SNP, CM (optional), POS, A1, A2.
    Chromosome codes are converted to numbers (see :func:`normalize_chr`), variants on other contigs or with a
    position that is not an integer are removed and missing variant IDs ('.') are set to NaN.

    Args:
        path (str): Path to the .bim or .pvar file.

    Returns:
        pd.DataFrame: Variants with the columns CHR, POS, SNP, A1, A2.
    """
    n_meta, first_line = 0, None
    with open(path, "r") as f:
        for line in f:
            if line.startswith("##"):
                n_meta += 1
                continue
            first_line = line
            break
    if first_line is None:
        return pd.DataFrame(columns=["CHR", "POS", "SNP", "A1", "A2"])

    if path.endswith(".pvar") and first_line.startswith("#CHROM"):
        variants = _read_variants_arrow(path, n_meta, None, tab_only=False)
        if variants is None:
            variants = pd.read_csv(
                path, sep="\t", skiprows=n_meta, usecols=["#CHROM", "POS", "ID", "REF", "ALT"],
                dtype={"#CHROM": "category", "ID": str, "REF": str, "ALT": str}, keep_default_na=False,
            )
        variants = variants.rename(columns={"#CHROM": "CHR", "ID": "SNP", "REF": "A1", "ALT": "A2"})
    else:
        n_columns = len(first_line.split())
        if n_columns == 6:
            names = ["CHR", "SNP", "CM", "POS", "A1", "A2"]
        elif n_columns == 5:
            names = ["CHR", "SNP", "POS", "A1", "A2"]
        else:
            raise ValueError(f"The variant file {path} should have 5 or 6 columns (.bim format) or a '#CHROM' header line (.pvar format).")
        variants = _read_variants_arrow(path, n_meta, names, tab_only=True)
        if variants is None:
            variants = pd.read_csv(
                path, sep=r"\s+", header=None, names=names, skiprows=n_meta,
                dtype={"CHR": "category", "SNP": str, "A1": str, "A2": str}, keep_default_na=False,
            )

    variants["CHR"] = normalize_chr(variants["CHR"])
    variants["POS"] = pd.to_numeric(variants["POS"], errors="coerce")
    invalid = variants["CHR"].isna() | variants["POS"].isna() | (variants["POS"] != np.floor(variants["POS"]))
    if invalid.any():
        print(f"{invalid.sum()} variants of {path} without a valid chromosome (1-22, X, Y, XY, MT) or position are ignored.")
        variants = variants[~invalid].reset_index(drop=True)
    variants["CHR"] = variants["CHR"].astype(int)
    variants["POS"] = variants["POS"].astype(int)
    missing_id = variants["SNP"] == "."
    if missing_id.any():
        variants["SNP"] = variants["SNP"].mask(missing_id)
    columns = ["CHR", "POS", "SNP", "A1", "A2"]
    return variants if list(variants.columns) == columns else variants[columns]

def check_reference_panel(df):
    """
    Check and standardize the reference panel DataFrame.

    Args:
        df (pd.DataFrame): Reference panel DataFrame with required columns
            ["CHR", "POS", "SNP", "A1", "A2"]

    Returns:
        pd.DataFrame: Processed reference panel with standardized columns

    Raises:
        ValueError: If required columns are missing or if data types are invalid
    """
    df = df.copy()

    # Convert allele columns to uppercase strings
    for col in ["A1", "A2"]:
        df[col] = df[col].astype(str).str.upper()

    # Convert chromosome codes to numbers ('chr' prefix removed, X/Y/XY/MT coded 23/24/25/26)
    chrom = normalize_chr(df["CHR"])
    if chrom.isna().any():
        raise ValueError("Chromosome (CHR) column of the reference panel must contain chromosome numbers or X, Y, XY, MT codes")
    df["CHR"] = chrom.astype(int)

    # Convert POS to integer
    try:
        df["POS"] = df["POS"].astype(int)
    except ValueError:
        raise ValueError("Position (POS) column of the reference panel must contain integer values")
        
    return df

def is_plink2_installed(plink_path):
        try:
            result = subprocess.run([plink_path, '--version'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True,
                                    check=True)
            if 'plink v2.0' in result.stdout.lower():
                return True
            else:
                return False
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False


def set_plink(path=""):
    """
    Set the path to the plink 2.0 executable, verify that it is the correct version and save it in the configuration file.

    Args:
        path (str): Path to the plink2 executable, or to the folder containing it.

    Raises:
        TypeError: If no path is provided, or if the path does not lead to an executable plink 2.0.
    """
    if not path:
        raise TypeError("You need to provide a path.")

    if not os.path.isfile(path):
        path = os.path.join(path, "plink2")

    try:
        process = subprocess.run(
            [path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60, text=True
        )

        if not os.access(path, os.X_OK):
            raise TypeError("You do not have permission to execute the plink executable.")

        if not "PLINK v2" in process.stdout:
            raise TypeError(
                "The path provided is an executable, but not the plink 2.0 executable. Check the path and plink version.\
                 Otherwise, run genal.install_plink() to install it automatically."
            )
    except Exception as e:
        raise TypeError(e)

    config = read_config()
    config["paths"]["plink2_path"] = path
    write_config(config)

    print(f"Path to plink 2.0 successfully set: '{path}'")
    return


def get_plink_path():
    """
    Return the path to the plink 2.0 executable saved in the configuration file.

    Returns:
        str: Path to the plink2 executable.

    Raises:
        ValueError: If no plink 2.0 path has been set (see :func:`set_plink` and :func:`install_plink`).
    """
    config = read_config()
    if not config["paths"]["plink2_path"]:
        raise ValueError(
            "The path to plink 2.0 has not been set yet. Use genal.set_plink(path_to_plink) first or run genal.install_plink() to install it automatically."
        )
    else:
        return config["paths"]["plink2_path"]


def install_plink(path=None):
    """
    Download plink 2.0 for the current operating system and save its path in the configuration file.

    Args:
        path (str, optional): Folder where plink 2.0 is installed. Defaults to ~/.genal/plink2/.
            If this folder already contains a working plink 2.0 executable, it is used without downloading.

    Raises:
        ValueError: If the operating system is not supported.
        RuntimeError: If the download, the extraction or the installation fails.
    """
    # Use default path if none provided
    if not path:
        path = os.path.join(CONFIG_DIR, "plink2")
        print(f"You have not specified a path for the installation of plink. The following directory will be used: {path}")

    # Set up paths based on OS
    system = platform.system()
    system_arch = platform.architecture()[0][:2]
    
    if system == 'Windows':
        plink_path = os.path.join(path, 'plink2.exe')
    else:
        plink_path = os.path.join(path, 'plink2')
        
    if is_plink2_installed(plink_path):
        print(f"Plink2 is already installed at {plink_path}")
        config = read_config()
        config["paths"]["plink2_path"] = plink_path
        write_config(config)
        print(f"Path set to: {plink_path}")
        return
    
    if not os.path.isdir(path):
        try:
            os.makedirs(path, exist_ok=True)
        except OSError:
            raise OSError(f"Unable to create '{path}' directory")

    # Determine correct PLINK2 build based on OS, architecture and CPU features
    if system == "Linux":
        if system_arch == "32":
            download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_i686_20241222.zip"
            build_type = "32-bit Linux"
        else:  # 64-bit
            # Check for AVX2 support and CPU vendor
            cpu_info = ""
            try:
                with open('/proc/cpuinfo', 'r') as f:
                    cpu_info = f.read().lower()
            except:
                print("Warning: Could not read CPU info, defaulting to standard Linux 64-bit build")
                
            has_avx2 = 'avx2' in cpu_info
            is_amd = 'amd' in cpu_info
            
            if has_avx2:
                if is_amd:
                    download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_amd_avx2_20241222.zip"
                    build_type = "64-bit Linux AMD with AVX2"
                else:  # Intel
                    download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20241222.zip"
                    build_type = "64-bit Linux Intel with AVX2"
            else:
                download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241222.zip"
                build_type = "64-bit Linux standard"
            
    elif system == "Darwin":  # macOS
        # Check for Apple Silicon
        is_arm = platform.machine().lower() == 'arm64'
        
        # Check for AVX2 support on Intel macs
        has_avx2 = False
        if not is_arm:
            try:
                sysctl = subprocess.check_output(['sysctl', '-a']).decode()
                has_avx2 = 'avx2' in sysctl.lower()
            except:
                print("Warning: Could not determine AVX2 support, defaulting to standard macOS 64-bit build")
        
        if is_arm:
            download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_arm64_20241222.zip"
            build_type = "macOS ARM64 (Apple Silicon)"
        elif has_avx2:
            download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_avx2_20241222.zip"
            build_type = "macOS Intel with AVX2"
        else:
            download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_20241222.zip"
            build_type = "macOS Intel standard"
            
    elif system == "Windows":
        if system_arch == "32":
            download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win32_20241222.zip"
            build_type = "32-bit Windows"
        else:  # 64-bit
            # Check for AVX2 support
            has_avx2 = False
            try:
                output = subprocess.check_output(['wmic', 'cpu', 'get', 'caption']).decode()
                has_avx2 = 'avx2' in output.lower()
            except:
                print("Warning: Could not determine AVX2 support, defaulting to standard Windows 64-bit build")
                
            if has_avx2:
                download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win_avx2_20241222.zip"
                build_type = "64-bit Windows with AVX2"
            else:
                download_url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win64_20241222.zip"
                build_type = "64-bit Windows standard"
    else:
        raise ValueError("Unsupported operating system for automatic PLINK 2 installation")

    print(f"Detected system configuration: {build_type}")
    
    # Download and extract PLINK2
    create_tmp()
    zip_path = os.path.join("tmp_GENAL", 'plink2.zip')
    
    print(f"Downloading PLINK 2 from {download_url}...")
    try:
        response = requests.get(download_url, stream=True)
        response.raise_for_status()
        with open(zip_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to download PLINK 2: {e}")
        
    print("Extracting PLINK 2...")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(path)
    except zipfile.BadZipFile as e:
        raise RuntimeError(f"Failed to extract PLINK 2: {e}")
        
    os.remove(zip_path)
    
    # Set executable permissions and update config
    try:
        os.chmod(plink_path, 0o755)
    except OSError as e:
        print(f"Warning: Could not set executable permissions: {e}")

    if is_plink2_installed(plink_path):
        print("PLINK 2 installed successfully")
        config = read_config()
        config["paths"]["plink2_path"] = plink_path
        write_config(config)
        print(f"Path set to: {plink_path}")
    else:
        raise RuntimeError("PLINK 2 installation failed")

    return

def run_plink_command(command):
    """
    Execute a PLINK command with proper error handling.

    Args:
        command (list or str): The PLINK command to execute. A list of arguments is run directly (no shell,
            so paths may contain spaces); a string is run through the shell.

    Returns:
        subprocess.CompletedProcess: The completed PLINK process (with captured stdout and stderr).

    Raises:
        RuntimeError: If the PLINK command fails because of insufficient memory.
        ValueError: If the PLINK command fails for another reason (the PLINK output is printed, and the message
            contains the PLINK error).
    """
    try:
        return subprocess.run(command, shell=isinstance(command, str), capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        if "Out of memory" in f"{e.stderr}{e.stdout}":
            raise RuntimeError(
                f"PLINK command failed due to insufficient memory.\n"
                f"Try increasing the RAM allocation. Example:\n"
                f"  # If geno_object is your Geno object\n"
                f"  geno_object.ram = 25000  # Allocates 25GB of RAM for plink commands\n"
            )
        else:
            print(f"Error running PLINK command: {e}")
            print(f"PLINK stdout: {e.stdout}")
            print(f"PLINK stderr: {e.stderr}")
            raise ValueError(plink_error_message(e.stderr))


def plink_error_message(stderr):
    """
    Message of a failed PLINK command: the PLINK error (from its 'Error' line to the end of stderr, line breaks removed),
    with a hint when the covariates are collinear (VIF_TOO_HIGH or CORR_TOO_HIGH).
    """
    lines = [line.strip() for line in (stderr or "").splitlines() if line.strip()]
    start = next((i for i, line in enumerate(lines) if line.startswith("Error")), None)
    error = " ".join(lines[start:] if start is not None else lines[-3:])
    message = f"PLINK command failed: {error}" if error else "PLINK command failed (see the PLINK output printed above)."
    if "VIF_TOO_HIGH" in error or "CORR_TOO_HIGH" in error:
        message += (
            " The covariates are collinear (rescaling does not help): center continuous covariates before squaring them "
            "or forming interactions, and remove redundant covariates."
        )
    return message
