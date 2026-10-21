from pyliftover import LiftOver
import os, subprocess
import numpy as np
import gzip
import shutil
import uuid
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

from .constants import CHR_CODES
from .tools import create_tmp, download_file, read_config

# Builds accepted for the start and end arguments (UCSC names) and their aliases
LIFT_BUILDS = ["hg17", "hg18", "hg19", "hg38"]
BUILD_ALIASES = {"grch37": "hg19", "grch38": "hg38"}
# Chain file chromosome names of the numeric codes 23 (X), 24 (Y) and 26 (MT)
CHAIN_CHR_NAMES = {23: "chrX", 24: "chrY", 26: "chrM"}
# Primary chromosomes of the target build: mappings to other contigs (alt, random, Un, fix) are treated as unlifted
PRIMARY_CONTIG_PATTERN = r"chr(\d+|X|Y|M)"
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def lift_data(
    data,
    start="hg19",
    end="hg38",
    extraction_file=False,
    chain_file=None,
    name=None,
    liftover_path=None,
    object_id="tmp_id",
):
    """
    Perform a liftover from one genetic build to another. If the chain file required for the liftover is not present, it will be downloaded. It"s also possible to manually provide the path to the chain file.
    If the dataset is large, it is suggested to use the liftOver executable (liftover_path argument).

    Args:
        data (pd.DataFrame): The input data containing at least "CHR" and "POS" columns (1-based positions).
        start (str, optional): The current build of the data: "hg17", "hg18", "hg19" or "hg38" ("GRCh37" and "GRCh38" are accepted for "hg19" and "hg38"). Defaults to "hg19".
        end (str, optional): The target build for liftover (same accepted values as start). Defaults to "hg38".
        extraction_file (bool, optional): If True, also prints a CHR POS SNP space-delimited file for extraction. Defaults to False.
        chain_file (str, optional): Path to a local chain file for the lift. Overrides the start and end arguments if provided.
        name (str, optional): Specify a filename or filepath (without extension) for saving. If not provided, the data is not saved.
        liftover_path (str, optional): Specify the path to the USCS liftover executable. If not provided, the lift will be done in python (slower for large amount of SNPs).
        object_id (str, optional): Specify the object id for tmp file writing (internal use only)

    Raises:
        ValueError: If required columns are missing, if the builds are not valid or if provided chain file path is incorrect.

    Returns:
        pd.DataFrame: Lifted data.

    Notes:
        Function for the :meth:`Geno.lift` method.
        Variants lifted to a non-primary contig (alt, random, Un, fix) are removed as unlifted.
        Chromosomes 23, 24 and 26 are lifted as chrX, chrY and chrM.
        The EA and NEA alleles of variants lifted to the opposite strand are reverse-complemented.
    """

    # Prepare chain file and get its path
    chain_path = prepare_chain_file(chain_file, start, end)

    # Prepare the data for lifting: handle missing values in CHR, POS columns
    nrows = data.shape[0]
    data.dropna(subset=["CHR", "POS"], inplace=True)
    data.reset_index(drop=True, inplace=True)  # unique labels for the label-based operations below
    # Remove absurd positions (positions are 1-based)
    data.drop(data[(data.POS < 1) | (data.POS >= 300_000_000)].index, inplace=True)
    data.reset_index(drop=True, inplace=True)
    data["CHR"] = data["CHR"].astype("int64")
    data["POS"] = data["POS"].astype("int64")
    n_na = nrows - data.shape[0]
    if n_na:
        print(
            f"Excluded {n_na} SNPs ({n_na/nrows*100:.3f}%) with NaN values in CHR or POS columns or with a position outside [1, 300,000,000)."
        )

    # Perform liftover with the liftover executable or in python
    if data.empty:
        print("No SNPs to lift.")
    elif liftover_path:
        data = lift_coordinates_liftover(data, object_id, chain_path, liftover_path)
    else:
        data = lift_coordinates_python(data, chain_path)

    # Handle post-liftover operations
    data = post_lift_operations(data, name, extraction_file, object_id)

    return data


def normalize_build(build, argument):
    """Return the UCSC name of a build ('GRCh37' -> 'hg19', 'GRCh38' -> 'hg38'). Raise a ValueError if it is not supported."""
    build_name = str(build).strip().lower()
    build_name = BUILD_ALIASES.get(build_name, build_name)
    if build_name not in LIFT_BUILDS:
        raise ValueError(
            f"Invalid {argument} build: {build!r}. Must be one of {', '.join(LIFT_BUILDS)} "
            f"('GRCh37' and 'GRCh38' are accepted for 'hg19' and 'hg38'), or provide a chain file with the chain_file argument."
        )
    return build_name


def prepare_chain_file(chain_file, start, end):
    """Handle chain file loading, downloading if necessary. Return its path."""
    if chain_file is not None:  # If a local chain file is provided
        if not os.path.isfile(chain_file):
            raise ValueError("The provided path does not lead to a valid file.")
        print(
            "You provided a path to a local chain path which will be used for the lift."
        )
        chain_path = chain_file
    else:  # Use the specified start and end builds to identify chain file
        start = normalize_build(start, "start")
        end = normalize_build(end, "end")
        if start == end:
            raise ValueError(f"The start and end builds are identical ({start}): there is nothing to lift.")
        # Construct chain filename
        chain_name = f"{start}To{end.capitalize()}.over.chain"
        config = read_config()
        ref_path = config["paths"]["ref_path"]
        chains_folder_path = os.path.join(ref_path, "chain_files")

        # Ensure directory for chain files exists
        if not os.path.exists(chains_folder_path):
            try:
                os.makedirs(chains_folder_path)
            except OSError:
                raise OSError(
                    f"Unable to create the '{chains_folder_path}' directory. Check permissions."
                )

        # Check for the chain file locally or download it if necessary
        chain_path = os.path.join(chains_folder_path, chain_name)
        if not os.path.isfile(chain_path):
            print(
                f"The chain file to lift from {start} to {end} was not found. Attempting to download it..."
            )
            # Download the chain file
            url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{start}/liftOver/{chain_name}.gz"
            gz_path = f"{chain_path}.gz"
            try:
                download_file(url, gz_path, description=f"{chain_name}.gz")
                # Decompress the downloaded file
                print(f"The download was successful. Unzipping...")
                with gzip.open(gz_path, "rb") as f_in, open(
                    chain_path, "wb"
                ) as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(gz_path)
            except Exception as e:
                if os.path.exists(gz_path):
                    os.remove(gz_path)
                print(f"The download was unsuccessful: {e}")
                print(
                    "Consider downloading the chain file manually from the UCSC website and providing its path via the chain_file argument."
                )
                raise FileNotFoundError("Chain file not found.")

    return chain_path


def _chain_chromosome(chrom):
    """Chromosome name used in chain files for a numeric CHR code (e.g. 1 -> 'chr1', 23 -> 'chrX')."""
    return CHAIN_CHR_NAMES.get(chrom, f"chr{chrom}")


def _reverse_complement(alleles):
    """Reverse-complement the alleles made of A, C, G, T letters (other values are returned unchanged)."""
    alleles = alleles.copy()
    is_dna = alleles.astype(str).str.fullmatch(r"[ACGTacgt]+") & alleles.notna()
    alleles[is_dna] = alleles[is_dna].str.translate(COMPLEMENT).str[::-1]
    return alleles


def _apply_lift(data, contigs, positions, strands):
    """
    Update CHR/POS with the lifted coordinates (1-based), remove unlifted variants and variants lifted to non-primary
    contigs, and reverse-complement the alleles of variants lifted to the opposite strand.

    Args:
        data (pd.DataFrame): Data being lifted.
        contigs, positions, strands (array-like): Target contig (None if unlifted), 1-based position and strand of each row.
    """
    nrows = data.shape[0]
    data.reset_index(drop=True, inplace=True)  # unique labels for the label-based operations below
    contigs = pd.Series(contigs, index=data.index, dtype=object)
    strands = pd.Series(strands, index=data.index, dtype=object)
    primary = contigs.str.fullmatch(PRIMARY_CONTIG_PATTERN, na=False).astype(bool)
    n_unmapped = int(contigs.isna().sum())
    n_other_contig = int((contigs.notna() & ~primary).sum())

    chr_names = contigs.where(primary).str[3:]
    chr_numbers = pd.to_numeric(chr_names, errors="coerce").fillna(chr_names.map(CHR_CODES).astype(float))
    data["CHR"] = chr_numbers.astype("Int64")
    data["POS"] = pd.to_numeric(pd.Series(positions, index=data.index), errors="coerce").astype("Int64")

    minus = primary & (strands == "-")
    n_minus = int(minus.sum())
    allele_columns = [col for col in ["EA", "NEA"] if col in data.columns]
    for col in allele_columns:
        data.loc[minus, col] = _reverse_complement(data.loc[minus, col])

    data.drop(index=data.index[~primary], inplace=True)
    data.reset_index(drop=True, inplace=True)

    n_lost = n_unmapped + n_other_contig
    if n_lost:
        print(f"{n_lost} SNPs ({n_lost/nrows*100:.3f}%) could not be lifted and are deleted from the final data.")
        if n_other_contig:
            print(
                f"Of these, {n_other_contig} SNPs were lifted to alternative, unplaced or patch contigs (e.g. chrUn_..., ..._alt, ..._random, ..._fix) and are counted as unlifted."
            )
    else:
        print("All SNPs have been lifted.")
    if data.empty:
        print("No SNPs could be lifted.")
    if n_minus:
        if allele_columns:
            print(
                f"{n_minus} SNPs were lifted to the opposite strand: their {' and '.join(allele_columns)} alleles have been reverse-complemented."
            )
        else:
            print(f"{n_minus} SNPs were lifted to the opposite strand (no EA/NEA columns to reverse-complement).")
    return data


def lift_coordinates_liftover(data, object_id, chain_path, liftover_path):
    """Lift data using the liftover executable and a chain file."""
    # Add the executable part if not there
    if not os.path.isfile(liftover_path):
        liftover_path = os.path.join(liftover_path, "liftOver")
    # Check that it is indeed the path to liftOver executable
    try:
        process = subprocess.run(
            [liftover_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=60,
            text=True,
        )
        if not process.stderr.startswith("liftOver"):
            raise TypeError(
                "The path provided is an executable, but not the liftOver executable. Check the path."
            )
    except Exception as e:
        raise TypeError(e)
    print("Lifting coordinates using liftOver.")

    # Write data in BED6 format for liftOver: 0-based half-open interval [POS-1, POS), row number as name, + strand
    create_tmp()
    to_lift_filename = os.path.join("tmp_GENAL", f"{object_id}.prelift")
    lifted_filename = os.path.join("tmp_GENAL", f"{object_id}.postlift")
    unmapped_filename = os.path.join("tmp_GENAL", f"{object_id}_unMapped")
    bed = pd.DataFrame(
        {
            "chrom": data["CHR"].map(_chain_chromosome),
            "start": data["POS"] - 1,
            "end": data["POS"],
            "name": np.arange(data.shape[0]),
            "score": 0,
            "strand": "+",
        }
    )
    bed.to_csv(to_lift_filename, sep="\t", index=False, header=False)

    # Call the liftOver software
    try:
        subprocess.run(
            [liftover_path, to_lift_filename, chain_path, lifted_filename, unmapped_filename],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running liftOver: {e}\n{e.stderr}")
        raise ValueError("Error running liftOver. Check error message for more details.")

    ## Read the output (the 1-based lifted position is the end of the interval) and match it to the input rows by name
    if os.path.getsize(lifted_filename) == 0:
        lifted = pd.DataFrame(columns=["chrom", "start", "end", "name", "score", "strand"])
    else:
        lifted = pd.read_csv(
            lifted_filename, sep="\t", header=None,
            names=["chrom", "start", "end", "name", "score", "strand"],
            dtype={"chrom": str, "strand": str},
        )
    lifted = lifted.drop_duplicates(subset=["name"], keep="first")
    rows = lifted["name"].astype(int).to_numpy()
    contigs = np.full(data.shape[0], None, dtype=object)
    positions = np.full(data.shape[0], np.nan)
    strands = np.full(data.shape[0], None, dtype=object)
    contigs[rows] = lifted["chrom"].to_numpy()
    positions[rows] = lifted["end"].to_numpy(dtype=float)
    strands[rows] = lifted["strand"].to_numpy()

    return _apply_lift(data, contigs, positions, strands)


def lift_coordinates_python(data, chain_path):
    """Perform liftover on data using the chain passed."""
    lo = LiftOver(chain_path)

    # Print message
    print("Lifting coordinates in python...")
    nrows = data.shape[0]
    if nrows > 500000:
        print("Your data is large, this can take a few minutes...")

    # Create a list of tuples and lift
    coordinates = list(zip(data["CHR"], data["POS"]))

    # Perform the lift: pyliftover uses 0-based positions, POS is 1-based
    def convert_coordinate(args):
        return lo.convert_coordinate(_chain_chromosome(args[0]), int(args[1]) - 1, "+")

    results = list(ThreadPoolExecutor().map(convert_coordinate, coordinates))

    # The first mapping (highest chain score) is used
    contigs = [res[0][0] if res else None for res in results]
    positions = [res[0][1] + 1 if res else np.nan for res in results]
    strands = [res[0][2] if res else None for res in results]
    return _apply_lift(data, contigs, positions, strands)


def post_lift_operations(data, name, extraction_file, object_id="tmp_id"):
    """Handle post-liftover operations like reporting, and saving results."""
    base_name = os.path.splitext(name)[0] if name else None
    if name:
        filename = base_name + ".txt"
        data.to_csv(f"{filename}", sep="\t", header=True, index=False)
        print(f"Lifted list of SNPs saved to {filename}")
    if extraction_file:
        if not ("SNP" in data.columns):
            data["SNP"] = data["CHR"].astype(str) + ":" + data["POS"].astype(str)
        extraction_filename = f"{base_name if base_name else object_id}_lifted_extraction.txt"
        data[["CHR", "POS", "SNP"]].to_csv(
            extraction_filename, sep=" ", header=False, index=False
        )
        print(f"Extraction file saved to {extraction_filename}")
    return data
