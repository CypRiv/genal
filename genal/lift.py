from pyliftover import LiftOver
import os, subprocess
import numpy as np
import gzip
import shutil
import uuid
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv

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

    Args:
        data (pd.DataFrame): The input data containing at least "CHR" and "POS" columns (1-based positions).
        start (str, optional): The current build of the data: "hg17", "hg18", "hg19" or "hg38" ("GRCh37" and "GRCh38" are accepted for "hg19" and "hg38"). Defaults to "hg19".
        end (str, optional): The target build for liftover (same accepted values as start). Defaults to "hg38".
        extraction_file (bool, optional): If True, also prints a CHR POS SNP space-delimited file for extraction. Defaults to False.
        chain_file (str, optional): Path to a local chain file for the lift. Overrides the start and end arguments if provided.
        name (str, optional): Specify a filename or filepath (without extension) for saving. If not provided, the data is not saved.
        liftover_path (str, optional): Path to the UCSC liftOver executable, to lift with it instead of in python. If not provided (default), the lift is done in python, which is faster.
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


def _chain_chromosomes(chrom):
    """Chain file chromosome names (see :func:`_chain_chromosome`) of an array of CHR codes, computed once per distinct code."""
    codes, values = pd.factorize(chrom)
    return np.array([_chain_chromosome(value) for value in values], dtype=object)[codes]


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
    # The contig names are processed once per distinct name
    contig_codes, names = pd.factorize(pd.Series(contigs, dtype=object))
    names = pd.Series(np.asarray(names, dtype=object), dtype=object)
    primary_names = names.str.fullmatch(PRIMARY_CONTIG_PATTERN, na=False).astype(bool)
    chr_names = names.where(primary_names).str[3:]
    chr_numbers = pd.to_numeric(chr_names, errors="coerce").fillna(chr_names.map(CHR_CODES).astype(float)).to_numpy()
    lifted = contig_codes >= 0
    primary = np.zeros(nrows, dtype=bool)
    primary[lifted] = primary_names.to_numpy()[contig_codes[lifted]]
    n_unmapped = int((~lifted).sum())
    n_other_contig = int((lifted & ~primary).sum())

    chromosomes = np.full(nrows, np.nan)
    chromosomes[primary] = chr_numbers[contig_codes[primary]]
    data["CHR"] = pd.Series(chromosomes, index=data.index).astype("Int64")
    data["POS"] = pd.to_numeric(pd.Series(positions, index=data.index), errors="coerce").astype("Int64")

    minus = primary & (pd.Series(strands, index=data.index, dtype=object) == "-").to_numpy()
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
    n = data.shape[0]
    bed = pa.table(
        {
            "chrom": _chain_chromosomes(data["CHR"]),
            "start": data["POS"].to_numpy() - 1,
            "end": data["POS"].to_numpy(),
            "name": np.arange(n),
            "score": np.zeros(n, dtype=np.int64),
            "strand": np.full(n, "+", dtype=object),
        }
    )
    pacsv.write_csv(bed, to_lift_filename, write_options=pacsv.WriteOptions(include_header=False, delimiter="\t", quoting_style="none"))

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


def read_chain(chain_path):
    """
    Read the aligned blocks of a UCSC chain file (gzip-compressed if its name ends with .gz).

    Returns:
        tuple: (blocks, chains). blocks: DataFrame with one row per aligned block of non-zero size: source contig,
            source interval [start, end) (0-based), start of the target interval and chain (row of chains).
            chains: DataFrame with the score, target contig, target size and target strand of each chain.

    Raises:
        ValueError: If the file is not a valid chain file.
    """
    opener = gzip.open if str(chain_path).lower().endswith(".gz") else open
    chains, blocks = [], []
    source_sizes, target_sizes = {}, {}  # Size of each contig, which must be the same in all chains
    in_chain = False  # Whether the last chain read has not reached its last block (a line with one number) yet
    with opener(chain_path, "rt") as f:
        for line in f:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] == "chain":
                if in_chain:
                    raise ValueError(f"A chain of {chain_path} does not end with its last block (a line with one number).")
                if len(fields) < 12 or fields[4] != "+" or fields[9] not in ("+", "-"):
                    raise ValueError(f"Invalid chain header in {chain_path}: {line.strip()}")
                for sizes, contig, size in ((source_sizes, fields[2], int(fields[3])), (target_sizes, fields[7], int(fields[8]))):
                    if sizes.setdefault(contig, size) != size:
                        raise ValueError(f"The chains of {chain_path} give different sizes to {contig} ({sizes[contig]} and {size}).")
                chains.append((int(fields[1]), fields[7], int(fields[8]), fields[9]))
                source, source_pos, source_end = fields[2], int(fields[5]), int(fields[6])
                target_pos, target_end = int(fields[10]), int(fields[11])
                in_chain = True
                continue
            if not in_chain or len(fields) not in (1, 3):
                raise ValueError(f"Invalid alignment line in {chain_path}: {line.strip()}")
            size = int(fields[0])
            if size > 0:
                blocks.append((source, source_pos, source_pos + size, target_pos, len(chains) - 1))
            if len(fields) == 3:
                source_pos += size + int(fields[1])
                target_pos += size + int(fields[2])
            elif source_pos + size != source_end or target_pos + size != target_end:
                raise ValueError(f"The alignment blocks of a chain of {chain_path} do not match its interval.")
            else:
                in_chain = False
    if in_chain:
        raise ValueError(f"A chain of {chain_path} does not end with its last block (a line with one number).")
    # Integer columns also when the chain file has no aligned blocks (nothing is lifted)
    blocks = pd.DataFrame(blocks, columns=["source", "start", "end", "target_start", "chain"]).astype(
        {"start": np.int64, "end": np.int64, "target_start": np.int64, "chain": np.int64}
    )
    chains = pd.DataFrame(chains, columns=["score", "target", "target_size", "strand"]).astype(
        {"score": np.int64, "target_size": np.int64}
    )
    return blocks, chains


def lift_positions(chrom, pos, chain_path):
    """
    Lift CHR/POS coordinates with a chain file, as the first mapping of pyliftover's LiftOver.convert_coordinate: the
    position is mapped through the aligned block that contains it, in the chain with the highest score. Positions in
    regions where blocks of several chains overlap (there are none in the UCSC hg19/hg38 chain files) are lifted by
    pyliftover itself.

    Args:
        chrom (np.ndarray): Chromosome codes (1-26, see :func:`_chain_chromosome`).
        pos (np.ndarray): 1-based positions.
        chain_path (str): Path to the chain file.

    Returns:
        tuple: Target contig (None if unlifted), 1-based target position (NaN if unlifted) and target strand of each position.
    """
    blocks, chains = read_chain(chain_path)
    n = len(pos)
    pos0 = np.asarray(pos, dtype=np.int64) - 1
    block_of_row = np.full(n, -1)
    overlapping = np.zeros(n, dtype=bool)
    contigs = _chain_chromosomes(chrom)
    contig_codes, contig_names = pd.factorize(contigs)
    blocks_by_contig = dict(tuple(blocks.groupby("source", sort=False)))
    for code, contig in enumerate(contig_names):
        contig_blocks = blocks_by_contig.get(contig)
        if contig_blocks is None:
            continue
        contig_blocks = contig_blocks.sort_values("start", kind="stable")
        start, end = contig_blocks["start"].to_numpy(), contig_blocks["end"].to_numpy()
        rows = np.flatnonzero(contig_codes == code)
        p = pos0[rows]
        i = np.searchsorted(start, p, side="right") - 1
        inside = (i >= 0) & (p < end[np.maximum(i, 0)])
        # Groups of overlapping blocks (of different chains): their positions are lifted by pyliftover
        new_group = np.ones(len(start), dtype=bool)
        new_group[1:] = start[1:] >= np.maximum.accumulate(end)[:-1]
        group_size = np.bincount(np.cumsum(new_group) - 1)
        if (group_size > 1).any():
            group_end = np.maximum.reduceat(end, np.flatnonzero(new_group))
            k = np.searchsorted(start[new_group], p, side="right") - 1
            k0 = np.maximum(k, 0)
            in_overlap = (k >= 0) & (p < group_end[k0]) & (group_size[k0] > 1)
            overlapping[rows[in_overlap]] = True
            inside &= ~in_overlap
        block_of_row[rows[inside]] = contig_blocks.index.to_numpy()[i[inside]]

    lifted = block_of_row >= 0
    block = block_of_row[lifted]
    chain = blocks["chain"].to_numpy()[block]
    target = blocks["target_start"].to_numpy()[block] + (pos0[lifted] - blocks["start"].to_numpy()[block])
    strand = chains["strand"].to_numpy()[chain]
    target = np.where(strand == "-", chains["target_size"].to_numpy()[chain] - 1 - target, target)
    target_contigs = np.full(n, None, dtype=object)
    target_positions = np.full(n, np.nan)
    target_strands = np.full(n, None, dtype=object)
    target_contigs[lifted] = chains["target"].to_numpy()[chain]
    target_positions[lifted] = target + 1
    target_strands[lifted] = strand
    if overlapping.any():
        lo = LiftOver(chain_path)
        for row in np.flatnonzero(overlapping):
            result = lo.convert_coordinate(contigs[row], int(pos0[row]), "+")
            if result:  # The first mapping (highest chain score) is used
                target_contigs[row], target_positions[row], target_strands[row] = result[0][0], result[0][1] + 1, result[0][2]
    return target_contigs, target_positions, target_strands


def lift_coordinates_python(data, chain_path):
    """Lift the CHR/POS coordinates of data in python with a chain file (see :func:`lift_positions`)."""
    print("Lifting coordinates in python...")
    contigs, positions, strands = lift_positions(data["CHR"].to_numpy(), data["POS"].to_numpy(), chain_path)
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
