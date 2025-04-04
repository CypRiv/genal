from pyliftover import LiftOver
import os, subprocess
import numpy as np
import wget
import gzip
import shutil
import uuid
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

from .tools import read_config, create_tmp


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
    If the dataset is large, it is suggested to use an alternate method (e.g., `lift_data_liftover`).

    Args:
        data (pd.DataFrame): The input data containing at least "CHR" and "POS" columns.
        start (str, optional): The current build of the data. Defaults to "hg19".
        end (str, optional): The target build for liftover. Defaults to "hg38".
        extraction_file (bool, optional): If True, also prints a CHR POS SNP space-delimited file for extraction. Defaults to False.
        chain_file (str, optional): Path to a local chain file for the lift. Overrides the start and end arguments if provided.
        name (str, optional): Specify a filename or filepath (without extension) for saving. If not provided, the data is not saved.
        liftover_path (str, optional): Specify the path to the USCS liftover executable. If not provided, the lift will be done in python (slower for large amount of SNPs).
        object_id (str, optional): Specify the object id for tmp file writing (internal use only)

    Raises:
        ValueError: If required columns are missing or if provided chain file path is incorrect.

    Returns:
        pd.DataFrame: Lifted data.

    Notes:
        Function for the :meth:`Geno.lift` method.
    """
    # Generate a default name if none is provided
    if name is None:
        name = str(uuid.uuid4())[:8]

    # Prepare chain file and get its path
    chain_path = prepare_chain_file(chain_file, start, end)

    # Prepare the data for lifting: handle missing values in CHR, POS columns
    nrows = data.shape[0]
    data.dropna(subset=["CHR", "POS"], inplace=True)
    data.reset_index(drop=True, inplace=True)
    n_na = nrows - data.shape[0]
    if n_na:
        print(
            f"Excluded {n_na} SNPs ({n_na/nrows*100:.3f}%) with NaN values in CHR or POS columns."
        )

    # Perform liftover with the liftover executable or in python
    if liftover_path:
        data = lift_coordinates_liftover(data, object_id, chain_path, liftover_path)
    else:
        data = lift_coordinates_python(data, chain_path)

    # Handle post-liftover operations
    data = post_lift_operations(data, name, extraction_file)

    return data


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
        # Construct chain filename
        chain_name = f"{start.lower()}To{end.capitalize()}.over.chain"
        config = read_config()
        ref_path = config["paths"]["ref_path"]
        chains_folder_path = os.path.join(ref_path, "chain_files")

        # Ensure directory for chain files exists
        if not os.path.exists(chains_folder_path):
            try:
                os.makedirs(chains_folder_path)
            except OSError:
                raise OSError(
                    "Unable to create the 'tmp_GENAL' directory. Check permissions."
                )

        # Check for the chain file locally or download it if necessary
        chain_path = os.path.join(chains_folder_path, chain_name)
        if not os.path.isfile(chain_path):
            print(
                f"The chain file to lift from {start} to {end} was not found. Attempting to download it..."
            )
            # Download the chain file
            url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{start.lower()}/liftOver/{chain_name}.gz"
            try:
                wget.download(url, out=chains_folder_path)
                # Decompress the downloaded file
                print(f"The download was successful. Unzipping...")
                with gzip.open(f"{chain_path}.gz", "rb") as f_in, open(
                    chain_path, "wb"
                ) as f_out:
                    shutil.copyfileobj(f_in, f_out)
            except Exception as e:
                print(f"The download was unsuccessful: {e}")
                print(
                    "Consider downloading the chain file manually from the UCSC website and providing its path via the chain_file argument."
                )
                raise FileNotFoundError("Chain file not found.")

    return chain_path


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
            timeout=5,
            text=True,
        )
        if not process.stderr.startswith("liftOver"):
            raise TypeError(
                "The path provided is an executable, but not the liftOver executable. Check the path."
            )
    except Exception as e:
        raise TypeError(e)
    print("Lifting coordinates using liftOver.")

    # Write data in correct format for liftOver
    create_tmp()
    data["CHR_liftover"] = "chr" + data.CHR.astype(str)
    to_lift_filename = os.path.join("tmp_GENAL", f"{object_id}.prelift")
    lifted_filename = os.path.join("tmp_GENAL", f"{object_id}.postlift")
    unmapped_filename = os.path.join("tmp_GENAL", f"{object_id}_unMapped")
    data[["CHR_liftover", "POS", "POS"]].to_csv(
        to_lift_filename, sep=" ", index=False, header=False
    )

    # Call the liftOver software
    command = f"{liftover_path} {to_lift_filename} \
    {chain_path} {lifted_filename} {unmapped_filename}"
    try:
        output = subprocess.run(
            command, shell=True, capture_output=True, text=True, check=True
        )
    except Exception as e:
        print(f"Error running liftOver: {e}")
        raise ValueError("Error running liftOver. Check error message for more details.")

    ## Read the output, print the number of unlifted SNPs and remove them from the prelift data.
    df_post = pd.read_csv(lifted_filename, sep="\t", header=None)
    unMapped = open(unmapped_filename, "r")
    Lines = unMapped.readlines()
    if len(Lines) > 0:
        print(f"{int(len(Lines)/2)} SNPs could not be lifted.")
    else:
        print(f"All SNPs have been lifted.")
    indices = list()
    for i in range(1, len(Lines), 2):
        c = Lines[i].strip()
        (chrom, pos, pos) = c.split("\t")
        indices.append(str(chrom) + ":" + str(pos))
    drop_indices = data[(data.CHR_liftover.astype(str) + ":" + data.POS.astype(str)).isin(indices)].index
    data.drop(index=drop_indices, inplace=True)
    data.reset_index(drop=True, inplace=True)

    # Check the length of files
    if len(data) != len(df_post):
        raise ValueError(
            "There was a problem lifting with liftOver. Try lifting in python (liftover_path = None)."
        )

    ## Merge prelift and postlift data. Unknown chr from the output of liftOver are assigned the value 99. SNPs mapped to unknown chr are deleted from the final data and their number printed.
    data["POS"] = df_post[1].astype(int)
    data["CHR"] = (
        df_post[0]
        .str.split("chr", expand=True)[1]
        .str.split("_", expand=True)[0]
        .replace({"X": 99, "Y": 99, "Un": 99})
        .astype(int)
    )
    nrow_before = data.shape[0]
    drop_chr_indices = data[data.CHR == 99].index
    data.drop(index=drop_chr_indices, inplace=True)
    nrow_diff = nrow_before - data.shape[0]
    if nrow_diff > 0:
        print(
            f"{nrow_diff} SNPs were lifted to an unknown chromosome and deleted from the final files."
        )
    data.drop(columns=["CHR_liftover"], inplace=True)
    return data


def lift_coordinates_python(data, chain_path):
    """Perform liftover on data using the chain passed."""
    lo = LiftOver(chain_path)

    # Print message
    print("Lifting coordinates in python...")
    nrows = data.shape[0]
    if nrows > 500000:
        print("Your data is large, this can take a few minutes...")

    # Perform the lift
    def convert_coordinate(args):
        return lo.convert_coordinate(f"chr{args[0]}", args[1], "-")

    args = data[["CHR", "POS"]].to_records(index=False)
    results = list(ThreadPoolExecutor().map(convert_coordinate, args))

    data["POS"] = [res[0][1] if res else np.nan for res in results]
    data["CHR"] = [res[0][0].split("chr")[1] if res else np.nan for res in results]
    nrows = data.shape[0]
    data.dropna(subset=["POS", "CHR"], inplace=True)
    data["POS"] = data["POS"].astype("Int32")
    data["CHR"] = data["CHR"].astype("Int32")
    data.reset_index(drop=True, inplace=True)
    n_na = nrows - data.shape[0]
    if n_na:
        print(f"{n_na} SNPs ({n_na/nrows*100:.3f}%) could not be lifted.")
    else:
        print("All SNPs have been lifted.")
    return data


def post_lift_operations(data, name, extraction_file):
    """Handle post-liftover operations like reporting, and saving results."""
    if name:
        filename = os.path.splitext(name)[0] + ".txt"
        data.to_csv(f"{filename}", sep="\t", header=True, index=False)
        print(f"Lifted list of SNPs saved to {filename}")
    if extraction_file:
        if not ("SNP" in data.columns):
            data["SNP"] = data["CHR"].astype(str) + ":" + data["POS"].astype(str)
        data[["CHR", "POS", "SNP"]].to_csv(
            f"{name + '_lifted'}_extraction.txt", sep=" ", header=False, index=False
        )
        print(f"Extraction file saved to {name+ '_lifted'}_extraction.txt")
    return data
