import os, subprocess
import pandas as pd
import json
import wget
import shutil
import tarfile

from .constants import REF_PANELS

config_path = os.path.join(os.path.expanduser("~/.genal/"), "config.json")
# default_ref_path = os.path.join(os.getcwd(), "tmp_GENAL", "Reference_files")
default_ref_path = os.path.join(os.path.expanduser("~/.genal/"), "Reference_files")


def default_config():
    """Returns default config values"""
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
    default_config = {
        "paths": {
            "plink19_path": "",
            "liftover_path": "",
            "geno_path": "",
            "ref_path": default_ref_path,
        }
    }
    return default_config


def read_config():
    """Get config file data"""
    with open(config_path, "r") as f:
        config = json.load(f)
    return config


def write_config(config):
    """Write data to config file"""
    with open(config_path, "w") as f:
        json.dump(config, f, indent=4)
    return

def setup_genetic_path(path):
    """Configure the genetic data path based on user input and saved configuration."""
    config = read_config()
    if path is None:
        if not "geno_path" in config["paths"]:
            raise TypeError("No path has been saved in the config file. Please provide one.")
        path = config["paths"]["geno_path"]
        print(f"Using path saved in config file: {path}")
        return path

    # Ensure correct path format
    if path.count("$") > 1:
        raise TypeError(
            "The path should contain at most 1 '$' sign. Use it to indicate the chromosome number if the data is split by chromosomes."
        )
    # Check if the path is valid
    path = os.path.splitext(path)[0]  # Remove file extension
    if path.count("$") == 0:
        if not check_bfiles(path):
            raise TypeError("The path does not lead to valid bed/bim/fam files.")
    else:
        check_split = [check_bfiles(path.replace("$", str(i))) for i in range(1, 23, 1)]
        if not any(check_split):
            raise TypeError("The path does not lead to valid bed/bim/fam files.")

    # Save path to config.json file
    config["paths"]["geno_path"] = path
    write_config(config)
    return path

def create_tmp():
    """Create the temporary folder if not present"""
    if not os.path.exists("tmp_GENAL"):
        try:
            os.makedirs("tmp_GENAL")
        except OSError:
            raise OSError(
                "Unable to create the 'tmp_GENAL' directory. Check permissions."
            )


def set_reference_folder(path=""):
    """
    Set a folder path to store reference data.

    This function allows users to specify a directory where reference data will be stored.
    If the directory doesn't exist, it will be created. If no path is provided, a default
    directory named 'tmp_GENAL' in the current working directory will be used.

    Parameters:
        path (str, optional): The desired directory path for storing reference data. Defaults to a temporary folder in the current working directory.

    Raises:
        OSError: If the directory cannot be created.

    Returns:
        None: The function prints messages to inform the user of the status and any errors.
    """

    # If no path is provided, set default path to 'tmp_GENAL' in the current directory
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


def get_reference_panel_path(reference_panel="eur"):
    """
    Retrieve the path of the specified reference panel.

    This function checks if the provided reference panel is a valid path to bed/bim/fam files.
    If not, it checks if the reference panel exists in the reference folder. If it doesn't exist,
    the function attempts to download it.

    Parameters:
        reference_panel (str, optional): The name of the reference panel or a path to bed/bim/fam files. Defaults to "eur".

    Raises:
        ValueError: If the provided reference panel is not recognized.
        OSError: If there's an issue creating the directory.
        FileNotFoundError: If the reference panel is not found.

    Returns:
        str: The path to the reference panel.
    """

    # Remove file extension and check if it's a path to a bed/bim/fam triple
    reference_panel = os.path.splitext(reference_panel)[0]
    if check_bfiles(reference_panel):
        ref_panel_path = reference_panel
        print(f"Using the provided path as the reference panel.")

    else:
        # If it's not a valid path, check if the reference panel is recognized
        reference_panel = reference_panel.lower()
        config = read_config()
        if reference_panel not in REF_PANELS:
            raise ValueError(
                f"The reference_panel argument can only take values in {REF_PANELS} or be a valid path to bed/bim/fam files."
            )

        ref_path = config["paths"]["ref_path"]

        # Create the reference path if it doesn't exist
        if not os.path.exists(ref_path):
            try:
                os.makedirs(ref_path)
            except OSError:
                raise OSError(
                    "Unable to create the 'tmp_GENAL' directory. Check permissions."
                )

        ref_panel_name = reference_panel.upper()
        ref_panel_path = os.path.join(ref_path, ref_panel_name)

        # If the reference panel files don't exist, attempt to download them
        if not check_bfiles(ref_panel_path):
            print(
                f"The {reference_panel.capitalize()} reference panel was not found. Attempting to download it..."
            )
            print(
                "If you have already downloaded it, use set_reference_folder(path) to avoid downloading again."
            )
            url = f"https://storage.googleapis.com/genal_files/1kg.v3.tgz"
            try:
                wget.download(url, out=os.path.join(ref_path, "1kg.v3.tgz"))
            except Exception as e:
                print(f"Download unsuccessful: {e}")
                print(
                    "Manually download the reference file and use set_reference_folder(path)."
                )
                raise FileNotFoundError(f"Reference panel {reference_panel} not found.")

            print("Download successful. Decompressing...")
            with tarfile.open(os.path.join(ref_path, "1kg.v3.tgz"), "r:gz") as tar_ref:
                tar_ref.extractall(ref_path)
        else:
            print(f"Using the {ref_panel_name} reference panel.")

    return ref_panel_path


## Need to do the multi option
def load_reference_panel(reference_panel="eur"):
    """Load the bim file from the reference panel specified."""
    
    # Check if it's a path to a .bim file
    reference_panel = os.path.splitext(reference_panel)[0]
    if os.path.exists(reference_panel + ".bim"):
        ref_panel_path = reference_panel
        print(f"Using the provided bim file as the reference dataset.")
        
    # Else, check if it's one of the reference datasets names and get the path
    else:
        reference_panel = reference_panel.lower()
        if reference_panel == "multi":
            raise ValueError("Multi reference dataset not implemented yet.") 
        else:
            ref_panel_path = get_reference_panel_path(reference_panel)
            
    #Load it and return it
    reference_panel_df = pd.read_csv(
        ref_panel_path + ".bim", sep ="\t", names=["CHR","SNP","F","POS","A1","A2"]
    )
    return reference_panel_df

def set_plink(path=""):
    """Set the plink 1.9 path and verify that it is the correct version."""
    if not path:
        raise TypeError("You need to provide a path.")

    if not os.path.isfile(path):
        path = os.path.join(path, "plink")

    try:
        process = subprocess.run(
            [path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5, text=True
        )
        if not process.stdout.startswith("PLINK v1.9"):
            raise TypeError(
                "The path provided is an executable, but not the plink 1.9 executable. Check the path and plink version."
            )
    except Exception as e:
        raise TypeError(e)

    # Change config file
    config = read_config()
    config["paths"]["plink19_path"] = path
    write_config(config)

    print(f"Path to plink 1.9 successfully set: '{path}'")
    return


def get_plink19_path():
    """Return the plink19 path if it exists in the config file."""
    config = read_config()
    if not config["paths"]["plink19_path"]:
        raise ValueError(
            "The path to plink 1.9 has not been set yet. Use set_plink(path_to_plink) first."
        )
    else:
        return config["paths"]["plink19_path"]


def check_bfiles(filepath):
    """Check if the path specified leads to a bed/bim/fam triple."""
    if (
        os.path.exists("{}.bed".format(filepath))
        and os.path.exists("{}.bim".format(filepath))
        and os.path.exists("{}.fam".format(filepath))
    ):
        return True
    return False


def delete_tmp():
    """Delete the tmp folder."""
    if os.path.isdir("tmp_GENAL"):
        shutil.rmtree("tmp_GENAL")
        print("The tmp_GENAL folder has been successfully deleted.")
    else:
        print("There is no tmp_GENAL folder to delete in the current directory.")
    return