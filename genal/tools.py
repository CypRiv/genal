import os, subprocess, sys
import pandas as pd
import json
import wget
import shutil
import tarfile
import platform
import requests
import zipfile

from .constants import REF_PANELS, REF_PANELS_URL, CONFIG_DIR, BUILDS, REF_PARQUET_URL

config_path = os.path.join(CONFIG_DIR, "config.json")
default_ref_path = os.path.join(CONFIG_DIR, "Reference_files")


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
        filetype = config["paths"].get("geno_filetype", "bed")  # Default to bed if not specified
        print(f"Using path saved in config file: {path}")
        return path, filetype

    # Ensure correct path format
    if path.count("$") > 1:
        raise TypeError(
            "The path should contain at most 1 '$' sign. Use it to indicate the chromosome number if the data is split by chromosomes."
        )
    # Check if the path is valid
    path = os.path.splitext(path)[0]  # Remove file extension
    if path.count("$") == 0:
        is_bed = check_bfiles(path)
        is_pgen = check_pfiles(path)
        if not (is_bed or is_pgen):
            raise TypeError("The path does not lead to valid bed/bim/fam or pgen/pvar/psam files.")
        filetype = "bed" if is_bed else "pgen"
    else:
        check_split_bed = [check_bfiles(path.replace("$", str(i))) for i in range(1, 23, 1)]
        check_split_pgen = [check_pfiles(path.replace("$", str(i))) for i in range(1, 23, 1)]
        if not (any(check_split_bed) or any(check_split_pgen)):
            raise TypeError("The path does not lead to valid bed/bim/fam or pgen/pvar/psam files.")
        filetype = "bed" if any(check_split_bed) else "pgen"

    # Save path and filetype to config.json file
    config["paths"]["geno_path"] = path
    config["paths"]["geno_filetype"] = filetype
    write_config(config)
    return path, filetype

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
    """Delete the tmp folder."""
    if os.path.isdir("tmp_GENAL"):
        shutil.rmtree("tmp_GENAL")
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
    If not, it checks if the reference panel exists in the reference folder. If it doesn't exist,
    the function attempts to download it.
    
    Args:
        reference_panel (str): Reference panel identifier (e.g., "EUR_38", "AFR_37") 
            or path to custom reference panel (bed/bim/fam or pgen/pvar/psam files)
    
    Returns:
        tuple: (path to reference panel, filetype ['bed' or 'pgen'])
    """

    # Remove file extension and check if it's a path to a bed/bim/fam or pgen/pvar/psam triple
    reference_panel = os.path.splitext(reference_panel)[0]
    if check_bfiles(reference_panel):
        print(f"Using the provided path as the reference panel (bed format).")
        return reference_panel, "bed"
    elif check_pfiles(reference_panel):
        print(f"Using the provided path as the reference panel (pgen format).")
        return reference_panel, "pgen"
    else:
        # Standardize panel name
        reference_panel = reference_panel.upper()
        
        # Validate reference panel name
        if reference_panel not in REF_PANELS:
            raise ValueError(
                f"Invalid reference panel: {reference_panel}. Must be one of: {', '.join(REF_PANELS)} "
                "or a path to a valid reference panel in bed/bim/fam or pgen/pvar/psam format."
            )
        
        # Get config and set paths
        config = read_config()
        ref_path = config["paths"]["ref_path"]
        panel_dir = os.path.join(ref_path, reference_panel)
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
            wget.download(url, tar_path)
            print(f"\nExtracting {reference_panel}...")
            with tarfile.open(tar_path) as tar:
                tar.extractall(panel_dir)
            os.remove(tar_path)
        except Exception as e:
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
    reference_panel = os.path.splitext(reference_panel)[0]
    if os.path.exists(reference_panel + ".bim"):
        print(f"Using the provided bim file as the reference dataset.")
        reference_panel_df = pd.read_csv(
            reference_panel + ".bim", sep="\t", names=["CHR","SNP","F","POS","A1","A2"]
        )
        reference_panel_df = check_reference_panel(reference_panel_df)

    elif os.path.exists(reference_panel + ".pvar"):
        print(f"Using the provided pvar file as the reference dataset.")
        reference_panel_df = pd.read_csv(
            reference_panel + ".pvar", sep="\t", comment="#",
            names=["CHR","POS","SNP","A1","A2"]
        )
        reference_panel_df = check_reference_panel(reference_panel_df)
        
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
        config = read_config()
        ref_path = config["paths"]["ref_path"]
        variants_dir = os.path.join(ref_path, "reference_variants")
        variants_file = os.path.join(variants_dir, f"reference_variants_{build}.parquet")
        
        # Download if not exists
        if not os.path.exists(variants_file):
            print(f"Downloading reference variants for build {build} in {variants_dir} ...")
            os.makedirs(variants_dir, exist_ok=True)
            
            # Download parquet file
            url = REF_PARQUET_URL.format(build=build)
            try:
                wget.download(url, variants_file)
                print("\nDownload complete.")
            except Exception as e:
                if os.path.exists(variants_file):
                    os.remove(variants_file)
                raise RuntimeError(f"Failed to download reference variants: {e}")
        
        # Load parquet file
        try:
            reference_panel_df = pd.read_parquet(variants_file, engine="fastparquet")
            print(f"Using reference variants from build {build}.")
        except Exception as e:
            raise RuntimeError(f"Failed to load reference variants: {e}")
        
    
    return reference_panel_df

def check_reference_panel(df):
    """
    Check and standardize the reference panel DataFrame.
    
    Args:
        reference_panel_df (pd.DataFrame): Reference panel DataFrame with required columns 
            ["CHR", "POS", "SNP", "A1", "A2"]
            
    Returns:
        pd.DataFrame: Processed reference panel with standardized columns
        
    Raises:
        ValueError: If required columns are missing or if data types are invalid
    """
    
    # Convert allele columns to uppercase strings
    for col in ["A1", "A2"]:
        df[col] = df[col].astype(str).str.upper()
        
    # Convert CHR to string and remove 'chr' prefix if present
    if str(df["CHR"][0]).startswith("chr"):
        df["CHR"] = df["CHR"].astype(str).str.replace("^chr", "", regex=True)
    # Convert numeric values to int
    try:
        df["CHR"] = df["CHR"].astype(int)
    except ValueError:
        raise ValueError("Chromosome (CHR) column of the reference panel must contain integer values")
    
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
    """Set the plink 2.0 path and verify that it is the correct version."""
    if not path:
        raise TypeError("You need to provide a path.")

    if not os.path.isfile(path):
        path = os.path.join(path, "plink2")

    try:
        process = subprocess.run(
            [path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5, text=True
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
    """Return the plink2 path if it exists in the config file."""
    config = read_config()
    if not config["paths"]["plink2_path"]:
        raise ValueError(
            "The path to plink 2.0 has not been set yet. Use genal.set_plink(path_to_plink) first or run genal.install_plink() to install it automatically."
        )
    else:
        return config["paths"]["plink2_path"]


def install_plink(path=None):
    """Install plink 2.0 for the current operating system."""
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
        command (str): The PLINK command to execute
        ram_param_name (str): Name of the RAM parameter in the calling method
        class_name (str): Name of the class containing the RAM parameter
        
    Raises:
        RuntimeError: If PLINK command fails with detailed error message
    """
    try:
        output = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        if output.returncode != 0:
            raise subprocess.CalledProcessError(output.returncode, command, output.stdout, output.stderr)
    except subprocess.CalledProcessError as e:
        if "Out of memory" in e.stderr:
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
            raise ValueError("PLINK command failed. Check the error messages above for details.")