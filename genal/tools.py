import os, subprocess
import pandas as pd
import json
import wget
import tarfile

config_path = os.path.join(os.path.expanduser("~/.genal/"), "config.json")

def default_config():
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
    tmp_path = os.path.join(os.getcwd(), "tmp_GENAL", "Reference_files")
    default_config = {
        "paths": {
            "plink19_path": "/gpfs/gibbs/pi/falcone/00_General/software/plink1.9/plink",
            "liftover_path": "",
            "geno_path": "/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Resources/UKB_geno_files/",
            "ref_path": tmp_path
        },
        "ref_panels": ["eur","sas","eas","amr"]
    }
    return default_config

def read_config():
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config
        
def write_config(config):
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=4)
    return
        
def set_reference_folder(path=""):
    """
    Choose a folder path that will be used to store reference data to avoid downloading them multiple times.
    path="tmp_GENAL" By default, the path is set to the temporary folder created in the local directory
    """
    if not path:
        path = os.path.join(os.getcwd(), "tmp_GENAL", "Reference_files")
        print ("No path provided, defaulting to the local tmp_GENAL directory.")        
    if not os.path.isdir(path):
        os.makedirs(path)
        print(f"Creating the '{path}' directory.")
    if not os.access(path, os.R_OK):
        print(f"Error: The directory '{path}' is not readable.")
        return False
    # Check write permission
    if not os.access(path, os.W_OK):
        print(f"Error: The directory '{path}' is not writable.")
        return False        
    # Change config file
    config = read_config()
    config["paths"]["ref_path"] = path
    write_config (config)   
    print(f"Reference files will be downloaded and stored in: '{path}'")
    return

def get_reference_panel_path(reference_panel="eur"):
    """
    If the reference_panel is a path, checks that it leads to valid bed/bim/fam files.
    Otherwise, check if the given reference panel exists in the reference folder. 
    If not, attempts to download it.
    Returns its path.
    """
    #Check if it's a path to a bed/bim/fam triple
    reference_panel = os.path.splitext(reference_panel)[0]
    if check_bfiles(reference_panel):
        ref_panel_path = reference_panel
        print(f"Using the provided path as the reference panel.")
    else:
        reference_panel = reference_panel.lower()
        config = read_config()
        if reference_panel not in config['ref_panels']:
            raise ValueError(f"The reference_panel argument can only take values in {config['ref_panels']} depending on the reference panel to be used. Or it should be a valid path to bed/bim/fam files.")    
        ref_path = config["paths"]["ref_path"]
        if not os.path.exists(ref_path):
            os.makedirs(ref_path)
        ref_panel_name = reference_panel.upper()
        ref_panel_path = os.path.join(ref_path, ref_panel_name)
        if not check_bfiles(ref_panel_path):
            print (f"The {reference_panel.capitalize()} reference panel was not found. Attempting to download it... ")
            print ("If you have already downloaded it, you can specify its directory with set_reference_folder(path) to avoid downloading it again.")
            url = f"https://storage.googleapis.com/genal_files/1kg.v3.tgz"
            try: #Attempts to download
                wget.download(url, out=os.path.join(ref_path,"1kg.v3.tgz"))
            except Exception as e:
                print(f"The download was unsuccessful: {e}")
                print(f"You can manually download the reference file at xxx and specify its directory with the set_reference_folder(path) function.")
                raise FileNotFoundError(f"Reference panel {reference_panel} not found.")
            print(f"The download was successful. Decompressing...")
            with tarfile.open(os.path.join(ref_path, "1kg.v3.tgz"), "r:gz") as tar_ref:
                tar_ref.extractall(ref_path)
        else:
            print(f"Using the {ref_panel_name} reference panel.")
    return ref_panel_path

def load_reference_panel(reference_panel="eur"):
    """
    % Need to do multi
    """
    #Check if it's a path to a .bim file
    reference_panel = os.path.splitext(reference_panel)[0]
    if os.path.exists(reference_panel+".bim"):
        ref_panel_path = reference_panel
        print(f"Using the provided bim file as the reference dataset.")
    #If it's one of the reference datasets names
    else:
        reference_panel = reference_panel.lower()
        if reference_panel=="multi":
            raise ValueError("Multi reference dataset not implemented yet.") 
        else:
            ref_panel_path = get_reference_panel_path(reference_panel)
    reference_panel_df = pd.read_csv(ref_panel_path + ".bim", sep ="\t", names=["CHR","SNP","F","POS","A1","A2"])
    return reference_panel_df

def set_plink(path=""):
    """
    Set the plink 1.9 path.
    """
    if not path:
        raise TypeError("You need to provide a path.")

    if not os.path.isfile(path):
        path = os.path.join(path, "plink19")
        
    try:
        process = subprocess.run([path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5, text=True)
        if not process.stdout.startswith("PLINK v1.9"):
            raise TypeError("The path provided is an executable, but not the plink 1.9 executable. Check the path and plink version.")
    except Exception as e:
        raise TypeError(e)
 
    # Change config file
    config = read_config()
    config["paths"]["plink19_path"] = path
    write_config (config)
    
    print(f"Path to plink 1.9 successfully set: '{path}'")
    return

def get_plink19_path():
    """
    Return the plink19 path if it exists in the config file.
    """
    config = read_config()
    if not config["paths"]["plink19_path"]:
        raise ValueError("The path to plink 1.9 has not been set yet. Please use set_plink(path_to_plink) first.")
    else:
        return config["paths"]["plink19_path"]

def check_bfiles(filepath):
    if os.path.exists("{}.bed".format(filepath)) and os.path.exists("{}.bim".format(filepath)) and os.path.exists("{}.fam".format(filepath)):
        return True
    return False
