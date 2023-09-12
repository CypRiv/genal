from pyliftover import LiftOver
import os
import numpy as np
import wget
import gzip
import shutil
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from .tools import *


def lift_data(data, start="hg19", end="hg38", extraction_file=False, chain_file=None, name=""):
    """
    Perform a liftover from one genetic build to another. If the chain file required for the liftover is not present, it will be downloaded. It's also possible to manually provide the path to the chain file. 
    If the dataset is large, it is suggested to use an alternate method (e.g., `lift_data_liftover`). 
    
    Args:
        data (pd.DataFrame): The input data containing at least "CHR" and "POS" columns.
        start (str, optional): The current build of the data. Defaults to "hg19".
        end (str, optional): The target build for liftover. Defaults to "hg38".
        extraction_file (bool, optional): If True, also prints a CHR POS SNP space-delimited file for extraction. Defaults to False.
        chain_file (str, optional): Path to a local chain file for the lift. Overrides the start and end arguments if provided.
        name (str, optional): Specify a filename or filepath (without extension) for saving. If not provided, uses [name]_lifted.txt.

    Raises:
        ValueError: If required columns are missing or if provided chain file path is incorrect.

    Returns:
        pd.DataFrame: Lifted data.
        
    Notes:
        Function for the :meth:`GENO.lift` method.
    """
    
    # Ensure mandatory columns are present in the input data
    for column in ["CHR", "POS"]:
        if column not in data.columns:
            raise ValueError(f"The column {column} is not found in the data!")

    # Create temporary directory if it doesn't exist
    if not os.path.exists("tmp_GENAL"):
        try:
            os.makedirs("tmp_GENAL")
        except OSError:
            raise OSError("Unable to create the 'tmp_GENAL' directory. Check permissions.")
    
    # Prepare chain file and get LiftOver object
    lo = prepare_chain_file(chain_file, start, end)
    
    # Perform liftover
    data = lift_coordinates(data, lo)
    
    # Handle post-liftover operations
    data = post_lift_operations(data, name, extraction_file)
    
    return data

def prepare_chain_file(chain_file, start, end):
    """Handle chain file loading, downloading if necessary. Return LiftOver object."""
    if chain_file is not None:  # If a local chain file is provided
        if not os.path.isfile(chain_file):
            raise ValueError("The provided path does not lead to a valid file.")
        print("You provided a path to a local chain path which will be used for the lift.")
        lo = LiftOver(chain_file)
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
                raise OSError("Unable to create the 'tmp_GENAL' directory. Check permissions.")

        # Check for the chain file locally or download it if necessary
        chain_path = os.path.join(chains_folder_path, chain_name)
        if not os.path.isfile(chain_path):
            # Download the chain file
            url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{start.lower()}/liftOver/{chain_name}.gz"
            try:
                wget.download(url, out=chains_folder_path)
                # Decompress the downloaded file
                print(f"The download was successful. Unzipping...")
                with gzip.open(f'{chain_path}.gz', 'rb') as f_in, open(chain_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            except Exception as e:
                print(f"The download was unsuccessful: {e}")
                print("Consider downloading the chain file manually from the UCSC website and providing its path via the chain_file argument.")
                raise FileNotFoundError("Chain file not found.")
        lo = LiftOver(chain_path)
    return lo

def lift_coordinates(data, lo):
    """Perform liftover on data using LiftOver object after handling missing values."""
    # Handle missing values in key columns
    nrows = data.shape[0]
    data.dropna(subset=["CHR", "POS"], inplace=True)
    data.reset_index(drop=True, inplace=True)
    n_na = nrows - data.shape[0]
    if n_na:
        print(f"Excluded {n_na} SNPs ({n_na/nrows*100:.3f}%) with NaN values in CHR or POS columns.")

    # Perform the lift
    def convert_coordinate(args):
        return lo.convert_coordinate(f'chr{args[0]}', args[1], '-')
    
    args = data[['CHR', 'POS']].to_records(index=False)
    results = list(ThreadPoolExecutor().map(convert_coordinate, args))
    
    data['POS'] = [res[0][1] if res else np.nan for res in results]
    nrows = data.shape[0]
    data.dropna(subset=["POS"], inplace=True)
    data['POS'] = data['POS'].astype("Int32")
    data.reset_index(drop=True, inplace=True)
    n_na = nrows - data.shape[0]
    if n_na:
        print(f"{n_na} SNPs ({n_na/nrows*100:.3f}%) could not be lifted.")
    else:
        print("All SNPs have been lifted successfully.")
    return data

def post_lift_operations(data, name, extraction_file):
    """Handle post-liftover operations like reporting, and saving results."""
    data.to_csv(f"{name + '_lifted'}.txt", sep="\t", header=True, index=False)
    print(f"Lifted list of SNPs saved to {name + '_lifted'}.txt")
    if extraction_file:
        if not("SNP" in data.columns):
            data["SNP"] = data["CHR"].astype(str) + ":" + data["POS"].astype(str)
        data[["CHR", "POS", "SNP"]].to_csv(f"{name + '_lifted'}_extraction.txt", sep=" ", header=False, index=False)
        print(f"Extraction file saved to {name+ '_lifted'}_extraction.txt")
    return data



### TO DO
def lift_data_liftover(data, start="hg19",end="hg38", replace=False, extraction_file=False, chain_file="", name="", genal_tools_path="", liftover_path=""):
    """Perform a liftover from a genetic build to another using the LiftOver software (requires Linux).
    The software must be installed on your system. It can be downloaded from https://genome-store.ucsc.edu
    If the chain file required to do the liftover is not present, download it first. It is also possible to manually provide the path to the chain file.
    start="hg19": current build of the data
    end="hg38": build to be lifted to
    replace=False: whether to change the dataframe inplace or return a new one
    extraction_file==True, also print a CHR POS SNP space delimited file for extraction in All of Us (WES data)
    chain_file="": path to a local chain file to be used for the lift. If provided, the start and end arguments are not considered.
    name="": can be used to specify a filename or filepath (without extension) to save the lifted dataframe. If not provided, will be saved in the current folder as [name]_lifted.txt
    genal_tools_path: path to the Genal_tools folder. If not provided, the tmp_GENAL folder will be used.
    liftover_path: path to LiftOver executable
        """
    ##Check that the mandatory columns are present, make sure the type is right, and create the tmp folder
    for column in ["CHR","POS"]:
        if not(column in self.data.columns):
            raise ValueError("The column {column} is not found in the data!".format(column=column))
    if not os.path.exists("tmp_GENAL"):
        os.makedirs("tmp_GENAL")

    ## Decide whether to lift the clumped data or the base data based on argument clumped and the presence of attribute data_clumped. Write the correct data in the format needed by liftOver.
    if clumped==False or not(hasattr(self,"data_clumped")):
        print("Lifting the unclumped data.")
        data=self.data.copy()
        data["CHR"]="chr"+data.CHR.astype(str)
        data[["CHR","POS","POS"]].to_csv(f"tmp_GENAL/{self.name}.prelift",sep=" ",index=False,header=False)
    elif clumped==True:
        print("Lifting the clumped data.")
        data=self.data_clumped.copy()
        data["CHR"]="chr"+data.CHR.astype(str)
        data[["CHR","POS","POS"]].to_csv(f"tmp_GENAL/{self.name}.prelift",sep=" ",index=False,header=False)


    ## Call the liftOver software.
    command=f'{liftover_path}liftOver tmp_GENAL/{self.name}.prelift \
    {liftover_path}hg19ToHg38.over.chain \
    tmp_GENAL/{self.name}.postlift tmp_GENAL/unMapped'
    output=subprocess.run(command, shell=True,capture_output=True,text=True,check=True)

    ## Read the output, print the number of unlifted SNPs and remove them from the prelift data. 
    df_post=pd.read_csv(f"tmp_GENAL/{self.name}.postlift",sep="\t",header=None)
    unMapped = open('tmp_GENAL/unMapped', 'r')
    Lines = unMapped.readlines()
    if len(Lines)>0:
        print(f"{int(len(Lines)/2)} SNPs could not been lifted.")
    indices=list()
    for i in range(1,len(Lines),2):
        c=Lines[i].strip()
        (chrom,pos,pos)=c.split("\t")
        indices.append(str(chrom)+":"+str(pos))
    data["SNP_IDS"]=data.CHR.astype(str)+":"+data.POS.astype(str)
    data=data[~(data.SNP_IDS.isin(indices))].drop(columns="SNP_IDS").reset_index(drop=True)

    ## Merge prelift and postlift data. Unknown chr from the output of liftOver are assigned the value 99. SNPs mapped to unknown chr are deleted from the final data and their number printed.
    data["POS"]=df_post[1].astype(int)
    data["CHR"]=df_post[0].str.split("chr",expand=True)[1].str.split("_",expand=True)[0].replace({"X":99,"Un":99}).astype(int)
    nrow_before=data.shape[0]
    data=data[data.CHR!=99]
    nrow_diff=nrow_before-data.shape[0]
    if nrow_diff>0:
        print(f"{nrow_diff} SNPs were lifted to an unknown chromosome and deleted from the final files.")

    ## Save files: whole data and also the extraction file to be used in All of Us if extraction_file=True
    data.to_csv(f"{self.name}_38.txt",sep="\t",header=True,index=False)
    if extraction_file:
        if not("SNP" in data.columns):
            data["SNP"]=data.CHR.astype(str)+":"+data.POS.astype(str)
        data[["CHR","POS","SNP"]].to_csv(f"{self.name}_38_extraction.txt",sep=" ",header=False,index=False)
    return data