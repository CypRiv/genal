from pyliftover import LiftOver
import os
import numpy as np
import wget
import gzip
import shutil
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from .tools import *

def lift_data(data, start="hg19",end="hg38", replace=False, extraction_file=False, chain_file="", name=""):
    """Perform a liftover from a genetic build to another.
    If the chain file required to do the liftover is not present, download it first. It is also possible to manually provide the path to the chain file.
    If the number of SNPs to be lifted is large (> a few millions), it is faster to use the lift_data_liftover function.
    start="hg19": current build of the data
    end="hg38": build to be lifted to
    replace=False: whether to change the dataframe inplace or return a new one
    extraction_file==True, also print a CHR POS SNP space delimited file for extraction in All of Us (WES data)
    chain_file="": path to a local chain file to be used for the lift. If provided, the start and end arguments are not considered.
    name="": can be used to specify a filename or filepath (without extension) to save the lifted dataframe. If not provided, will be saved in the current folder as [name]_lifted.txt
    """
    ##Check that the mandatory columns are present, make sure the type is right, and create the tmp folder
    for column in ["CHR","POS"]:
        if not(column in data.columns):
            raise ValueError("The column {column} is not found in the data!".format(column=column))

    if not os.path.exists("tmp_GENAL"):
        os.makedirs("tmp_GENAL")     

    # Loading the appropriate chain file.
    if chain_file!="": #If user provides path to a local chain path
        print("You provided a path to a local chain path. This will be used for the lift.")
        if not os.path.isfile(chain_file):
            raise ValueError("The path you provided does not lead to a file.")
        else:
            lo = LiftOver(chain_file)
    else: #Interpret the start and end build
        chain_name=f"{start.lower()}To{end.capitalize()}.over.chain"
        config = read_config()
        ref_path = config["paths"]["ref_path"]
        chains_folder_path = os.path.join(ref_path, "chain_files")
        if not os.path.exists(chains_folder_path): #Create the directory if necessary
            os.makedirs(chains_folder_path)
        chain_path = os.path.join(chains_folder_path, chain_name)
        if os.path.isfile(chain_path): #If file already exists
            print(f"Found chain file to lift from {start} to {end}.")
        else: #If not: attempts to download it
            print(f"The chain file to lift from {start} to {end} was not found. Attempting to download it to {chains_folder_path} (You can change this with set_reference_folder(path).)")
            url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{start.lower()}/liftOver/{chain_name}.gz"
            try: #Attempts to download
                wget.download(url,out=chains_folder_path) 
            except Exception as e:
                print(f"The download was unsuccessful: {e}")
                print (f"You can manually download the appropriate chain file at http://hgdownload.cse.ucsc.edu/downloads.html#liftover and specify its local path with the chain_file argument.")
                raise FileNotFoundError(f"Chain file not found.")
            print(f"The download was successful. Unzipping...")
            with gzip.open(f'{chain_path}.gz', 'rb') as f_in:
                with open(chain_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        lo = LiftOver(chain_path)

    #Deleting NAs in CHR and POS columns
    nrows=data.shape[0]
    data.dropna(subset=["CHR","POS"],inplace=True)
    data.reset_index(drop=True,inplace=True)
    n_na=nrows-data.dropna().shape[0]
    if n_na>0:
        print(f"Excluding {n_na}({n_na/nrows*100:.3f}%) SNPs containing NaN values in columns CHR or POS.")

    #Lifting
    if nrows>1000000:
        print("Your data is large, this can take a few minutes...")
    
   # print("top")
   # def convert_coordinate_wrapper(i):
   #     return lo.convert_coordinate(f'chr{data.loc[i,"CHR"]}', data.loc[i,"POS"], '-')
   # with ThreadPoolExecutor() as executor: 
   #     results = list(tqdm(executor.map(convert_coordinate_wrapper, range(len(data))), total=len(data)))
    
    
    def convert_coordinate_wrapper(args):
        return lo.convert_coordinate(f'chr{args[0]}', args[1], '-')
    args = data[['CHR', 'POS']].to_records(index=False)
    with ThreadPoolExecutor() as executor: 
        results = list(executor.map(convert_coordinate_wrapper, args))
    
    data['POS'] = [res[0][1] if res else np.nan for res in results]
    nrows=data.shape[0]
    data.dropna(subset=["POS"],inplace=True)
    data['POS'] = data['POS'].astype("Int32")
    data.reset_index(drop=True,inplace=True)
    n_na=nrows-data.dropna().shape[0]
    if n_na > 0:
        print(f"Lift completed. {n_na}({n_na/nrows*100:.3f}%) SNPs could not been lifted.")
    else:
        print("Lift completed. All the SNPs have been lifted.")

    ## Save files: whole data and also the extraction file if extraction_file=True
    data.to_csv(f"{name + '_lifted'}.txt",sep="\t",header=True,index=False)
    print(f"Lifted list of SNPs saved to {name + '_lifted'}.txt")
    if extraction_file:
        if not("SNP" in data.columns):
            data["SNP"]=data.CHR.astype(str)+":"+data.POS.astype(str)
        data[["CHR","POS","SNP"]].to_csv(f"{name+ '_lifted'}_extraction.txt", sep=" ", header=False, index=False)
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