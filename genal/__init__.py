from .GENO import *
from .MR import *
from .MR_tools import *
from .MRpresso import *
from .proxy import *
from .clump import *
from .lift import *
from .tools import *
from .geno_tools import *
from .association import *

import os
import json

config_dir = os.path.expanduser("~/.genal/") # Pas oublier de changer le config_path dans tools.py
config_path = os.path.join(config_dir, "config.json")

if not os.path.exists(config_dir):
    os.makedirs(config_dir)
    
if not os.path.exists(config_path):
    default_config = default_config()
    with open(config_path, "w") as f:
        json.dump(default_config, f)
    print (f"Configuration file for genal placed at '{config_path}'")