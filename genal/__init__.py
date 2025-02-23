import os
import json
from .tools import default_config, write_config, set_plink, install_plink, delete_tmp, get_reference_panel_path, get_plink_path
from .geno_tools import Combine_Geno
from .constants import CONFIG_DIR

__version__ = "1.2.9"

config_path = os.path.join(CONFIG_DIR, "config.json")

if not os.path.exists(CONFIG_DIR):
    os.makedirs(CONFIG_DIR)


if not os.path.exists(config_path):
    write_config(default_config())
    print(f"Configuration file for genal placed at '{config_path}'")

from .Geno import Geno
