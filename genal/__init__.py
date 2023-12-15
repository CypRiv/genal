import os
import json
from .tools import default_config, write_config, set_plink
from .geno_tools import delete_tmp

config_dir = os.path.expanduser(
    "~/.genal/"
)  # Don't forget to change the config_path dans tools.py
config_path = os.path.join(config_dir, "config.json")

if not os.path.exists(config_dir):
    os.makedirs(config_dir)

if not os.path.exists(config_path):
    write_config(default_config())
    print(f"Configuration file for genal placed at '{config_path}'")

from .Geno import Geno
