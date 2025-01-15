import os

STANDARD_COLUMNS = ["CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P"]
REF_PANELS = ["eur", "sas", "eas", "amr", "afr"]
REF_PANEL_COLUMNS = ["CHR", "SNP", "POS", "A1", "A2"]
REF_PANELS_URL = "https://storage.googleapis.com/genal_files/1kg.v3.tgz"
CONFIG_DIR = os.path.expanduser("~/.genal/")
CHECKS_DICT = {
    "CHR": False,
    "POS": False,
    "P": False,
    "EA": False,
    "NEA": False,
    "BETA": False,
    "SNP": False,
    "NA_removal": False,
}
MR_METHODS_NAMES = {
    "IVW": "Inverse-Variance Weighted",
    "IVW-RE": "Inverse Variance Weighted (Random Effects)",
    "IVW-FE": "Inverse Variance Weighted (Fixed Effects)",
    "UWR": "Unweighted Regression",
    "WM": "Weighted Median",
    "WM-pen": "Penalised Weighted Median",
    "Simple-median": "Simple Median",
    "Sign": "Sign concordance test",
    "Egger": ("MR Egger", "Egger Intercept"),
    "Egger-boot": ("MR Egger bootstrap", "Egger Intercept bootstrap"),
    "Simple-mode": "Simple mode",
    "Weighted-mode": "Weighted mode",
}