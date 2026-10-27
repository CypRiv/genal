import os

STANDARD_COLUMNS = ["CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P"]
CHR_CODES = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}
NULL_SNP_IDS = {"", ".", "NA", "NAN", "NONE"}
# Values read as missing by pandas.read_csv (default na_values), used to read PLINK outputs with pyarrow
READ_CSV_NA_VALUES = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>",
                      "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null"]
BUILDS = ["37", "38"]
POPULATIONS = ["EUR", "AFR", "EAS", "AMR", "SAS"]
REF_PANELS = [f"{pop}_{build}" for pop in POPULATIONS for build in BUILDS]
REF_PANEL_COLUMNS = ["CHR", "SNP", "POS", "A1", "A2"]
REF_DATA_VERSION = "v2"
BUCKET_URL = f"https://storage.googleapis.com/genal_files/{REF_DATA_VERSION}/"
REF_PANELS_URL = BUCKET_URL + "{panel}.tar.gz"
REF_PARQUET_URL = BUCKET_URL + "reference_variants_{build}.parquet"
GENE_INFO_URL = BUCKET_URL + "gene_id_mapping_filtered.parquet"
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
    "Wald": "Wald ratio",
}