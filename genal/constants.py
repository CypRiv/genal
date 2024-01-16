STANDARD_COLUMNS = ["CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P"]
REF_PANELS = ["eur", "sas", "eas", "amr", "afr"]
REF_PANEL_COLUMNS = ["CHR", "SNP", "POS", "A1", "A2"]
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
    "Egger": ("MR Egger", "Egger Intercept"),
    "Egger-boot": ("MR Egger bootstrap", "Egger Intercept bootstrap"),
    "WM": "Weighted Median",
    "WM-pen": "Penalised Weighted Median",
    "Simple-median": "Simple Median",
    "IVW": "Inverse-Variance Weighted",
    "IVW-RE": "Inverse Variance Weighted (Random Effects)",
    "IVW-FE": "Inverse Variance Weighted (Fixed Effects)",
    "UWR": "Unweighted Regression",
    "Sign": "Sign concordance test",
}