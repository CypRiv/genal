from pydantic import BaseModel
import pandas as pd

from .constants import REQUIRED_GWAS_COLUMNS


class GwasData(BaseModel):
    data = pd.DataFrame

    def validate_data(self):
        # Warn if essential columns are missing
        for column in REQUIRED_GWAS_COLUMNS:
            if not (column in self.data.columns):
                print(
                    f"Warning: the data doesn't include a {column} column. This may become an issue later on."
                )


class ClumpedData(BaseModel):
    data = pd.DataFrame
    kb = int = 250
    r2 = float = 0.1
    p1 = float = 5e-8
    p2 = float = 0.01
    reference = str = "eur"
