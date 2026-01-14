from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


@pytest.mark.slow
@pytest.mark.plink
def test_clump_data_plink2_smoke_if_plink_available(
    reset_genal_config, studyA_genetic_dir: Path
) -> None:
    """
    Optional smoke test: exercises the PLINK2 integration if `plink2` is available.

    This stays skipped in environments without PLINK2 installed.
    """
    plink2 = shutil.which("plink2")
    if not plink2:
        pytest.skip("plink2 not available in PATH")

    from genal.tools import read_config, write_config, create_tmp, load_reference_panel
    from genal.clump import clump_data_plink2

    cfg = read_config()
    cfg["paths"]["plink2_path"] = plink2
    write_config(cfg)

    create_tmp()

    ref_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(ref_prefix)
    snps = ref["SNP"].astype(str).head(50)

    # Use permissive thresholds so this reliably produces at least one clump.
    df = pd.DataFrame({"SNP": snps, "P": np.linspace(1e-10, 0.9, len(snps))})
    out = clump_data_plink2(
        df,
        reference_panel=ref_prefix,
        kb=250,
        r2=0.5,
        p1=1,
        p2=1,
        ram=2000,
        name="pytest_plink_smoke",
    )

    assert out is not None
    assert out.shape[0] >= 1
    assert set(out["SNP"]).issubset(set(df["SNP"]))
