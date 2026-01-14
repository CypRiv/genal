from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _set_plink2_path(reset_genal_config, tests_dir: Path) -> str:
    from genal.tools import read_config, write_config

    plink2 = shutil.which("plink2")
    if not plink2:
        candidate = tests_dir / ".genal_test_home" / ".genal" / "plink2" / "plink2"
        if candidate.exists():
            plink2 = str(candidate)
    if not plink2:
        pytest.skip("plink2 not available (needed for clumping integration tests)")

    cfg = read_config()
    cfg["paths"]["plink2_path"] = plink2
    write_config(cfg)
    return plink2


def _link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        dst.symlink_to(src)
    except OSError:
        shutil.copyfile(src, dst)


@pytest.mark.slow
@pytest.mark.plink
def test_clump_data_plink2_basic(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.clump import clump_data_plink2
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    ref_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(ref_prefix)
    snps = ref["SNP"].astype(str).head(100)

    df = pd.DataFrame({"SNP": snps, "P": np.linspace(1e-10, 0.9, len(snps))})
    out = clump_data_plink2(
        df,
        reference_panel=ref_prefix,
        kb=250,
        r2=0.5,
        p1=1,
        p2=1,
        ram=2000,
        name="pytest_clump_basic",
    )
    assert out is not None
    assert out.shape[0] >= 1
    assert set(out["SNP"]).issubset(set(df["SNP"]))


@pytest.mark.slow
@pytest.mark.plink
def test_clump_data_plink2_no_significant_results(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.clump import clump_data_plink2
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    ref_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(ref_prefix)
    snps = ref["SNP"].astype(str).head(50)

    df = pd.DataFrame({"SNP": snps, "P": np.full(len(snps), 1e-6)})
    out = clump_data_plink2(
        df,
        reference_panel=ref_prefix,
        p1=1e-300,  # too stringent: no top variants
        p2=1e-300,
        ram=2000,
        name="pytest_clump_none",
    )
    assert out is None


@pytest.mark.slow
@pytest.mark.plink
def test_clump_data_plink2_custom_parameters(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.clump import clump_data_plink2
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    ref_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(ref_prefix)
    snps = ref["SNP"].astype(str).head(100)

    df = pd.DataFrame({"SNP": snps, "P": np.linspace(1e-20, 1e-4, len(snps))})
    out = clump_data_plink2(
        df,
        reference_panel=ref_prefix,
        kb=500,
        r2=0.2,
        p1=1e-4,
        p2=1e-4,
        ram=2000,
        name="pytest_clump_custom",
    )
    assert out is not None
    assert out.shape[0] >= 1


@pytest.mark.slow
@pytest.mark.plink
def test_clump_data_plink2_bed_vs_pgen(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path, tmp_path: Path
) -> None:
    from genal.clump import clump_data_plink2
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    bed_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(bed_prefix)
    snps = ref["SNP"].astype(str).head(100)
    df = pd.DataFrame({"SNP": snps, "P": np.linspace(1e-10, 0.9, len(snps))})

    out_bed = clump_data_plink2(
        df,
        reference_panel=bed_prefix,
        p1=1,
        p2=1,
        ram=2000,
        name="pytest_clump_bed",
    )
    assert out_bed is not None

    # Make a pgen-only copy to force `get_reference_panel_path()` to choose pgen format.
    pgen_dir = tmp_path / "pgen_only"
    pgen_prefix = pgen_dir / "panel9"
    src = studyA_genetic_dir / "plink_sampled_9"
    for ext in [".pgen", ".pvar", ".psam"]:
        _link_or_copy(src.with_suffix(ext), pgen_prefix.with_suffix(ext))

    out_pgen = clump_data_plink2(
        df,
        reference_panel=str(pgen_prefix),
        p1=1,
        p2=1,
        ram=2000,
        name="pytest_clump_pgen",
    )
    assert out_pgen is not None
