from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest


def test_check_bfiles_and_pfiles(studyA_genetic_dir: Path) -> None:
    from genal.tools import check_bfiles, check_pfiles

    bed_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    missing_prefix = str(studyA_genetic_dir / "does_not_exist")

    assert check_bfiles(bed_prefix) is True
    assert check_pfiles(bed_prefix) is True  # StudyA fixtures include both formats
    assert check_bfiles(missing_prefix) is False
    assert check_pfiles(missing_prefix) is False


def test_setup_genetic_path_writes_config(reset_genal_config, studyA_genetic_dir: Path) -> None:
    from genal.tools import setup_genetic_path, read_config

    # Pass a path with extension; function should strip it.
    bed_file = str(studyA_genetic_dir / "plink_sampled_9.bed")
    path, filetype = setup_genetic_path(bed_file)

    assert filetype == "bed"
    assert path.endswith("plink_sampled_9")

    cfg = read_config()
    assert cfg["paths"]["geno_path"].endswith("plink_sampled_9")
    assert cfg["paths"]["geno_filetype"] == "bed"


def test_setup_genetic_path_split_placeholder(reset_genal_config, studyA_genetic_dir: Path) -> None:
    from genal.tools import setup_genetic_path

    path_with_dollar = str(studyA_genetic_dir / "plink_sampled_$")
    path, filetype = setup_genetic_path(path_with_dollar)

    assert path.endswith("plink_sampled_$")
    assert filetype in {"bed", "pgen"}


def test_create_and_delete_tmp_genal_folder() -> None:
    from genal.tools import create_tmp, delete_tmp

    create_tmp()
    assert Path("tmp_GENAL").is_dir()
    delete_tmp()
    assert not Path("tmp_GENAL").exists()


def test_get_reference_panel_path_local_prefix(studyA_genetic_dir: Path) -> None:
    from genal.tools import get_reference_panel_path

    prefix = str(studyA_genetic_dir / "plink_sampled_13")
    out_path, out_type = get_reference_panel_path(prefix)

    assert out_path == prefix
    assert out_type == "bed"  # if both are present, bed is preferred by get_reference_panel_path()


def test_load_reference_panel_from_bim(studyA_genetic_dir: Path) -> None:
    from genal.tools import load_reference_panel

    prefix = str(studyA_genetic_dir / "plink_sampled_9")
    df = load_reference_panel(prefix)

    assert {"CHR", "POS", "SNP", "A1", "A2"}.issubset(df.columns)
    assert df["CHR"].dtype.kind in {"i", "u"}
    assert df["POS"].dtype.kind in {"i", "u"}
    alleles = df[["A1", "A2"]].astype(str)
    assert (alleles.apply(lambda col: col.str.upper().eq(col))).all().all()


def test_run_plink_command_error_branches() -> None:
    from genal.tools import run_plink_command

    with pytest.raises(ValueError):
        cmd_fail = f"PYTHONDONTWRITEBYTECODE=1 {sys.executable} -c \"import sys; sys.exit(1)\""
        run_plink_command(cmd_fail)

    cmd_oom = (
        f"PYTHONDONTWRITEBYTECODE=1 {sys.executable} -c "
        "\"import sys; sys.stderr.write('Out of memory'); sys.exit(1)\""
    )
    with pytest.raises(RuntimeError):
        run_plink_command(cmd_oom)
