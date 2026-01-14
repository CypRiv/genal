from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _set_plink2_path(reset_genal_config, tests_dir: Path) -> str:
    """
    Ensure `genal.tools.get_plink_path()` works in tests.

    Preference order:
    1) `plink2` in PATH (CI)
    2) `tests/.genal_test_home/.genal/plink2/plink2` (local, isolated install)
    """
    from genal.tools import read_config, write_config

    plink2 = shutil.which("plink2")
    if not plink2:
        candidate = (
            tests_dir
            / ".genal_test_home"
            / ".genal"
            / "plink2"
            / ("plink2.exe" if shutil.which("cmd") else "plink2")
        )
        if candidate.exists():
            plink2 = str(candidate)

    if not plink2:
        pytest.skip("plink2 not available (needed for extraction/PRS integration tests)")

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


def _make_minimal_split_bed_dataset(
    studyA_genetic_dir: Path, tmp_path: Path, *, chrom: int = 9
) -> str:
    """
    Build a '$'-style split dataset containing only a single chromosome.

    This keeps split-data code paths fast: missing chromosomes are detected and skipped
    without running PLINK on 22 inputs.
    """
    base = tmp_path / "geno_split_onechr"
    prefix = base / "studyA_chr$"
    src_prefix = studyA_genetic_dir / f"plink_sampled_{chrom}"
    dst_prefix = Path(str(prefix).replace("$", str(chrom)))
    for ext in [".bed", ".bim", ".fam"]:
        _link_or_copy(src_prefix.with_suffix(ext), dst_prefix.with_suffix(ext))
    return str(prefix)


def test_extract_snps_func_empty_list() -> None:
    from genal.extract_prs import extract_snps_func

    empty = pd.Series([], dtype="object")
    assert extract_snps_func(empty, path="unused") == "FAILED"


def test_extract_snps_func_validates_path(tmp_path: Path) -> None:
    from genal.extract_prs import extract_snps_func

    snps = pd.Series(["rs1"])
    with pytest.raises(TypeError, match="does not lead to valid"):
        extract_snps_func(snps, path=str(tmp_path / "missing_prefix"))


@pytest.mark.slow
@pytest.mark.plink
def test_extract_snps_func_handles_split_data(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.extract_prs import extract_snps_func
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    split_path = _make_minimal_split_bed_dataset(studyA_genetic_dir, Path.cwd(), chrom=9)
    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9"))
    snps = ref["SNP"].astype(str).head(25)

    out_prefix = extract_snps_func(snps, name="pytest_split_extract", path=split_path, ram=4000, cpus=1)
    assert out_prefix.endswith("tmp_GENAL/pytest_split_extract_allchr")
    assert Path(out_prefix + ".pgen").is_file()
    assert Path(out_prefix + ".pvar").is_file()
    assert Path(out_prefix + ".psam").is_file()


def test_prs_func_validates_empty_dataframe(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.extract_prs import prs_func

    _set_plink2_path(reset_genal_config, tests_dir)

    empty = pd.DataFrame({"SNP": pd.Series(dtype=str), "EA": pd.Series(dtype=str), "BETA": pd.Series(dtype=float)})
    with pytest.raises(ValueError, match="No SNPs were extracted|can't be computed"):
        prs_func(empty, path=str(studyA_genetic_dir / "plink_sampled_9"), name="empty_prs", cpus=1, ram=2000)


def test_prs_func_unweighted_sets_beta_to_one(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from genal.extract_prs import prs_func
    from genal.tools import create_tmp

    create_tmp()

    # Avoid running PLINK: stub extraction + PLINK execution, but keep file I/O behavior.
    monkeypatch.setattr("genal.extract_prs.setup_genetic_path", lambda path: (str(tmp_path / "geno"), "bed"))
    monkeypatch.setattr("genal.extract_prs.get_plink_path", lambda: "plink2")
    monkeypatch.setattr("genal.extract_prs.extract_snps_func", lambda snps, name, path, ram=0, cpus=0: "tmp_GENAL/fake_extract")

    def _fake_run(cmd, shell, capture_output, text, check):  # noqa: ANN001
        # Create the files prs_func expects.
        out_prefix = Path("tmp_GENAL") / "tprs_prs"
        (out_prefix.parent).mkdir(parents=True, exist_ok=True)
        (out_prefix.with_suffix(".log")).write_text("--score: 3 variants processed\n")
        (out_prefix.with_suffix(".sscore")).write_text(
            "#FID IID SCORE1_SUM SCORE1_AVG\n0 1 1.0 1.0\n"
        )
        return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr("genal.extract_prs.subprocess.run", _fake_run)

    df = pd.DataFrame({"SNP": ["rs1", "rs2", "rs3"], "EA": ["A", "C", "G"], "BETA": [0.2, -0.1, 0.4]})
    out = prs_func(df, weighted=False, path=str(tmp_path / "geno"), name="tprs", cpus=1, ram=2000)

    score_path = Path("tmp_GENAL") / "tprs_to_prs.txt"
    scored = pd.read_csv(score_path, sep="\t")
    assert (scored["BETA"] == 1).all()
    assert {"FID", "IID"}.issubset(out.columns)


def test_prs_func_with_missing_beta_column_raises_clear_error(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.extract_prs import prs_func

    _set_plink2_path(reset_genal_config, tests_dir)

    df = pd.DataFrame({"SNP": ["rs1"], "EA": ["A"]})
    with pytest.raises(ValueError, match="BETA"):
        prs_func(df, path=str(studyA_genetic_dir / "plink_sampled_9"), name="missing_beta", cpus=1, ram=2000)


@pytest.mark.slow
@pytest.mark.plink
def test_prs_computation_with_real_data(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.extract_prs import prs_func
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    split_path = _make_minimal_split_bed_dataset(studyA_genetic_dir, Path.cwd(), chrom=9)
    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9"))
    subset = ref.head(50)

    weights = pd.DataFrame(
        {
            "SNP": subset["SNP"].astype(str),
            "EA": subset["A1"].astype(str),
            "BETA": np.linspace(-0.2, 0.2, len(subset)),
        }
    )

    scores = prs_func(weights, weighted=True, path=split_path, name="pytest_prs_real", cpus=1, ram=4000)
    assert scores.shape[0] > 0
    assert "IID" in scores.columns
    assert any(col.startswith("SCORE") for col in scores.columns)
