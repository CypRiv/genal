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
        pytest.skip("plink2 not available (needed for association integration tests)")

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
    base = tmp_path / "geno_split_onechr"
    prefix = base / "studyA_chr$"
    src_prefix = studyA_genetic_dir / f"plink_sampled_{chrom}"
    dst_prefix = Path(str(prefix).replace("$", str(chrom)))
    for ext in [".bed", ".bim", ".fam"]:
        _link_or_copy(src_prefix.with_suffix(ext), dst_prefix.with_suffix(ext))
    return str(prefix)


def test_prepare_psam_file_merges_phenotype(tmp_path: Path) -> None:
    from genal.association import _prepare_psam_file

    genetic_path = str(tmp_path / "geno")
    psam = pd.DataFrame({"#FID": [0, 0], "IID": [1, 2], "SEX": ["", None]})
    psam.to_csv(genetic_path + ".psam", sep="\t", index=False)

    pheno = pd.DataFrame({"FID": [0, 0], "IID": [1, 2], "PHENO": [0, 1]})
    out = _prepare_psam_file(genetic_path, pheno, pheno_type="binary", standardize=False)

    assert "PHENO1" in out.columns
    assert set(out["PHENO1"]) == {"1", "2"}
    assert set(out["SEX"]) == {"NA"}


def test_prepare_psam_validates_id_mismatch(tmp_path: Path) -> None:
    from genal.association import _prepare_psam_file

    genetic_path = str(tmp_path / "geno")
    psam = pd.DataFrame({"#FID": [0, 0], "IID": [1, 2], "SEX": [1, 2]})
    psam.to_csv(genetic_path + ".psam", sep="\t", index=False)

    pheno = pd.DataFrame({"FID": [0], "IID": [999], "PHENO": [0]})
    with pytest.raises(ValueError, match="inconsistent"):
        _prepare_psam_file(genetic_path, pheno, pheno_type="binary", standardize=False)


def test_handle_covariates_removes_single_value_and_non_numeric(tmp_path: Path) -> None:
    from genal.association import _handle_covariates
    from genal.tools import create_tmp

    create_tmp()
    df = pd.DataFrame(
        {
            "FID": [0, 0, 0],
            "IID": [1, 2, 3],
            "age": [10, 11, 12],
            "sex": ["male", "female", "male"],
            "const": [1, 1, 1],
        }
    )

    covars, covar_path = _handle_covariates(["age", "sex", "const"], df, name="pytest_covar")
    assert covars == ["age"]
    assert covar_path is not None
    assert Path(covar_path).is_file()


def test_process_results_plink2_allele_alignment(tmp_path: Path) -> None:
    from genal.association import _process_results_plink2

    out_prefix = str(tmp_path / "assoc")
    results_path = out_prefix + ".PHENO1.glm.linear"

    assoc = pd.DataFrame(
        {
            "#CHROM": [1, 1, 1],
            "POS": [10, 20, 30],
            "A1": ["A", "G", "C"],
            "TEST": ["ADD", "ADD", "ADD"],
            "BETA": [0.5, 0.5, 0.5],
            "SE": [0.1, 0.1, 0.1],
            "P": [0.1, 0.1, 0.1],
        }
    )
    assoc.to_csv(results_path, sep="\t", index=False)

    data = pd.DataFrame(
        {
            "CHR": [1, 1, 1],
            "POS": [10, 20, 30],
            "EA": ["A", "A", "A"],
            "NEA": ["G", "G", "G"],
        }
    )
    out = _process_results_plink2(out_prefix, data, pheno_type="quant")
    assert out["POS"].tolist() == [10, 20]
    beta = dict(zip(out["POS"], out["BETA"]))
    assert beta[10] == pytest.approx(0.5)
    assert beta[20] == pytest.approx(-0.5)


@pytest.mark.slow
@pytest.mark.plink
def test_association_test_quantitative_phenotype(
    reset_genal_config,
    tests_dir: Path,
    studyA_genetic_dir: Path,
    phenotype_csv_path: Path,
) -> None:
    from genal.association import association_test_func_plink2, set_phenotype_func
    from genal.extract_prs import extract_snps_func
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    split_path = _make_minimal_split_bed_dataset(studyA_genetic_dir, Path.cwd(), chrom=9)
    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9"))
    subset = ref.head(50)
    data = pd.DataFrame(
        {
            "CHR": subset["CHR"].astype(int),
            "POS": subset["POS"].astype(int),
            "SNP": subset["SNP"].astype(str),
            "EA": subset["A1"].astype(str),
            "NEA": subset["A2"].astype(str),
            "BETA": np.linspace(-0.1, 0.1, len(subset)),
            "SE": 0.1,
            "P": 0.5,
        }
    )

    name = "pytest_assoc_quant"
    _ = extract_snps_func(data["SNP"], name=name, path=split_path, ram=4000, cpus=1)

    pheno_raw = pd.read_csv(phenotype_csv_path)
    pheno, pheno_type = set_phenotype_func(pheno_raw, PHENO="continuous_column", PHENO_type=None, IID="IID", FID="FID")
    assert pheno_type == "quant"

    res = association_test_func_plink2(
        data,
        covar_list=[],
        standardize=False,
        name=name,
        data_pheno=pheno,
        pheno_type=pheno_type,
    )
    assert res.shape[0] > 0
    assert {"CHR", "POS", "BETA", "SE", "P"}.issubset(res.columns)


@pytest.mark.slow
@pytest.mark.plink
def test_association_test_binary_phenotype(
    reset_genal_config,
    tests_dir: Path,
    studyA_genetic_dir: Path,
    phenotype_csv_path: Path,
) -> None:
    from genal.association import association_test_func_plink2, set_phenotype_func
    from genal.extract_prs import extract_snps_func
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    split_path = _make_minimal_split_bed_dataset(studyA_genetic_dir, Path.cwd(), chrom=9)
    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9"))
    subset = ref.head(50)
    data = pd.DataFrame(
        {
            "CHR": subset["CHR"].astype(int),
            "POS": subset["POS"].astype(int),
            "SNP": subset["SNP"].astype(str),
            "EA": subset["A1"].astype(str),
            "NEA": subset["A2"].astype(str),
            "BETA": np.linspace(-0.1, 0.1, len(subset)),
            "SE": 0.1,
            "P": 0.5,
        }
    )

    name = "pytest_assoc_bin"
    _ = extract_snps_func(data["SNP"], name=name, path=split_path, ram=4000, cpus=1)

    pheno_raw = pd.read_csv(phenotype_csv_path)
    pheno, pheno_type = set_phenotype_func(pheno_raw, PHENO="binary_column", PHENO_type=None, IID="IID", FID="FID")
    assert pheno_type == "binary"

    res = association_test_func_plink2(
        data,
        covar_list=[],
        standardize=False,
        name=name,
        data_pheno=pheno,
        pheno_type=pheno_type,
    )
    assert res.shape[0] > 0
    assert {"CHR", "POS", "BETA", "SE", "P"}.issubset(res.columns)


@pytest.mark.slow
@pytest.mark.plink
def test_association_test_with_covariates(
    reset_genal_config,
    tests_dir: Path,
    studyA_genetic_dir: Path,
    phenotype_csv_path: Path,
) -> None:
    from genal.association import association_test_func_plink2, set_phenotype_func
    from genal.extract_prs import extract_snps_func
    from genal.tools import create_tmp, load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    create_tmp()

    split_path = _make_minimal_split_bed_dataset(studyA_genetic_dir, Path.cwd(), chrom=9)
    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9"))
    subset = ref.head(50)
    data = pd.DataFrame(
        {
            "CHR": subset["CHR"].astype(int),
            "POS": subset["POS"].astype(int),
            "SNP": subset["SNP"].astype(str),
            "EA": subset["A1"].astype(str),
            "NEA": subset["A2"].astype(str),
            "BETA": np.linspace(-0.1, 0.1, len(subset)),
            "SE": 0.1,
            "P": 0.5,
        }
    )

    name = "pytest_assoc_covar"
    _ = extract_snps_func(data["SNP"], name=name, path=split_path, ram=4000, cpus=1)

    pheno_raw = pd.read_csv(phenotype_csv_path)
    pheno, pheno_type = set_phenotype_func(pheno_raw, PHENO="continuous_column", PHENO_type=None, IID="IID", FID="FID")
    assert pheno_type == "quant"

    # Include one numeric covariate and one non-numeric; association helper should drop the latter.
    res = association_test_func_plink2(
        data,
        covar_list=["age", "sex"],
        standardize=False,
        name=name,
        data_pheno=pheno,
        pheno_type=pheno_type,
    )
    assert res.shape[0] > 0
