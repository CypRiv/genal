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
        pytest.skip("plink2 not available for Geno method tests")

    cfg = read_config()
    cfg["paths"]["plink2_path"] = plink2
    write_config(cfg)
    return plink2


def _write_gene_info(ref_path: Path) -> None:
    ref_path.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        [
            {
                "CHR": "1",
                "symbol": "GENE1",
                "HGNC_id": "HGNC:1",
                "name": "gene1",
                "gene_id": "ENSG00000000001",
                "NCBI_id": "1",
                "UCSC_id": "uc000000.1",
                "Vega_id": "OTT000000001",
                "gene_start_37": 1000,
                "gene_end_37": 2000,
                "gene_start_38": 1000,
                "gene_end_38": 2000,
            }
        ]
    )
    df.to_parquet(ref_path / "gene_id_mapping_filtered.parquet", engine="pyarrow", index=False)


def test_geno_lift_method(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from genal.Geno import Geno

    class _StubLiftOver:
        def __init__(self, chain_path: str) -> None:  # noqa: ARG002
            return

        def convert_coordinate(self, chrom: str, pos: int, strand: str):  # noqa: ARG002
            chr_num = chrom.replace("chr", "")
            return [(f"chr{chr_num}", int(pos) + 100, "+")]

    monkeypatch.setattr("genal.lift.LiftOver", _StubLiftOver)

    chain = tmp_path / "dummy.chain"
    chain.write_text("dummy")

    df = pd.DataFrame(
        {
            "CHR": [1, 1],
            "POS": [1000, 2000],
            "SNP": ["rs1", "rs2"],
            "EA": ["A", "C"],
            "NEA": ["G", "T"],
            "BETA": [0.1, 0.2],
            "SE": [0.1, 0.1],
            "P": [0.5, 0.5],
        }
    )
    g = Geno(df)
    lifted = g.lift(chain_file=str(chain), liftover_path=None, replace=False)
    assert lifted["POS"].tolist() == [1100, 2100]
    assert g.data["POS"].tolist() == [1000, 2000]  # replace=False keeps original


def test_geno_filter_by_gene_method(reset_genal_config, tmp_path: Path) -> None:
    from genal.Geno import Geno
    from genal.tools import read_config, write_config

    ref_path = tmp_path / "ref"
    _write_gene_info(ref_path)

    cfg = read_config()
    cfg["paths"]["ref_path"] = str(ref_path)
    write_config(cfg)

    df = pd.DataFrame(
        {
            "CHR": [1, 1, 2],
            "POS": [900, 1500, 1500],
            "SNP": ["rs1", "rs2", "rs3"],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "A"],
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0.1, 0.1, 0.1],
            "P": [0.5, 0.5, 0.5],
        }
    )
    g = Geno(df)
    out = g.filter_by_gene("GENE1", window_size=1000, build="37", replace=False)
    assert out.data["CHR"].unique().tolist() == [1]
    assert "Distance" in out.data.columns


@pytest.mark.slow
@pytest.mark.plink
def test_geno_clump_method(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.Geno import Geno
    from genal.tools import load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    ref_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(ref_prefix).head(200)

    df = pd.DataFrame(
        {
            "SNP": ref["SNP"].astype(str),
            "P": np.linspace(1e-10, 0.9, len(ref)),
            "EA": ref["A1"].astype(str),
            "NEA": ref["A2"].astype(str),
            "BETA": np.linspace(-0.2, 0.2, len(ref)),
        }
    )
    g = Geno(df)
    clumped = g.clump(reference_panel=ref_prefix, kb=250, r2=0.5, p1=1, p2=1)
    assert clumped is not None
    assert clumped.data.shape[0] >= 1


@pytest.mark.slow
@pytest.mark.plink
def test_geno_extract_snps_method(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from genal.Geno import Geno
    import importlib

    _set_plink2_path(reset_genal_config, tests_dir)

    # Make the generated tmp name deterministic.
    geno_mod = importlib.import_module("genal.Geno")
    monkeypatch.setattr(geno_mod.uuid, "uuid4", lambda: "deadbeefdeadbeef")

    df = pd.DataFrame({"SNP": ["rs141734683", "rs56377469"], "P": [1e-6, 1e-4]})
    g = Geno(df)
    g.extract_snps(path=str(studyA_genetic_dir / "plink_sampled_9"))

    out_prefix = Path("tmp_GENAL") / "deadbeef_allchr"
    assert out_prefix.with_suffix(".pgen").is_file()
    assert out_prefix.with_suffix(".pvar").is_file()
    assert out_prefix.with_suffix(".psam").is_file()


@pytest.mark.slow
@pytest.mark.plink
def test_geno_prs_method(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.Geno import Geno
    from genal.tools import load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9")).head(50)
    df = pd.DataFrame(
        {
            "SNP": ref["SNP"].astype(str),
            "EA": ref["A1"].astype(str),
            "BETA": np.linspace(-0.1, 0.1, len(ref)),
        }
    )
    g = Geno(df)
    g.prs(name="geno_prs_out", path=str(studyA_genetic_dir / "plink_sampled_9"), weighted=True, proxy=False)

    out_csv = Path("geno_prs_out.csv")
    assert out_csv.is_file()
    scored = pd.read_csv(out_csv)
    assert "IID" in scored.columns


@pytest.mark.slow
@pytest.mark.plink
def test_geno_association_test_method(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path, phenotype_csv_path: Path
) -> None:
    from genal.Geno import Geno
    from genal.tools import load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)

    ref = load_reference_panel(str(studyA_genetic_dir / "plink_sampled_9")).head(50)
    df = pd.DataFrame(
        {
            "CHR": ref["CHR"].astype(int),
            "POS": ref["POS"].astype(int),
            "SNP": ref["SNP"].astype(str),
            "EA": ref["A1"].astype(str),
            "NEA": ref["A2"].astype(str),
            "P": np.linspace(1e-10, 0.9, len(ref)),
        }
    )
    g = Geno(df)
    g.set_phenotype(pd.read_csv(phenotype_csv_path), IID="IID", FID="FID", PHENO="continuous_column")
    g.association_test(path=str(studyA_genetic_dir / "plink_sampled_9"), covar=["age"], standardize=False)

    assert {"BETA", "SE", "P", "FSTAT"}.issubset(g.data.columns)
    assert g.data.shape[0] > 0


@pytest.mark.network
def test_geno_query_gwas_catalog_method() -> None:
    from genal.Geno import Geno

    g = Geno(pd.DataFrame({"SNP": ["rs7412"]}))
    out = g.query_gwas_catalog(p_threshold=1, max_associations=1, timeout=30)
    assert "ASSOC" in out.columns


def test_geno_mrpresso_method(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.Geno import Geno
    from gwas_examples import make_mr_outcome_from_exposure

    # Patch MR-PRESSO to run inline (no process spawning).
    class _InlineProcessPool:
        def __init__(self, max_workers: int | None = None) -> None:  # noqa: ARG002
            return

        def __enter__(self):  # noqa: ANN001
            return self

        def __exit__(self, exc_type, exc, tb) -> bool:  # noqa: ANN001
            return False

        def map(self, func, iterable):  # noqa: ANN001
            return list(map(func, iterable))

    monkeypatch.setattr("genal.MRpresso.ProcessPoolExecutor", _InlineProcessPool)
    monkeypatch.setattr("genal.MRpresso.tqdm", lambda it, **kwargs: it)

    seed_base = 123
    counter = {"n": 0}

    def _seeded_default_rng():  # noqa: ANN001
        counter["n"] += 1
        return np.random.default_rng(seed_base + counter["n"])

    monkeypatch.setattr("genal.MRpresso.default_rng", _seeded_default_rng)

    exposure_df = pd.DataFrame(
        {
            "SNP": [f"rs{i}" for i in range(10)],
            "BETA": np.linspace(-0.2, 0.2, 10),
            "SE": 0.1,
            "EA": ["A"] * 10,
            "NEA": ["G"] * 10,
            "EAF": 0.2,
        }
    )
    outcome_df = make_mr_outcome_from_exposure(exposure_df, causal_beta=2.0, noise_seed=0, noise_sd=0.01)

    exp = Geno(exposure_df)
    out = Geno(outcome_df)
    exp.query_outcome(out, proxy=False)

    mod_table, global_test, outlier_test, bias_test = exp.MRpresso(
        action=1,
        n_iterations=200,
        outlier_test=False,
        distortion_test=False,
        cpus=1,
    )
    assert not mod_table.empty
    assert "global_test_p" in global_test
    assert outlier_test.empty
    assert bias_test == {}
