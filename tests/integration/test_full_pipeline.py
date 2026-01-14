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
        pytest.skip("plink2 not available for pipeline tests")

    cfg = read_config()
    cfg["paths"]["plink2_path"] = plink2
    write_config(cfg)
    return plink2


def _install_gene_mapping_for_chr10_region(reset_genal_config, tmp_path: Path) -> None:
    from genal.tools import read_config, write_config

    ref_path = tmp_path / "ref"
    ref_path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "CHR": "10",
                "symbol": "TESTGENE",
                "HGNC_id": "HGNC:TEST",
                "name": "testgene",
                "gene_id": "ENSGTEST",
                "NCBI_id": "TEST",
                "UCSC_id": "ucTEST",
                "Vega_id": "OTTTEST",
                "gene_start_37": 100_000_000,
                "gene_end_37": 100_020_000,
                "gene_start_38": 100_000_000,
                "gene_end_38": 100_020_000,
            }
        ]
    ).to_parquet(ref_path / "gene_id_mapping_filtered.parquet", engine="pyarrow", index=False)

    cfg = read_config()
    cfg["paths"]["ref_path"] = str(ref_path)
    write_config(cfg)


def test_gwas_to_mr_pipeline() -> None:
    from genal.Geno import Geno
    from gwas_examples import _filter_reasonable_rows, load_gwas_example1_b37, make_mr_outcome_from_exposure

    sample = load_gwas_example1_b37(nrows=20_000).df
    df = _filter_reasonable_rows(sample).drop_duplicates(subset=["SNP"]).head(200)

    exposure = Geno(df[["SNP", "BETA", "SE", "EA", "NEA", "EAF"]].copy())
    outcome_df = make_mr_outcome_from_exposure(
        exposure.data,
        causal_beta=2.0,
        noise_seed=0,
        noise_sd=0.01,
    )
    outcome = Geno(outcome_df)

    exposure.query_outcome(outcome, proxy=False)
    res = exposure.MR(methods=["IVW"], action=1, heterogeneity=False, nboot=10, cpus=1)
    assert not res.empty
    assert res.loc[0, "method"] == "Inverse-Variance Weighted"


def test_gene_filter_then_mr_pipeline(reset_genal_config, tmp_path: Path) -> None:
    from genal.Geno import Geno
    from gwas_examples import _filter_reasonable_rows, load_gwas_example1_b37, make_mr_outcome_from_exposure

    _install_gene_mapping_for_chr10_region(reset_genal_config, tmp_path)

    sample = load_gwas_example1_b37(nrows=20_000).df
    df = _filter_reasonable_rows(sample).drop_duplicates(subset=["SNP"]).head(2_000)
    g = Geno(df)

    g.filter_by_gene("TESTGENE", window_size=1_000_000, build="37", replace=True)
    assert g.data.shape[0] > 0

    exposure_df = g.data[["SNP", "BETA", "SE", "EA", "NEA", "EAF"]].head(200).copy()
    exposure = Geno(exposure_df)
    outcome = Geno(make_mr_outcome_from_exposure(exposure_df, causal_beta=1.5, noise_seed=1, noise_sd=0.01))
    exposure.query_outcome(outcome, proxy=False)
    res = exposure.MR(methods=["IVW"], action=1, heterogeneity=False, nboot=10, cpus=1)
    assert not res.empty


def test_liftover_then_mr_pipeline(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from genal.Geno import Geno
    from gwas_examples import make_mr_outcome_from_exposure

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
            "CHR": [1, 1, 1, 1],
            "POS": [1000, 1100, 1200, 1300],
            "EA": ["A", "C", "G", "T"],
            "NEA": ["G", "T", "A", "C"],
            "BETA": [0.1, 0.2, -0.1, 0.05],
            "SE": [0.1, 0.1, 0.1, 0.1],
        }
    )
    g = Geno(df)
    lifted = g.lift(
        chain_file=str(chain),
        liftover_path=None,
        extraction_file=True,
        name=str(tmp_path / "lifted_pipe"),
        replace=False,
    )
    assert "SNP" in lifted.columns

    exposure_df = lifted[["SNP", "BETA", "SE", "EA", "NEA"]].copy()
    outcome_df = make_mr_outcome_from_exposure(exposure_df, causal_beta=2.0, noise_seed=0, noise_sd=0.01)
    exposure = Geno(exposure_df)
    outcome = Geno(outcome_df)
    exposure.query_outcome(outcome, proxy=False)
    res = exposure.MR(methods=["IVW"], action=1, heterogeneity=False, nboot=10, cpus=1)
    assert not res.empty


@pytest.mark.slow
@pytest.mark.plink
def test_gwas_to_prs_pipeline(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path
) -> None:
    from genal.Geno import Geno
    from genal.tools import load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    geno_prefix = str(studyA_genetic_dir / "plink_sampled_9")

    ref = load_reference_panel(geno_prefix).head(500)
    df = pd.DataFrame(
        {
            "SNP": ref["SNP"].astype(str),
            "P": np.linspace(1e-10, 0.9, len(ref)),
            "EA": ref["A1"].astype(str),
            "BETA": np.linspace(-0.2, 0.2, len(ref)),
        }
    )
    g = Geno(df)
    clumped = g.clump(reference_panel=geno_prefix, kb=250, r2=0.5, p1=1, p2=1)
    assert clumped is not None

    clumped.extract_snps(path=geno_prefix)
    clumped.prs(name="pipeline_prs", path=geno_prefix, weighted=True, proxy=False)

    out_csv = Path("pipeline_prs.csv")
    assert out_csv.is_file()
    scored = pd.read_csv(out_csv)
    assert "IID" in scored.columns


@pytest.mark.slow
@pytest.mark.plink
def test_gwas_to_association_pipeline(
    reset_genal_config, tests_dir: Path, studyA_genetic_dir: Path, phenotype_csv_path: Path
) -> None:
    from genal.Geno import Geno
    from genal.tools import load_reference_panel

    _set_plink2_path(reset_genal_config, tests_dir)
    geno_prefix = str(studyA_genetic_dir / "plink_sampled_9")
    ref = load_reference_panel(geno_prefix).head(300)

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
    clumped = g.clump(reference_panel=geno_prefix, kb=250, r2=0.5, p1=1, p2=1)
    assert clumped is not None

    clumped.set_phenotype(pd.read_csv(phenotype_csv_path), IID="IID", FID="FID", PHENO="continuous_column")
    clumped.association_test(path=geno_prefix, covar=["age"], standardize=False)
    assert {"BETA", "SE", "P", "FSTAT"}.issubset(clumped.data.columns)


@pytest.mark.slow
def test_mr_with_mrpresso_outlier_correction(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.Geno import Geno
    from gwas_examples import make_mr_outcome_from_exposure

    # Make MR-PRESSO deterministic/fast for this pipeline test.
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
            "SNP": [f"rs{i}" for i in range(20)],
            "BETA": np.linspace(-0.2, 0.2, 20),
            "SE": 0.1,
            "EA": ["A"] * 20,
            "NEA": ["G"] * 20,
            "EAF": 0.2,
        }
    )
    outcome_df = make_mr_outcome_from_exposure(exposure_df, causal_beta=2.0, noise_seed=0, noise_sd=0.01)
    # Inject an extreme outlier to force outlier removal.
    outcome_df.loc[5, "BETA"] = 10.0

    exp = Geno(exposure_df)
    out = Geno(outcome_df)
    exp.query_outcome(out, proxy=False)

    before = exp.MR(methods=["IVW"], action=1, heterogeneity=False, nboot=10, cpus=1)
    exp.MRpresso(action=1, n_iterations=300, outlier_test=True, distortion_test=False, cpus=1)
    after = exp.MR(methods=["IVW"], action=1, heterogeneity=False, nboot=10, cpus=1, use_mrpresso_data=True)

    assert not before.empty and not after.empty
    assert after.loc[0, "nSNP"] <= before.loc[0, "nSNP"]
