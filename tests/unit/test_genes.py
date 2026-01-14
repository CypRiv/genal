from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


def _write_gene_info(ref_path: Path, rows: list[dict]) -> Path:
    ref_path.mkdir(parents=True, exist_ok=True)
    out = ref_path / "gene_id_mapping_filtered.parquet"
    pd.DataFrame(rows).to_parquet(out, engine="pyarrow", index=False)
    return out


def _set_ref_path(reset_genal_config, ref_path: Path) -> None:
    from genal.tools import read_config, write_config

    cfg = read_config()
    cfg["paths"]["ref_path"] = str(ref_path)
    write_config(cfg)


def test_filter_by_gene_func_validates_id_type() -> None:
    from genal.genes import filter_by_gene_func

    df = pd.DataFrame({"CHR": [1], "POS": [100]})
    with pytest.raises(ValueError, match="Invalid id_type"):
        filter_by_gene_func(df, gene_identifier="APOE", id_type="BAD")


def test_filter_by_gene_func_validates_build() -> None:
    from genal.genes import filter_by_gene_func

    df = pd.DataFrame({"CHR": [1], "POS": [100]})
    with pytest.raises(ValueError, match="Invalid build"):
        filter_by_gene_func(df, gene_identifier="APOE", build="39")


def test_filter_by_gene_func_converts_id_types_and_finds_gene(
    reset_genal_config, tmp_path: Path
) -> None:
    from genal.genes import filter_by_gene_func

    ref_path = tmp_path / "ref"
    _set_ref_path(reset_genal_config, ref_path)

    _write_gene_info(
        ref_path,
        [
            {
                "CHR": "19",
                "symbol": "APOE",
                "HGNC_id": "HGNC:613",
                "name": "apolipoprotein E",
                "gene_id": "ENSG00000130203",
                "NCBI_id": "348",
                "UCSC_id": "uc002mbe.4",
                "Vega_id": "OTTHUMG00000019505",
                "gene_start_37": 1000,
                "gene_end_37": 2000,
                "gene_start_38": 1100,
                "gene_end_38": 2100,
            }
        ],
    )

    data = pd.DataFrame({"CHR": [19, 19], "POS": [1500, 10_000]})

    out_hgnc = filter_by_gene_func(data, "HGNC:613", id_type="HGNC", window_size=1000, build="37")
    assert out_hgnc["POS"].tolist() == [1500]

    out_ensembl = filter_by_gene_func(
        data, "ENSG00000130203", id_type="Ensembl", window_size=1000, build="37"
    )
    assert out_ensembl["POS"].tolist() == [1500]


def test_filter_by_gene_func_gene_not_found_raises(reset_genal_config, tmp_path: Path) -> None:
    from genal.genes import filter_by_gene_func

    ref_path = tmp_path / "ref"
    _set_ref_path(reset_genal_config, ref_path)
    _write_gene_info(
        ref_path,
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
                "gene_start_37": 100,
                "gene_end_37": 200,
                "gene_start_38": 100,
                "gene_end_38": 200,
            }
        ],
    )

    df = pd.DataFrame({"CHR": [1], "POS": [150]})
    with pytest.raises(ValueError, match="not found"):
        filter_by_gene_func(df, gene_identifier="MISSING", id_type="symbol", build="37")


def test_filter_by_gene_func_window_and_distance_calculation(reset_genal_config, tmp_path: Path) -> None:
    from genal.genes import filter_by_gene_func

    ref_path = tmp_path / "ref"
    _set_ref_path(reset_genal_config, ref_path)
    _write_gene_info(
        ref_path,
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
        ],
    )

    # window_size=1000 -> +/-500 around gene, i.e. [500, 2500]
    df = pd.DataFrame({"CHR": [1] * 6, "POS": [499, 500, 900, 1500, 2100, 2501]})
    out = filter_by_gene_func(df, gene_identifier="GENE1", window_size=1000, build="37")

    assert out["POS"].tolist() == [500, 900, 1500, 2100]
    # Distances: before=-100, inside=0, after=+100
    dist = dict(zip(out["POS"], out["Distance"]))
    assert dist[900] == -100
    assert dist[1500] == 0
    assert dist[2100] == 100


def test_filter_by_gene_func_handles_chromosome_x(reset_genal_config, tmp_path: Path) -> None:
    from genal.genes import filter_by_gene_func

    ref_path = tmp_path / "ref"
    _set_ref_path(reset_genal_config, ref_path)
    _write_gene_info(
        ref_path,
        [
            {
                "CHR": "X",
                "symbol": "GENEX",
                "HGNC_id": "HGNC:X",
                "name": "genex",
                "gene_id": "ENSGX",
                "NCBI_id": "X",
                "UCSC_id": "ucX",
                "Vega_id": "OTTX",
                "gene_start_37": 100,
                "gene_end_37": 200,
                "gene_start_38": 100,
                "gene_end_38": 200,
            }
        ],
    )

    df = pd.DataFrame({"CHR": [23, 23, 22], "POS": [150, 9999, 150]})
    out = filter_by_gene_func(df, gene_identifier="GENEX", window_size=1000, build="37")
    assert out["CHR"].unique().tolist() == [23]
    assert out["POS"].tolist() == [150]


def test_filter_by_gene_func_multiple_entries_uses_first(reset_genal_config, tmp_path: Path) -> None:
    from genal.genes import filter_by_gene_func

    ref_path = tmp_path / "ref"
    _set_ref_path(reset_genal_config, ref_path)
    _write_gene_info(
        ref_path,
        [
            {
                "CHR": "1",
                "symbol": "DUPGENE",
                "HGNC_id": "HGNC:DUP",
                "name": "dup",
                "gene_id": "ENSGDUP",
                "NCBI_id": "DUP",
                "UCSC_id": "ucDUP",
                "Vega_id": "OTTDUP",
                "gene_start_37": 1000,
                "gene_end_37": 1100,
                "gene_start_38": 1000,
                "gene_end_38": 1100,
            },
            {
                "CHR": "1",
                "symbol": "DUPGENE",
                "HGNC_id": "HGNC:DUP",
                "name": "dup-alt",
                "gene_id": "ENSGDUP2",
                "NCBI_id": "DUP2",
                "UCSC_id": "ucDUP2",
                "Vega_id": "OTTDUP2",
                "gene_start_37": 5000,
                "gene_end_37": 5100,
                "gene_start_38": 5000,
                "gene_end_38": 5100,
            },
        ],
    )

    df = pd.DataFrame({"CHR": [1, 1], "POS": [1050, 5050]})
    out = filter_by_gene_func(df, gene_identifier="DUPGENE", window_size=500, build="37")
    assert out["POS"].tolist() == [1050]


def test_filter_by_gene_func_gene_on_chromosome_boundary_clamps_window_start(
    reset_genal_config, tmp_path: Path
) -> None:
    from genal.genes import filter_by_gene_func

    ref_path = tmp_path / "ref"
    _set_ref_path(reset_genal_config, ref_path)
    _write_gene_info(
        ref_path,
        [
            {
                "CHR": "1",
                "symbol": "EDGE",
                "HGNC_id": "HGNC:EDGE",
                "name": "edge",
                "gene_id": "ENSGEDGE",
                "NCBI_id": "EDGE",
                "UCSC_id": "ucEDGE",
                "Vega_id": "OTTEDGE",
                "gene_start_37": 100,
                "gene_end_37": 200,
                "gene_start_38": 100,
                "gene_end_38": 200,
            }
        ],
    )

    # window_size=1000 -> gene_start - 500 would be negative, should clamp to 0.
    df = pd.DataFrame({"CHR": [1, 1], "POS": [0, 10_000]})
    out = filter_by_gene_func(df, gene_identifier="EDGE", window_size=1000, build="37")
    assert out["POS"].tolist() == [0]


@pytest.mark.slow
def test_filter_by_gene_with_real_gwas_data(reset_genal_config, tests_data_dir: Path, tmp_path: Path) -> None:
    """
    Slow integration-like unit test: extracts a small region from the large GWAS example and
    ensures gene window filtering works on real-looking coordinates.
    """
    from genal.genes import filter_by_gene_func

    ref_path = tmp_path / "ref"
    _set_ref_path(reset_genal_config, ref_path)

    # APOE coordinates (GRCh37) are ~44.905–44.909 Mb on chr19.
    _write_gene_info(
        ref_path,
        [
            {
                "CHR": "19",
                "symbol": "APOE",
                "HGNC_id": "HGNC:613",
                "name": "apolipoprotein E",
                "gene_id": "ENSG00000130203",
                "NCBI_id": "348",
                "UCSC_id": "uc002mbe.4",
                "Vega_id": "OTTHUMG00000019505",
                "gene_start_37": 44_905_754,
                "gene_end_37": 44_909_393,
                "gene_start_38": 44_905_754,
                "gene_end_38": 44_909_393,
            }
        ],
    )

    gwas_path = tests_data_dir / "GWAS_example1_b37.txt"
    rows: list[dict] = []
    with gwas_path.open("r") as f:
        header = next(f)
        for line in f:
            if not line.startswith("19:4490"):
                continue
            marker = line.split(None, 1)[0]
            chr_str, pos_str, _ = marker.split(":", 2)
            rows.append({"CHR": int(chr_str), "POS": int(pos_str)})
            if len(rows) >= 500:
                break

    df = pd.DataFrame(rows)
    assert df.shape[0] > 0

    out = filter_by_gene_func(df, gene_identifier="APOE", window_size=1_000_000, build="37")
    assert out.shape[0] > 0
    assert (out["CHR"] == 19).all()
    assert "Distance" in out.columns

