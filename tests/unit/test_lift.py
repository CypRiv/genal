from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest


def test_prepare_chain_file_with_valid_local_file(tmp_path: Path) -> None:
    from genal.lift import prepare_chain_file

    chain = tmp_path / "local.chain"
    chain.write_text("dummy")

    out = prepare_chain_file(str(chain), start="hg19", end="hg38")
    assert out == str(chain)


def test_prepare_chain_file_raises_on_invalid_path(tmp_path: Path) -> None:
    from genal.lift import prepare_chain_file

    with pytest.raises(ValueError, match="does not lead to a valid file"):
        prepare_chain_file(str(tmp_path / "missing.chain"), start="hg19", end="hg38")


def test_prepare_chain_file_constructs_correct_chain_name_and_downloads_if_missing(
    reset_genal_config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from genal.lift import prepare_chain_file
    from genal.tools import read_config, write_config

    cfg = read_config()
    cfg["paths"]["ref_path"] = str(tmp_path)
    write_config(cfg)

    def _fake_download(url: str, out: str) -> str:
        gz_path = Path(out) / "hg19ToHg38.over.chain.gz"
        gz_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(gz_path, "wb") as f:
            f.write(b"not-a-real-chain-but-enough-for-prepare_chain_file")
        return str(gz_path)

    monkeypatch.setattr("genal.lift.wget.download", _fake_download)

    chain_path = prepare_chain_file(chain_file=None, start="hg19", end="hg38")
    assert chain_path.endswith("hg19ToHg38.over.chain")
    assert Path(chain_path).is_file()


def test_lift_coordinates_python_basic_conversion(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.lift import lift_coordinates_python

    class _StubLiftOver:
        def __init__(self, chain_path: str) -> None:
            self.chain_path = chain_path

        def convert_coordinate(self, chrom: str, pos: int, strand: str):
            chr_num = chrom.replace("chr", "")
            return [(f"chr{chr_num}", int(pos) + 1000, strand)]

    monkeypatch.setattr("genal.lift.LiftOver", _StubLiftOver)

    df = pd.DataFrame({"CHR": [1, 2], "POS": [100, 200]})
    out = lift_coordinates_python(df, chain_path="ignored.chain")

    assert out["CHR"].tolist() == [1, 2]
    assert out["POS"].tolist() == [1100, 1200]


def test_lift_coordinates_python_handles_unlifted_snps(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.lift import lift_coordinates_python

    class _StubLiftOver:
        def __init__(self, chain_path: str) -> None:
            self.chain_path = chain_path

        def convert_coordinate(self, chrom: str, pos: int, strand: str):
            if pos == 200:
                return []
            chr_num = chrom.replace("chr", "")
            return [(f"chr{chr_num}", int(pos) + 1, strand)]

    monkeypatch.setattr("genal.lift.LiftOver", _StubLiftOver)

    df = pd.DataFrame({"CHR": [1, 1, 1], "POS": [100, 200, 300]})
    out = lift_coordinates_python(df, chain_path="ignored.chain")

    assert out["POS"].tolist() == [101, 301]


def test_lift_coordinates_python_selects_first_mapping_when_multiple(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.lift import lift_coordinates_python

    class _StubLiftOver:
        def __init__(self, chain_path: str) -> None:
            self.chain_path = chain_path

        def convert_coordinate(self, chrom: str, pos: int, strand: str):
            chr_num = chrom.replace("chr", "")
            return [
                (f"chr{chr_num}", int(pos) + 10, strand),
                (f"chr{chr_num}", int(pos) + 20, strand),
            ]

    monkeypatch.setattr("genal.lift.LiftOver", _StubLiftOver)

    df = pd.DataFrame({"CHR": [1], "POS": [100]})
    out = lift_coordinates_python(df, chain_path="ignored.chain")
    assert out.loc[0, "POS"] == 110


def test_lift_data_removes_absurd_positions_and_nan(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from genal.lift import lift_data

    class _StubLiftOver:
        def __init__(self, chain_path: str) -> None:
            self.chain_path = chain_path

        def convert_coordinate(self, chrom: str, pos: int, strand: str):
            chr_num = chrom.replace("chr", "")
            return [(f"chr{chr_num}", int(pos) + 1, strand)]

    monkeypatch.setattr("genal.lift.LiftOver", _StubLiftOver)

    chain = tmp_path / "dummy.chain"
    chain.write_text("dummy")

    df = pd.DataFrame(
        {
            "CHR": [1, None, 1, 1],
            "POS": [100, 100, None, 400_000_000],
            "extra": [1, 2, 3, 4],
        }
    )
    out = lift_data(df.copy(), chain_file=str(chain), liftover_path=None)
    assert out.shape[0] == 1
    assert out.loc[0, "POS"] == 101


def test_post_lift_operations_generates_snp_column_and_saves_files(tmp_path: Path) -> None:
    from genal.lift import post_lift_operations

    df = pd.DataFrame({"CHR": [1, 2], "POS": [10, 20]})
    name = str(tmp_path / "lift_out")
    out = post_lift_operations(df, name=name, extraction_file=True)

    assert "SNP" in out.columns
    assert (tmp_path / "lift_out.txt").is_file()
    assert (tmp_path / "lift_out_lifted_extraction.txt").is_file()


@pytest.mark.slow
def test_lift_data_with_real_gwas_example1_b37(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Exercises lift_data end-to-end on a real GWAS sample (without requiring real chain files).
    """
    from gwas_examples import load_gwas_example1_b37
    from genal.lift import lift_data

    class _StubLiftOver:
        def __init__(self, chain_path: str) -> None:
            self.chain_path = chain_path

        def convert_coordinate(self, chrom: str, pos: int, strand: str):
            chr_num = chrom.replace("chr", "")
            return [(f"chr{chr_num}", int(pos) + 123, strand)]

    monkeypatch.setattr("genal.lift.LiftOver", _StubLiftOver)

    chain = tmp_path / "dummy.chain"
    chain.write_text("dummy")

    sample = load_gwas_example1_b37(nrows=20_000).df
    out = lift_data(sample.copy(), chain_file=str(chain), liftover_path=None)

    assert out.shape[0] > 0
    assert out["CHR"].notna().all()
    assert out["POS"].notna().all()


def test_lift_coordinates_python_raises_on_invalid_chain_file_format(tmp_path: Path) -> None:
    from genal.lift import lift_coordinates_python

    bad_chain = tmp_path / "bad.chain"
    # Malformed chain header: should fail during parsing.
    bad_chain.write_text(
        "chain score chr1 X + 0 1000 chr1 1000 + 0 1000 1\n"
        "1000\n"
    )

    df = pd.DataFrame({"CHR": [1], "POS": [100]})

    with pytest.raises(ValueError):
        lift_coordinates_python(df, chain_path=str(bad_chain))
