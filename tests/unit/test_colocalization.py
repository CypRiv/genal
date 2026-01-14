from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_combine_abf_probabilities_sum_to_one() -> None:
    from genal.colocalization import combine_abf

    l1 = np.array([0.1, 0.2, 0.3])
    l2 = np.array([0.2, 0.1, 0.0])
    pp = combine_abf(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5)
    s = sum(pp.values())
    assert s == pytest.approx(1.0, rel=1e-12)
    assert all(0 <= v <= 1 for v in pp.values())


def test_coloc_abf_func_runs_on_small_overlap() -> None:
    from genal.colocalization import coloc_abf_func

    d1 = pd.DataFrame(
        {
            "CHR": [1, 1, 1],
            "POS": [100, 200, 300],
            "SNP": ["rs1", "rs2", "rs3"],
            "EA": ["A", "A", "C"],
            "NEA": ["C", "G", "T"],
            "BETA": [0.1, 0.2, -0.1],
            "SE": [0.05, 0.05, 0.05],
            "EAF": [0.1, 0.2, 0.3],
        }
    )
    d2 = pd.DataFrame(
        {
            "CHR": [1, 1, 1],
            "POS": [100, 200, 300],
            "SNP": ["rs1", "rs2", "rs3"],
            "EA": ["A", "G", "C"],  # rs2 inverted => BETA2 should flip
            "NEA": ["C", "A", "T"],
            "BETA": [0.1, -0.2, -0.1],
            "SE": [0.05, 0.05, 0.05],
            "EAF": [0.1, 0.2, 0.3],
        }
    )

    out = coloc_abf_func(d1.copy(), d2.copy(), sdY1=1, sdY2=1, merge_on_snp=False)
    assert out["nsnps"] == 3
    assert all(k in out for k in ["PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"])
    assert sum(out[k] for k in ["PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"]) == pytest.approx(1.0, rel=1e-12)


def test_coloc_abf_func_raises_on_no_overlap() -> None:
    from genal.colocalization import coloc_abf_func

    d1 = pd.DataFrame({"CHR": [1], "POS": [100], "BETA": [0.1], "SE": [0.1]})
    d2 = pd.DataFrame({"CHR": [1], "POS": [200], "BETA": [0.1], "SE": [0.1]})

    with pytest.raises(ValueError, match="No overlapping variants"):
        coloc_abf_func(d1, d2, sdY1=1, sdY2=1)

