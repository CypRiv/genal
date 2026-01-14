from __future__ import annotations

import json
import os
import sys
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent

# Ensure local sources are importable even when running from another CWD.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

# Prevent bytecode writes anywhere (repo/site-packages).
sys.dont_write_bytecode = True
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

# If plink2 is vendored under tests/, make it available to tests that only check PATH.
_vendored_plink2_dir = TESTS_DIR / ".genal_test_home" / ".genal" / "plink2"
_vendored_plink2 = _vendored_plink2_dir / "plink2"
if not shutil.which("plink2") and _vendored_plink2.exists():
    os.environ["PATH"] = str(_vendored_plink2_dir) + os.pathsep + os.environ.get("PATH", "")

# Force genal's "~/.genal" config and any matplotlib caches to stay under tests/.
# When running under xdist, isolate per worker to avoid cross-test config races.
worker_id = os.environ.get("PYTEST_XDIST_WORKER", "local")
TEST_HOME = TESTS_DIR / ".genal_test_home" / worker_id
TEST_HOME.mkdir(parents=True, exist_ok=True)
os.environ["HOME"] = str(TEST_HOME)


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--run-network",
        action="store_true",
        default=False,
        help="Run tests marked with @pytest.mark.network (real external API calls).",
    )
    parser.addoption(
        "--keep-tmp",
        action="store_true",
        default=False,
        help="Keep temporary test artifacts under tests/.tmp/ and any repo-root tmp_GENAL/ for debugging.",
    )


def _is_xdist_worker() -> bool:
    return bool(os.environ.get("PYTEST_XDIST_WORKER"))


def pytest_sessionstart(session: pytest.Session) -> None:
    """
    Ensure temporary folders from prior runs don't accumulate.

    Only the controller process performs cleanup to avoid xdist races.
    """
    if _is_xdist_worker():
        return
    if session.config.getoption("--keep-tmp"):
        return

    tmp_root = TESTS_DIR / ".tmp"
    if tmp_root.exists():
        shutil.rmtree(tmp_root, ignore_errors=True)
    tmp_root.mkdir(parents=True, exist_ok=True)

    repo_tmp = REPO_ROOT / "tmp_GENAL"
    if repo_tmp.exists():
        shutil.rmtree(repo_tmp, ignore_errors=True)


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:  # noqa: ARG001
    """Clean up temporary folders created during the run (unless --keep-tmp)."""
    if _is_xdist_worker():
        return
    if session.config.getoption("--keep-tmp"):
        return

    shutil.rmtree(TESTS_DIR / ".tmp", ignore_errors=True)
    shutil.rmtree(REPO_ROOT / "tmp_GENAL", ignore_errors=True)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if not config.getoption("--run-network"):
        skip_network = pytest.mark.skip(reason="network tests skipped by default (use --run-network)")
        for item in items:
            if "network" in item.keywords:
                item.add_marker(skip_network)

    # Markers for PLINK2-backed tests: skip when plink2 isn't available anywhere.
    plink2 = shutil.which("plink2")
    if not plink2:
        candidate = TESTS_DIR / ".genal_test_home" / ".genal" / "plink2" / "plink2"
        if candidate.exists():
            plink2 = str(candidate)
    if not plink2:
        skip_plink = pytest.mark.skip(reason="plink tests skipped (plink2 not available)")
        for item in items:
            if "plink" in item.keywords:
                item.add_marker(skip_plink)


@pytest.fixture(autouse=True)
def _isolate_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep any relative-path writes (tmp_GENAL, etc.) inside tests/."""
    monkeypatch.chdir(tmp_path)


@pytest.fixture(scope="session")
def tests_dir() -> Path:
    return TESTS_DIR


@pytest.fixture(scope="session")
def tests_data_dir(tests_dir: Path) -> Path:
    return tests_dir / "tests_data"


@pytest.fixture(scope="session")
def studyA_genetic_dir(tests_data_dir: Path) -> Path:
    return tests_data_dir / "StudyA_genetic_files"


@pytest.fixture(scope="session")
def phenotype_csv_path(tests_data_dir: Path) -> Path:
    return tests_data_dir / "StudyA_phenotype_example.csv"


@pytest.fixture(scope="session")
def gwas_example1_b37() -> pd.DataFrame:
    """Load GWAS example 1 (GRCh37) for testing."""
    from gwas_examples import load_gwas_example1_b37, _filter_reasonable_rows

    sample = load_gwas_example1_b37(nrows=20_000).df
    return _filter_reasonable_rows(sample)


@pytest.fixture(scope="session")
def gwas_example3_b38() -> pd.DataFrame:
    """Load GWAS example 3 (GRCh38, Parquet) for testing."""
    from gwas_examples import load_gwas_example3_b38_parquet, _filter_reasonable_rows

    sample = load_gwas_example3_b38_parquet(nrows=20_000).df
    return _filter_reasonable_rows(sample)


@pytest.fixture
def synthetic_mr_data() -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """Generate a synthetic exposure/outcome pair for MR tests."""
    from gwas_examples import make_mr_outcome_from_exposure

    exposure = pd.DataFrame(
        {
            "SNP": [f"rs{i}" for i in range(50)],
            "BETA": np.linspace(-0.2, 0.2, 50),
            "SE": 0.1,
            "EA": ["A"] * 50,
            "NEA": ["G"] * 50,
            "EAF": 0.2,
        }
    )
    causal_beta = 2.0
    outcome = make_mr_outcome_from_exposure(exposure, causal_beta=causal_beta, noise_seed=0, noise_sd=0.01)
    return exposure, outcome, causal_beta


@pytest.fixture()
def reset_genal_config() -> Path:
    """
    Reset genal's config.json to defaults inside the test HOME.

    Many genal helpers read/write ~/.genal/config.json; this fixture keeps tests isolated
    and deterministic.
    """
    from genal.tools import default_config
    from genal.constants import CONFIG_DIR

    config_dir = Path(CONFIG_DIR)
    config_dir.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "config.json"
    config_path.write_text(json.dumps(default_config(), indent=4))
    return config_path
