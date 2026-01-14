from __future__ import annotations

from pathlib import Path


def test_import_genal_creates_config_under_test_home(reset_genal_config, tests_dir: Path) -> None:
    import os
    import genal
    from genal.constants import CONFIG_DIR

    config_dir = Path(CONFIG_DIR).resolve()
    home_dir = Path(os.environ["HOME"]).resolve()

    assert config_dir.is_relative_to(home_dir)
    assert config_dir.is_dir()
    assert (config_dir / "config.json").is_file()
    assert genal.__version__


def test_import_is_from_repo_sources() -> None:
    import genal

    # Should resolve to this workspace (not a site-packages install).
    genal_path = Path(genal.__file__).resolve()
    assert "site-packages" not in genal_path.parts

