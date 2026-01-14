from __future__ import annotations

import argparse
import datetime as _dt
import os
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def _parse_junit_counts(junit_path: Path) -> dict[str, int]:
    if not junit_path.exists():
        return {"tests": 0, "failures": 0, "errors": 0, "skipped": 0}

    tree = ET.parse(junit_path)
    root = tree.getroot()

    # pytest may output <testsuite> or <testsuites>.
    if root.tag == "testsuite":
        suites = [root]
    elif root.tag == "testsuites":
        suites = list(root.findall("testsuite"))
    else:
        suites = list(root.findall(".//testsuite"))

    counts = {"tests": 0, "failures": 0, "errors": 0, "skipped": 0}
    for suite in suites:
        for key in counts:
            try:
                counts[key] += int(suite.attrib.get(key, 0))
            except ValueError:
                pass
    return counts


def _extract_warnings_block(pytest_output: str) -> str | None:
    lines = pytest_output.splitlines()
    start = None
    for idx, line in enumerate(lines):
        if line.startswith("=============================== warnings summary"):
            start = idx
            break
    if start is None:
        return None

    end = None
    for idx in range(start + 1, len(lines)):
        if lines[idx].startswith("=========================== short test summary"):
            end = idx
            break
    if end is None:
        end = len(lines)

    return "\n".join(lines[start:end]).strip() or None


def _extract_short_summary_lines(pytest_output: str) -> list[str]:
    return [
        line
        for line in pytest_output.splitlines()
        if line.startswith(("SKIPPED", "XFAIL", "XPASS", "XFAILED", "XPASSED"))
    ]


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the genal pytest suite and write artifacts under tests/reports/.")
    parser.add_argument(
        "--run-network",
        action="store_true",
        help="Include @pytest.mark.network tests (real external API calls).",
    )
    parser.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep temporary artifacts under tests/.tmp/ (passed through to pytest).",
    )
    parser.add_argument(
        "--pytest-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Additional args to pass to pytest (everything after --pytest-args).",
    )
    return parser.parse_args(argv)


def main() -> int:
    args = _parse_args(sys.argv[1:])

    repo_root = Path(__file__).resolve().parents[1]
    tests_dir = repo_root / "tests"
    reports_dir = tests_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    test_home = tests_dir / ".genal_test_home"
    test_home.mkdir(parents=True, exist_ok=True)

    # Keep all side effects (genal config, matplotlib caches, pyc files) under tests/.
    sys.dont_write_bytecode = True
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["HOME"] = str(test_home)
    os.environ.setdefault("MPLCONFIGDIR", str(test_home / ".matplotlib"))

    env = os.environ.copy()

    junit_path = reports_dir / "junit.xml"
    stdout_path = reports_dir / "pytest_stdout.txt"

    started = _dt.datetime.now(_dt.timezone.utc)
    cmd = [
        sys.executable,
        "-m",
        "pytest",
        "-c",
        str(tests_dir / "pytest.ini"),
        str(tests_dir),
        "--junitxml",
        str(junit_path),
    ]
    if args.run_network:
        cmd.append("--run-network")
    if args.keep_tmp:
        cmd.append("--keep-tmp")
    if args.pytest_args:
        cmd.extend(args.pytest_args)

    proc = subprocess.run(
        cmd,
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
    )
    finished = _dt.datetime.now(_dt.timezone.utc)

    stdout_path.write_text(
        proc.stdout + ("\n\n[stderr]\n" + proc.stderr if proc.stderr else "")
    )

    pytest_output = stdout_path.read_text()
    warnings_block = _extract_warnings_block(pytest_output)
    short_summary = _extract_short_summary_lines(pytest_output)

    counts = _parse_junit_counts(junit_path)
    passed = counts["tests"] - counts["failures"] - counts["errors"] - counts["skipped"]

    report_path = reports_dir / f"run_{started.strftime('%Y%m%d_%H%M%S')}_UTC.md"
    duration_s = (finished - started).total_seconds()

    try:
        if str(repo_root) not in sys.path:
            sys.path.insert(0, str(repo_root))
        import genal
        genal_version = getattr(genal, "__version__", "unknown")
    except Exception:
        genal_version = "unknown (import failed)"

    skips_section = "None."
    if short_summary:
        skips_section = "```\n" + "\n".join(short_summary) + "\n```"

    warnings_section = "None."
    if warnings_block:
        warnings_section = "```\n" + warnings_block + "\n```"

    report_text = "\n".join(
        [
            "# genal-python test run",
            "",
            f"- Started (UTC): {started.isoformat()}",
            f"- Finished (UTC): {finished.isoformat()}",
            f"- Duration: {duration_s:.2f}s",
            f"- Python: {sys.version.split()[0]}",
            f"- Platform: {platform.platform()}",
            f"- genal: {genal_version}",
            "",
            "## Command",
            "",
            "```",
            " ".join(cmd),
            "```",
            "",
            "## Results",
            "",
            f"- Passed: {passed}",
            f"- Skipped: {counts['skipped']}",
            f"- Failed: {counts['failures']}",
            f"- Errors: {counts['errors']}",
            f"- Total: {counts['tests']}",
            f"- Exit code: {proc.returncode}",
            "",
            "## Artifacts",
            "",
            f"- JUnit XML: `{junit_path.relative_to(repo_root)}`",
            f"- Pytest output: `{stdout_path.relative_to(repo_root)}`",
            "",
            "## Skips / Expected Failures",
            "",
            skips_section,
            "",
            "## Warnings",
            "",
            warnings_section,
            "",
            "## QA Notes / Suggested Follow-ups (code changes not made here)",
            "",
            "- `genal/__init__.py` writes to `~/.genal/` on import; consider lazy/explicit initialization to avoid import-time side effects.",
            "- `genal/tools.py` and PLINK wrappers rely on a relative `tmp_GENAL/` folder; consider using `tempfile` or a configurable working directory.",
            "- `genal/clump.py:clump_data_plink2()` writes into `tmp_GENAL/` without creating it; calling `create_tmp()` internally would make it safer when used outside `Geno`.",
            "- `genal/Geno.py:Geno.MR()` raises `KeyError` when `query_outcome()` yields <2 instruments because it slices columns on an empty `res`; return a correctly-shaped empty DataFrame or early-return before slicing.",
            "- `genal/MR_tools.py:MR_func()` always expects `Q/Q_df/Q_pval`, but several bootstrap/mode methods omit those keys and trigger a `KeyError`; add those fields (NaN) for all methods or make the final column selection conditional.",
            '- Deprecation warnings show invalid escape sequences in regex strings; use raw strings (e.g. `r"(\\\\d+)"`) to silence them.',
            "- `genal/MR_tools.py:apply_action_2()` triggers `SettingWithCopyWarning`; using `.copy()` after filtering or stricter `.loc[...]` assignment would avoid chained-assignment pitfalls.",
            "",
        ]
    )
    report_path.write_text(report_text)

    print(f"Wrote report: {report_path}")
    print(f"Wrote junit:  {junit_path}")
    print(f"Wrote output: {stdout_path}")

    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
