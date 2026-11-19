"""Mandatory: version is defined once in VERSION."""

from __future__ import annotations

from pathlib import Path

import pytest

from metametro import __version__

pytestmark = pytest.mark.mandatory


def test_version_file_matches_package(repo_root: Path) -> None:
    """Package __version__ must equal the VERSION file."""
    written = (repo_root / "VERSION").read_text(encoding="utf-8").strip()
    assert written == __version__
    assert written == "0.0.1" or len(written.split(".")) == 3


def test_readme_mandatory_count_matches_collection(repo_root: Path) -> None:
    """README states the number of tests that ``pytest -m mandatory`` collects."""
    import re
    import subprocess
    import sys

    readme = (repo_root / "README.md").read_text(encoding="utf-8")
    match = re.search(r"The mandatory marker selects (\d+) tests", readme)
    assert match, "README must state the mandatory test count"
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-m", "mandatory", "--collect-only", "-q", "-o", "addopts=--strict-markers", "-p", "no:cacheprovider"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    collected = sum(1 for line in result.stdout.splitlines() if "::" in line)
    assert collected > 0
    assert int(match.group(1)) == collected


def test_readme_python_example_runs(repo_root: Path) -> None:
    """The python block in README "How it works" runs without errors."""
    import re

    readme = (repo_root / "README.md").read_text(encoding="utf-8")
    blocks = re.findall(r"```python\n(.*?)```", readme, flags=re.S)
    assert blocks, "README must contain a python example"
    for block in blocks:
        exec(compile(block, "README.md", "exec"), {})
