"""Mandatory: CLI entry works."""

from __future__ import annotations

import json
import subprocess
import sys

import pytest

from metametro import __version__
from metametro.cli import main

pytestmark = pytest.mark.mandatory


def test_cli_version(capsys: pytest.CaptureFixture[str]) -> None:
    """--version prints the package version."""
    assert main(["--version"]) == 0
    out = capsys.readouterr().out.strip()
    assert out


def test_cli_baseline_stdout(capsys: pytest.CaptureFixture[str]) -> None:
    """Default CLI run prints JSON with baseline status."""
    assert main([]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "baseline"
    assert payload["ok"] is True


def test_python_m_metametro_version() -> None:
    """``python -m metametro --version`` prints the package version."""
    done = subprocess.run([sys.executable, "-m", "metametro", "--version"], capture_output=True, text=True, check=True)
    assert done.stdout.strip() == __version__


def test_cli_help_mentions_benchbuild(capsys: pytest.CaptureFixture[str]) -> None:
    """--help names the benchbuild subcommand."""
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert "benchbuild" in capsys.readouterr().out
