"""Mandatory: the phage10-solve CI checker trusts only the phage_10 row. No network."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.mandatory

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_phage10_solve.py"


def _checker():
    spec = importlib.util.spec_from_file_location("check_phage10_solve", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("row", "code"),
    [
        ({"bench": "phage_10", "build_status": "built", "solved": True}, 0),
        ({"bench": "phage_10", "build_status": "skipped", "solved": False}, 1),
        ({"bench": "phage_10", "build_status": "built", "solved": "true"}, 1),
    ],
)
def test_solve_row_decides_not_the_ok_flag(tmp_path: Path, row: dict, code: int) -> None:
    summary = tmp_path / "solve.json"
    other = {"bench": "mock_bubble", "build_status": "built", "solved": True}
    summary.write_text(json.dumps({"ok": code != 0, "benches": [other, row]}), encoding="utf-8")
    assert _checker().check_solve(summary) == code


def test_identity_check_uses_canonical_identity(tmp_path: Path) -> None:
    from metametro.bench import canonical_identity

    (tmp_path / "identity.sha256").write_text(f"{canonical_identity('phage_10')}  identity\n", encoding="utf-8")
    assert _checker().check_identity(tmp_path) == 0
    (tmp_path / "identity.sha256").write_text(f"{'0' * 64}  identity\n", encoding="utf-8")
    assert _checker().check_identity(tmp_path) == 1
