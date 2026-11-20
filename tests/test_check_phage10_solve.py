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


def test_identity_check_uses_canonical_identity_and_gold(tmp_path: Path) -> None:
    from metametro.bench import canonical_gold, canonical_identity
    from metametro.tables import write_yaml

    gold = tmp_path / "ground_truth" / "taxon_nodes.tsv"
    gold.parent.mkdir()
    gold.write_text("node_id\ttaxon\n", encoding="utf-8")
    from metametro.bench.ground_truth import write_gold_digest

    digest = write_gold_digest(tmp_path)
    (tmp_path / "identity.sha256").write_text(f"{canonical_identity('phage_10')}  identity\n", encoding="utf-8")
    write_yaml(tmp_path / "manifest.yaml", {"gold_sha256": digest})
    assert _checker().check_identity(tmp_path) == 1
    write_yaml(tmp_path / "manifest.yaml", {"gold_sha256": canonical_gold("phage_10")})
    assert _checker().check_identity(tmp_path) == 1


def test_identity_check_passes_on_the_built_phage_10() -> None:
    import os

    folder = os.environ.get("METAMETRO_PHAGE10_DIR")
    if not folder or not (Path(folder) / "identity.sha256").is_file():
        pytest.skip("set METAMETRO_PHAGE10_DIR to a built phage_10 (the phage10-solve CI job does)")
    assert _checker().check_identity(Path(folder)) == 0
