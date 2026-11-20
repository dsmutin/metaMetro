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


def _row(**changes: object) -> dict:
    from metametro.bench import canonical_identity

    row = {
        "bench": "phage_10",
        "canonical_name": "phage_species_5_x10",
        "build_status": "built",
        "solved": True,
        "colourings": ["as_built", "composition_kmeans", "kraken2", "decaying"],
        "kraken2_input": "metametro_pin",
        "metrics": {"f1": 0.0},
        "identity_as_built": canonical_identity("phage_10"),
    }
    row.update(changes)
    return {key: value for key, value in row.items() if value is not None or key == "metrics"}


def _summary(tmp_path: Path, row: dict, ok: bool) -> Path:
    path = tmp_path / "solve.json"
    other = {"bench": "mock_bubble", "build_status": "built", "solved": True}
    path.write_text(json.dumps({"ok": ok, "benches": [other, row]}), encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("changes", "code"),
    [
        ({}, 0),
        ({"build_status": "skipped", "solved": False}, 1),
        ({"solved": "true"}, 1),
        ({"colourings": ["composition_kmeans", "decaying"], "colourings_skipped": ["kraken2"]}, 1),
        ({"colourings_skipped": ["kaiju"]}, 1),
        ({"kraken2_input": "given"}, 1),
        ({"metrics": None}, 1),
        ({"identity_as_built": "0" * 64}, 1),
        ({"identity_canonical": "0" * 64}, 1),
    ],
)
def test_solve_row_decides_not_the_ok_flag(tmp_path: Path, changes: dict, code: int) -> None:
    assert _checker().check_solve(_summary(tmp_path, _row(**changes), ok=code != 0)) == code


def test_solve_without_identity_as_built_reads_the_solved_build(tmp_path: Path) -> None:
    from metametro.bench import canonical_gold
    from metametro.tables import write_yaml

    checker = _checker()
    colourings = ["kraken2", "composition_kmeans", "decaying"]
    row = _row(identity_as_built=None, colourings=colourings)
    summary = _summary(tmp_path, row, ok=True)
    assert checker.check_solve(summary) == 1
    bench = tmp_path / "work" / "phage_species_5_x10"
    bench.mkdir(parents=True)
    expected = checker.SOLVE_IDENTITIES[tuple(colourings)]
    (bench / "identity.sha256").write_text(f"{expected}  identity\n", encoding="utf-8")
    write_yaml(bench / "manifest.yaml", {"gold_sha256": canonical_gold("phage_10")})
    assert checker.check_solve(summary, tmp_path / "work") == 0
    (bench / "identity.sha256").write_text(f"{'0' * 64}  identity\n", encoding="utf-8")
    assert checker.check_solve(summary, tmp_path / "work") == 1


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
