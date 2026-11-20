"""Mandatory: gold.sha256 covers ground_truth/ and is separate from identity.sha256. No network."""

from __future__ import annotations

from pathlib import Path

import pytest

from metametro.bench.colourings import BuildContext
from metametro.bench.ground_truth import gold_digest
from metametro.bench.materialize import materialize
from metametro.bench.registry import PREPARE
from metametro.tables import read_yaml

pytestmark = pytest.mark.mandatory


def _build(outdir: Path, extra: str) -> str:
    def write(path: Path) -> None:
        (path / "ground_truth" / "extra.tsv").write_text(extra, encoding="utf-8")

    ctx = BuildContext()
    graph = PREPARE["bubble_strain_2"](outdir, ctx)
    return materialize(graph, outdir, ctx, before_identity=write).identity


def test_gold_digest_changes_with_a_gold_file_and_identity_does_not(tmp_path: Path) -> None:
    first = _build(tmp_path / "a", "value\nA\n")
    second = _build(tmp_path / "b", "value\nB\n")
    assert first == second
    gold_a = read_yaml(tmp_path / "a" / "manifest.yaml")["gold_sha256"]
    gold_b = read_yaml(tmp_path / "b" / "manifest.yaml")["gold_sha256"]
    assert gold_a != gold_b
    assert gold_a == gold_digest(tmp_path / "a")
    lines = (tmp_path / "a" / "gold.sha256").read_text().splitlines()
    assert [line.split("  ")[1] for line in lines] == sorted(
        f"ground_truth/{path.name}" for path in (tmp_path / "a" / "ground_truth").iterdir()
    )
    (tmp_path / "a" / "ground_truth" / "taxon_nodes.tsv").write_text("node_id\ttaxon\n", encoding="utf-8")
    assert gold_digest(tmp_path / "a") != gold_a
    assert (tmp_path / "a" / "identity.sha256").read_text().split()[0] == first
