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


def test_loading_gold_checks_gold_sha256(tmp_path: Path) -> None:
    from metametro.bench import load_bench_gold, load_bench_gold_bubbles
    from metametro.errors import ContractError

    _build(tmp_path, "value\nA\n")
    assert load_bench_gold(tmp_path).nodes
    assert load_bench_gold_bubbles(tmp_path)
    (tmp_path / "ground_truth" / "extra.tsv").write_text("value\nB\n", encoding="utf-8")
    with pytest.raises(ContractError, match="extra.tsv.*rebuild the bench"):
        load_bench_gold(tmp_path)
    with pytest.raises(ContractError, match="does not match gold.sha256"):
        load_bench_gold_bubbles(tmp_path)
    (tmp_path / "gold.sha256").unlink()
    assert load_bench_gold(tmp_path).nodes


def _stale_phage_10(folder: Path) -> Path:
    """A self-consistent phage_10 ground_truth/ as 0.25.0 wrote it: 94 intra_species bubbles."""
    from metametro.bench.bubble_gold import BUBBLE_COLUMNS
    from metametro.bench.ground_truth import write_gold_digest
    from metametro.tables import write_tsv, write_yaml

    rows = [
        {column: "" for column in BUBBLE_COLUMNS}
        | {"bubble_id": f"b{i:02d}", "type": "intra_species", "class": "intra_species", "branch_ids": "n1,n2"}
        for i in range(94)
    ]
    write_tsv(folder / "ground_truth" / "bubbles.tsv", list(BUBBLE_COLUMNS), rows)
    write_tsv(folder / "ground_truth" / "taxon_nodes.tsv", ["node_id", "taxon"], [{"node_id": "n1", "taxon": "10665"}])
    write_tsv(folder / "ground_truth" / "taxon_edges.tsv", ["edge_id", "source", "target", "taxon"], [])
    digest = write_gold_digest(folder)
    write_yaml(folder / "manifest.yaml", {"bench": "phage_species_5_x10", "status": "built", "gold_sha256": digest})
    (folder / "identity.sha256").write_text(
        "78121c1b05a5814bfc17ac3b7650ed5223a0f584cfada24e9224b7ba0364564c  identity\n", encoding="utf-8"
    )
    return folder


def test_stale_self_consistent_phage_gold_is_refused_on_load_and_reuse(tmp_path: Path) -> None:
    from metametro.bench import load_bench_gold, load_bench_gold_bubbles, resolve
    from metametro.bench.build import build_all
    from metametro.bench.ground_truth import check_gold_digest
    from metametro.bench.paths import default_outdir
    from metametro.errors import ContractError

    stale = _stale_phage_10(tmp_path / "p10")
    for load in (check_gold_digest, load_bench_gold, load_bench_gold_bubbles):
        with pytest.raises(ContractError, match="not the canonical phage_species_5_x10 gold"):
            load(stale)
    (stale / "gold.sha256").unlink()
    with pytest.raises(ContractError, match="has no gold.sha256"):
        load_bench_gold_bubbles(stale)
    reused = _stale_phage_10(default_outdir(resolve("phage_10"), root=tmp_path / "root"))
    with pytest.raises(ContractError, match="not the canonical"):
        build_all(root=tmp_path / "root")
    assert reused.is_dir()
