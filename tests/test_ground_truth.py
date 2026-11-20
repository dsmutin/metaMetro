"""Mandatory checks that gold taxon labels are ground truth, not a training colour."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metametro.bench import load_bench_cgt, load_bench_gold, load_bench_gold_cgt
from metametro.bench.build import build
from metametro.bench.data.universal.bubbles import strain_bubble
from metametro.bench.leakage import assert_no_target_leak
from metametro.contracts.colouring import paint_namespace
from metametro.contracts.decaying import colour_decaying
from metametro.errors import ContractError
from metametro.formats.cfa.io import load_cfa

pytestmark = pytest.mark.mandatory


def test_guard_rejects_taxon_colour() -> None:
    """A CFA that still carries the gold ``taxon`` namespace fails the leakage guard."""
    with pytest.raises(ContractError, match="taxon"):
        assert_no_target_leak(strain_bubble())


def test_decaying_never_derives_from_taxon() -> None:
    """Decaying reads composition, not gold taxon, and refuses taxon as a source."""
    graph = paint_namespace(
        strain_bubble(),
        {"S": ["cluster_0"], "A": ["cluster_1"], "B": ["cluster_1"], "T": ["cluster_0"]},
        namespace="composition",
        operation="merge",
    )
    leaked = colour_decaying(graph, decay=0.5, iterations=4, threshold=0.2)
    values = {row["value"] for row in leaked.colors or [] if row["namespace"] == "decaying"}
    assert values and values <= {"cluster_0", "cluster_1"}
    with pytest.raises(ContractError, match="taxon"):
        colour_decaying(graph, source="taxon")
    with pytest.raises(ContractError, match="non-gold"):
        colour_decaying(strain_bubble())


def test_bubble_bench_decaying_uses_composition(tmp_path: Path) -> None:
    """On a bubble bench (no sample layer) decaying values are composition clusters."""
    result = build("bubble_strain_2", outdir=tmp_path / "out")
    cfa = load_cfa(result.outdir / "cfa")
    assert "taxon" not in {row["namespace"] for row in cfa.colors or []}
    decaying = {row["value"] for row in cfa.colors or [] if row["namespace"] == "decaying"}
    assert decaying and all(value.startswith("cluster_") for value in decaying)


def test_gold_labels_available_from_ground_truth(tmp_path: Path) -> None:
    """Gold taxon labels are read from ``ground_truth/`` per CFA id and aligned with the CGT."""
    result = build("bubble_strain_2", outdir=tmp_path / "out")
    assert (result.outdir / "ground_truth" / "taxon_nodes.tsv").is_file()
    gold = load_bench_gold(result.outdir)
    assert gold.values == ("taxon_1", "taxon_2")
    assert gold.nodes == {"S": ("taxon_1", "taxon_2"), "A": ("taxon_1",), "B": ("taxon_2",), "T": ("taxon_1", "taxon_2")}
    assert gold.edges["eSB"] == ("taxon_2",)
    labels = load_bench_gold_cgt(result.outdir)
    cgt = load_bench_cgt(result.outdir)
    assert labels.num_nodes == cgt.num_nodes
    assert np.array_equal(labels.indices, cgt.indices)
    assert labels.node_colors.shape == (cgt.num_nodes, 2)
    assert labels.edge_colors.shape == (cgt.edge_colors.shape[0], 2)
    assert labels.node_features.shape[1] == 0
    assert labels.node_colors.sum(axis=0).tolist() == [3, 3]
