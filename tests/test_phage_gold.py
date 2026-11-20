"""Mandatory: phage gold taxon labels come from the pinned genomes and land in ground_truth/. No network."""

from __future__ import annotations

from pathlib import Path

import pytest

from metametro.bench import load_bench_gold
from metametro.bench.ground_truth import split_gold, write_gold
from metametro.bench.leakage import assert_no_target_leak
from metametro.bench.phage_gold import node_genomes, paint_phage_gold, write_taxon_values
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import read_tsv

pytestmark = pytest.mark.mandatory

_GENOMES = {
    "T1": "AAACCCGGGTTTACGA",
    "T3": "CATGCATGCAAGGTTC",
    "T4": "GTCAGTCATTGACCAA",
    "T5": "TTGGCCAATCGATCGA",
    "T7": "AGCTAGCTTTCAGGCA",
}


def _graph() -> CfaGraph:
    return CfaGraph(
        metadata={
            "schema_version": "1.0",
            "graph_id": "phage_toy",
            "graph_type": "de_bruijn",
            "k": 5,
            "features": {"node": {"color_set": "color_set"}, "edge": {"color_set": "color_set"}},
        },
        sequences={"a": "AAACCCGG", "b": "CCGGGTTT", "c": "GCATGCAA", "d": "NNNNNNNN"},
        nodes=[
            {"node_id": "a", "color_set": "0"},
            {"node_id": "b", "color_set": "0"},
            {"node_id": "c", "color_set": "1"},
            {"node_id": "d", "color_set": ""},
        ],
        edges=[
            {"edge_id": "ab", "source": "a", "target": "b", "color_set": "0"},
            {"edge_id": "bc", "source": "b", "target": "c", "color_set": ""},
        ],
        node_header=["node_id", "color_set"],
        edge_header=["edge_id", "source", "target", "color_set"],
        colors=[
            {"color_id": "0", "namespace": "sample", "value": "sample_1"},
            {"color_id": "1", "namespace": "sample", "value": "sample_2"},
        ],
    )


def _genome_dir(tmp_path: Path) -> Path:
    root = tmp_path / "genomes"
    root.mkdir()
    for name, sequence in _GENOMES.items():
        (root / f"{name}.fna").write_text(f">{name}\n{sequence}\n", encoding="utf-8")
    return root


def test_node_genomes_needs_a_strict_majority() -> None:
    owners = node_genomes(_graph().sequences, _GENOMES, 5)
    assert owners == {"a": "T1", "b": "T1", "c": "T3", "d": ""}


def test_phage_gold_is_species_taxid_in_ground_truth_only(tmp_path: Path) -> None:
    painted = paint_phage_gold(_graph(), "phage_10", _genome_dir(tmp_path), 5)
    stripped, gold = split_gold(painted)
    assert_no_target_leak(stripped)
    assert stripped.colors == _graph().colors
    assert [row["color_set"] for row in stripped.nodes] == ["0", "0", "1", ""]
    write_gold(tmp_path / "bench", gold, stripped)
    labels = load_bench_gold(tmp_path / "bench")
    # Species taxids of T1 (Tunavirus T1) and T3 (Teetrevirus T3) from genomes.tsv.
    assert labels.nodes == {"a": ("1921008",), "b": ("1921008",), "c": ("2732706",), "d": ()}
    assert labels.edges == {"ab": ("1921008",), "bc": ()}
    path = write_taxon_values(tmp_path / "bench", "phage_10")
    header, rows = read_tsv(path)
    assert header == ["taxon", "rank", "name", "genome", "accession", "record_taxid"]
    by_taxon = {row["taxon"]: row for row in rows}
    assert by_taxon["1921008"]["rank"] == "species"
    assert by_taxon["1921008"]["genome"] == "T1"
    assert by_taxon["1921008"]["record_taxid"] == "2492962"
    assert len(rows) == 5


def test_membership_is_deterministic_and_a_shared_node_gets_both_genomes() -> None:
    from metametro.bench.phage_gold import node_genomes, node_membership

    genomes = {"T7": "AAAACCCCGGGG", "T3": "CCCCGGGGTATA", "T4": "ACGTTGCA"}
    sequences = {"shared": "CCCCGGGG", "t7": "AAAACC", "short": "AC", "none": "GAGAGA"}
    first = node_membership(sequences, genomes, 4)
    again = node_membership(dict(reversed(sequences.items())), dict(reversed(genomes.items())), 4)
    assert first == again
    assert first.genomes == ("T3", "T4", "T7")
    assert first.members("shared") == ("T3", "T7")
    assert first.counts["shared"] == (5, 0, 5) and first.positions["shared"] == 5
    assert first.owner("shared") == ""
    assert first.members("t7") == ("T7",) and first.owner("t7") == "T7"
    assert first.members("short") == () and first.positions["short"] == 0
    assert first.members("none") == ()
    assert node_genomes(sequences, genomes, 4) == {"shared": "", "t7": "T7", "short": "", "none": ""}


def test_membership_file_round_trip(tmp_path) -> None:
    from metametro.bench.phage_gold import node_membership, read_taxon_membership, write_taxon_membership

    membership = node_membership({"n1": "CCCCGGGG", "n2": "AAAACC"}, {"T7": "AAAACCCCGGGG", "T3": "CCCCGGGGTATA"}, 4)
    path = write_taxon_membership(tmp_path, membership)
    assert path.read_text().splitlines() == [
        "node_id\tkmers\tinterior_kmers\tT3\tT3_fraction\tT7\tT7_fraction",
        "n1\t5\t5\t1\t1.000000\t1\t1.000000",
        "n2\t3\t3\t0\t0.000000\t1\t1.000000",
    ]
    assert read_taxon_membership(tmp_path) == {"n1": ("T3", "T7"), "n2": ("T7",)}
