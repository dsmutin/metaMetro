"""Mandatory: strand-aware gold bubbles, their classes, and load_bench_gold_bubbles. No network."""

from __future__ import annotations

from pathlib import Path

import pytest

from metametro.bench import load_bench_gold_bubbles
from metametro.bench.bubble_gold import BUBBLE_COLUMNS, bubble_rows, strand_bubbles
from metametro.errors import ContractError
from metametro.formats.cdbg.model import Cdbg, Link, Unitig
from metametro.tables import write_tsv

pytestmark = pytest.mark.mandatory

_SEQ = {
    "S": "AAAA", "A": "CCCA", "B": "CCCT", "T": "GGGG",
    "C": "ACAC", "D": "AGAG", "E": "ATAT", "F": "TCTC",
    "G": "TGTG", "H": "CACC", "I": "GAGG", "J": "CTCC",
}


def _cdbg() -> Cdbg:
    arcs = [
        ("S", "A", "++"), ("S", "B", "++"), ("A", "T", "+-"), ("B", "T", "+-"),
        # The reverse strand of the first bubble, as MEGAHIT stores it.
        ("T", "A", "+-"), ("T", "B", "+-"), ("A", "S", "--"), ("B", "S", "--"),
        ("C", "D", "++"), ("C", "E", "++"), ("D", "F", "++"), ("E", "F", "++"),
        ("G", "H", "++"), ("G", "I", "++"), ("H", "J", "++"), ("I", "J", "++"),
    ]
    return Cdbg(
        metadata={},
        k=3,
        unitigs=[Unitig(f"u{name}", sequence, [f"n{name}"], []) for name, sequence in _SEQ.items()],
        links=[Link(f"e{index:02d}", f"u{a}", f"u{b}", orientation) for index, (a, b, orientation) in enumerate(arcs)],
        mapping=[],
    )


_SPECIES = {"nA": ("10665",), "nB": ("10726",), "nD": ("10665",), "nE": ("10665",), "nH": ("10760",)}


def test_each_bubble_is_found_once_per_strand_pair() -> None:
    bubbles = strand_bubbles(_cdbg())
    assert len(bubbles) == 3
    assert {frozenset(node[0] for path in bubble.branches for node in path) for bubble in bubbles} == {
        frozenset({"uA", "uB"}),
        frozenset({"uD", "uE"}),
        frozenset({"uH", "uI"}),
    }
    assert [bubble.bubble_id for bubble in bubbles] == [bubble.bubble_id for bubble in strand_bubbles(_cdbg())]


def test_bubble_classes_and_loader(tmp_path: Path) -> None:
    rows = bubble_rows(_cdbg(), _SPECIES, {"10665": "T4", "10726": "T5", "10760": "T7"})
    by_branches = {frozenset(row["branch_ids"].split(",")): row for row in rows}
    inter = by_branches[frozenset({"nA", "nB"})]
    assert inter["class"] == "inter_species" and inter["type"] == "strain" and inter["expected_action"] == "retain"
    assert inter["strain_ids"] == "10665,10726"
    assert by_branches[frozenset({"nD", "nE"})]["class"] == "intra_species"
    unresolved = by_branches[frozenset({"nH", "nI"})]
    assert unresolved["class"] == "unresolved" and unresolved["expected_action"] == ""
    write_tsv(tmp_path / "ground_truth" / "bubbles.tsv", list(BUBBLE_COLUMNS), rows)
    loaded = {bubble.branches: bubble for bubble in load_bench_gold_bubbles(tmp_path)}
    bubble = loaded[tuple((node,) for node in inter["branch_ids"].split(","))]
    assert bubble.bubble_class == "inter_species"
    assert set(bubble.branch_genomes) == {"T4", "T5"}


def test_loader_reads_the_older_bubble_file_and_needs_one(tmp_path: Path) -> None:
    with pytest.raises(ContractError, match="no gold bubbles"):
        load_bench_gold_bubbles(tmp_path)
    header = ["bubble_id", "type", "strain_ids", "source_id", "sink_id", "branch_ids", "expected_action"]
    write_tsv(
        tmp_path / "ground_truth" / "bubbles.tsv",
        header,
        [{"bubble_id": "B1", "type": "strain", "strain_ids": "taxon_1,taxon_2", "source_id": "S",
          "sink_id": "T", "branch_ids": "A,B", "expected_action": "retain"}],
    )
    (bubble,) = load_bench_gold_bubbles(tmp_path)
    assert bubble.branches == (("A",), ("B",))
    assert bubble.bubble_class == "" and bubble.branch_species == ()


def _toy(arcs: list[tuple[str, str]]) -> Cdbg:
    names = sorted({name for arc in arcs for name in arc})
    sequences = ["AAAC", "AAAG", "AACA", "AACC", "AACG", "AAGA", "AAGC", "AAGG"]
    return Cdbg(
        metadata={},
        k=3,
        unitigs=[Unitig(f"u{name}", sequences[index], [f"n{name}"], []) for index, name in enumerate(names)],
        links=[Link(f"e{index:02d}", f"u{a}", f"u{b}", "++") for index, (a, b) in enumerate(arcs)],
        mapping=[],
    )


def test_only_the_inner_bubble_of_a_nested_pair_is_gold() -> None:
    arcs = [("S", "A"), ("S", "B"), ("A", "X"), ("A", "Y"), ("X", "Z"), ("Y", "Z"), ("Z", "T"), ("B", "T")]
    bubbles = strand_bubbles(_toy(arcs))
    assert [(bubble.source, bubble.sink) for bubble in bubbles] == [(("uA", "+"), ("uZ", "+"))]
    assert {bubble.branches for bubble in bubbles} == {((("uX", "+"),), (("uY", "+"),))}


def test_a_direct_link_is_an_empty_unresolved_branch() -> None:
    cdbg = _toy([("S", "A"), ("S", "T"), ("A", "T")])
    [bubble] = strand_bubbles(cdbg)
    assert bubble.branches == ((("uA", "+"),), ())
    [row] = bubble_rows(cdbg, {"nA": ("10665",)}, {"10665": "T4"})
    assert row["branch_ids"] == "nA,"
    assert row["class"] == "unresolved" and row["type"] == "unresolved" and row["expected_action"] == ""
    assert row["branch_species"] == "10665,"


def test_bubble_class_comes_from_the_branch_membership() -> None:
    cdbg = _toy([("S", "A"), ("S", "B"), ("A", "T"), ("B", "T")])
    genome_of = {"10665": "T4", "10760": "T7", "10726": "T5"}
    [shared] = bubble_rows(cdbg, {"nA": ("10665", "10760"), "nB": ("10665", "10760")}, genome_of)
    assert shared["class"] == "intra_species" and shared["branch_species"] == "10665;10760,10665;10760"
    assert shared["branch_genomes"] == "T4;T7,T4;T7"
    [differs] = bubble_rows(cdbg, {"nA": ("10665", "10760"), "nB": ("10665",)}, genome_of)
    assert differs["class"] == "inter_species" and differs["expected_action"] == "retain"
    assert differs["strain_ids"] == "10665,10760"
    [missing] = bubble_rows(cdbg, {"nA": ("10665", "10760")}, genome_of)
    assert missing["class"] == "unresolved" and missing["branch_species"] == "10665;10760,"


def _overlapping_bubble(tmp_path: Path) -> Path:
    """Genomes A = P+X+Q and B = P+Y+Q; CFA nodes overlap by 21 bp, as MEGAHIT's do (k = 22)."""
    import random

    from metametro.bench.phage_gold import node_membership, overlapped_ends, write_taxon_membership
    from metametro.converters.cfa_to_cdbg import cfa_to_cdbg
    from metametro.formats.cdbg.io import dump_cdbg
    from metametro.formats.cfa.model import CfaGraph

    rng = random.Random(7)
    P, X, Y, Q = ("".join(rng.choice("ACGT") for _ in range(n)) for n in (40, 30, 30, 40))
    genomes = {"A": P + X + Q, "B": P + Y + Q}
    sequences = {"S": P, "a": P[-21:] + X + Q[:21], "b": P[-21:] + Y + Q[:21], "T": Q}
    edges = [
        {"edge_id": f"e{i}", "source": s, "target": t, "orientation": "++"}
        for i, (s, t) in enumerate([("S", "a"), ("S", "b"), ("a", "T"), ("b", "T")])
    ]
    graph = CfaGraph(
        metadata={
            "schema_version": "1.0",
            "graph_id": "toy",
            "graph_type": "de_bruijn",
            "k": 22,
            "features": {"edge": {"orientation": "orientation"}},
        },
        sequences=sequences,
        nodes=[{"node_id": name} for name in sequences],
        edges=edges,
        node_header=["node_id"],
        edge_header=["edge_id", "source", "target", "orientation"],
    )
    dump_cdbg(cfa_to_cdbg(graph), tmp_path / "cdbg")
    write_tsv(
        tmp_path / "ground_truth" / "taxon_values.tsv",
        ["taxon", "genome"],
        [{"taxon": "1", "genome": "A"}, {"taxon": "2", "genome": "B"}],
    )
    old = node_membership(sequences, genomes, 21)
    assert old.members("a") == ("A", "B") and old.members("b") == ("A", "B")
    new = node_membership(sequences, genomes, 21, ends=overlapped_ends(edges), overlap=21)
    assert new.members("a") == ("A",) and new.members("b") == ("B",)
    assert new.interior_positions["a"] == old.positions["a"] - 2
    assert new.owner("a") == old.owner("a") == "A"
    write_taxon_membership(tmp_path, new)
    return tmp_path


def test_overlap_kmers_do_not_make_a_branch_a_member_of_its_neighbours_genome(tmp_path: Path) -> None:
    from metametro.bench.bubble_gold import write_phage_bubbles
    from metametro.tables import read_tsv

    write_phage_bubbles(_overlapping_bubble(tmp_path))
    _, rows = read_tsv(tmp_path / "ground_truth" / "bubbles.tsv")
    assert [(row["class"], row["branch_genomes"], row["expected_action"]) for row in rows] == [
        ("inter_species", "A,B", "retain")
    ]


def test_a_non_genomic_branch_beside_a_genomic_one_is_an_error() -> None:
    cdbg = _toy([("S", "A"), ("S", "B"), ("A", "T"), ("B", "T")])
    [row] = bubble_rows(cdbg, {"nA": ("10665",)}, {"10665": "T4"}, informative={"nA", "nB"})
    assert (row["class"], row["type"], row["expected_action"]) == ("error", "error", "pop")
    assert row["branch_ids"] == "nA,nB" and row["branch_genomes"] == "T4,"
    [row] = bubble_rows(cdbg, {"nB": ("10665",)}, {"10665": "T4"}, informative={"nA", "nB"})
    assert row["class"] == "error"
    assert (row["branch_ids"], row["branch_genomes"], row["branch_species"]) == ("nB,nA", "T4,", "10665,")
    [row] = bubble_rows(cdbg, {"nA": ("10665",)}, {"10665": "T4"}, informative={"nA"})
    assert row["class"] == "unresolved"
    [row] = bubble_rows(cdbg, {}, {}, informative={"nA", "nB"})
    assert row["class"] == "unresolved"


def _variant_bubble(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, branches: dict[str, str], locus: str
) -> Path:
    """One genome P+locus+Q; each branch node is P[-21:]+variant+Q[:21] (21-bp overlaps, k = 22).

    Membership comes from ``phage_membership`` as in a build, with the pinned
    genome list and the genome files replaced by ``G``.
    """
    import random
    from types import SimpleNamespace

    from metametro.bench import phage_gold
    from metametro.converters.cfa_to_cdbg import cfa_to_cdbg
    from metametro.formats.cdbg.io import dump_cdbg
    from metametro.formats.cfa.model import CfaGraph

    rng = random.Random(11)
    P, Q = ("".join(rng.choice("ACGT") for _ in range(40)) for _ in range(2))
    genomes = {"G": P + locus + Q}
    sequences = {"S": P, "T": Q, **{name: P[-21:] + variant + Q[:21] for name, variant in branches.items()}}
    arcs = [("S", name) for name in branches] + [(name, "T") for name in branches]
    edges = [{"edge_id": f"e{i}", "source": s, "target": t, "orientation": "++"} for i, (s, t) in enumerate(arcs)]
    graph = CfaGraph(
        metadata={
            "schema_version": "1.0",
            "graph_id": "toy",
            "graph_type": "de_bruijn",
            "k": 22,
            "features": {"edge": {"orientation": "orientation"}},
        },
        sequences=sequences,
        nodes=[{"node_id": name} for name in sequences],
        edges=edges,
        node_header=["node_id"],
        edge_header=["edge_id", "source", "target", "orientation"],
    )
    dump_cdbg(cfa_to_cdbg(graph), tmp_path / "cdbg")
    write_tsv(tmp_path / "ground_truth" / "taxon_values.tsv", ["taxon", "genome"], [{"taxon": "1", "genome": "G"}])
    monkeypatch.setattr(phage_gold, "pinned_genomes", lambda bench: [SimpleNamespace(name="G")])
    monkeypatch.setattr(phage_gold, "load_genomes", lambda folder: genomes)
    membership = phage_gold.phage_membership(graph, "phage_species_5_x10", tmp_path, 21)
    phage_gold.write_taxon_membership(tmp_path, membership)
    return tmp_path


def _mutate(sequence: str, *positions: int) -> str:
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    letters = list(sequence)
    for position in positions:
        letters[position] = swap[letters[position]]
    return "".join(letters)


def test_two_error_variants_of_one_locus_are_unresolved_not_intra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import random

    from metametro.bench.bubble_gold import write_phage_bubbles
    from metametro.bench.phage_gold import read_taxon_membership
    from metametro.tables import read_tsv

    locus = "".join(random.Random(3).choice("ACGT") for _ in range(30))
    bench = _variant_bubble(tmp_path, monkeypatch, {"a": _mutate(locus, 5, 25), "b": _mutate(locus, 8, 22)}, locus)
    assert read_taxon_membership(bench)["a"] == ("G",) and read_taxon_membership(bench)["b"] == ("G",)
    write_phage_bubbles(bench)
    _, rows = read_tsv(bench / "ground_truth" / "bubbles.tsv")
    assert [(row["class"], row["branch_genomes"], row["expected_action"]) for row in rows] == [("unresolved", ",", "")]


def test_the_error_branch_is_listed_last_as_bubbleblower_reads_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import random

    from metametro.bench.bubble_gold import write_phage_bubbles
    from metametro.tables import read_tsv

    locus = "".join(random.Random(3).choice("ACGT") for _ in range(30))
    # "a" is an error variant; "b" is the genome sequence. Sorted branch order would put "a" first.
    bench = _variant_bubble(tmp_path, monkeypatch, {"a": _mutate(locus, 5, 25), "b": locus}, locus)
    write_phage_bubbles(bench)
    _, rows = read_tsv(bench / "ground_truth" / "bubbles.tsv")
    [row] = rows
    assert (row["class"], row["type"], row["expected_action"]) == ("error", "error", "pop")
    assert row["branch_ids"] == "b,a"
    assert row["branch_genomes"] == "G,"


def test_phage_membership_leaves_out_the_overlaps_of_a_built_graph(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A SNP branch of a 1-bp locus shares only its two overlap k-mers with the genome."""
    from metametro.bench.bubble_gold import write_phage_bubbles
    from metametro.bench.phage_gold import read_taxon_membership
    from metametro.tables import read_tsv

    bench = _variant_bubble(tmp_path, monkeypatch, {"a": "C", "b": "A"}, "A")
    assert read_taxon_membership(bench)["a"] == () and read_taxon_membership(bench)["b"] == ("G",)
    write_phage_bubbles(bench)
    _, rows = read_tsv(bench / "ground_truth" / "bubbles.tsv")
    assert [(row["class"], row["branch_ids"], row["branch_genomes"]) for row in rows] == [("error", "b,a", "G,")]
