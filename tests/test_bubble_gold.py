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
