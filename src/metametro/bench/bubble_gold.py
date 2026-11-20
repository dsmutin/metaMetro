"""Gold bubble labels of an assembled bench, written to ``ground_truth/bubbles.tsv``.

A bubble is a simple bubble on the strand-aware CDBG: two or more internally
disjoint non-branching paths that leave one oriented unitig and enter one
oriented unitig, at most ``max_internal`` unitigs each (the definition of
BubbleBlower ``detect_bubbles``, with each link read as an arc between oriented
unitigs). MEGAHIT graphs store both strands of every adjacency as separate
links, so every bubble is found twice, once per strand; the copy with the
smaller ``bubble_id`` is kept. ``bubble_id`` is BubbleBlower's sequence
signature (the first 16 hex digits of SHA-256 over the source, sink, and branch
sequence hashes), computed on the oriented sequences.

A branch is labelled with the gold species shared by all its CFA nodes, or is
unlabelled. ``class`` is ``inter_species`` when the labelled branches map to
different species, ``intra_species`` when they map to one, and ``unresolved``
when any branch is unlabelled. ``type`` is ``strain`` and ``expected_action``
``retain`` for ``inter_species``; otherwise ``type`` is the class and
``expected_action`` is empty (the hand-built benches use ``strain``/``error``
and ``retain``/``pop``). Nested or overlapping bubbles and branches longer than
``max_internal`` unitigs are not gold.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from metametro.formats.cdbg.model import Cdbg
from metametro.tables import read_tsv, write_tsv

BUBBLE_COLUMNS = (
    "bubble_id",
    "type",
    "strain_ids",
    "source_id",
    "sink_id",
    "branch_ids",
    "expected_action",
    "class",
    "branch_genomes",
    "branch_species",
)
_FLIP = {"+": "-", "-": "+"}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

Oriented = tuple[str, str]


@dataclass(frozen=True)
class StrandBubble:
    """One simple bubble between oriented unitigs."""

    bubble_id: str
    source: Oriented
    sink: Oriented
    branches: tuple[tuple[Oriented, ...], ...]


def _oriented_sequence(sequences: dict[str, str], node: Oriented) -> str:
    sequence = sequences[node[0]]
    return sequence if node[1] == "+" else sequence.translate(_COMPLEMENT)[::-1]


def _hash(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def _signature(sequences: dict[str, str], source: Oriented, sink: Oriented, branches) -> str:
    tokens = [_hash(_oriented_sequence(sequences, source)), _hash(_oriented_sequence(sequences, sink))]
    for path, links in branches:
        text = "".join(_oriented_sequence(sequences, node) for node in path) if path else links[0]
        tokens.append(_hash(text))
    return hashlib.sha256("|".join(tokens).encode("ascii")).hexdigest()[:16]


def strand_bubbles(cdbg: Cdbg, *, max_internal: int = 32) -> list[StrandBubble]:
    """Return the simple bubbles of ``cdbg``, one per strand pair, sorted by ``bubble_id``."""
    sequences = {unitig.unitig_id: unitig.sequence for unitig in cdbg.unitigs}
    outgoing: dict[Oriented, list[tuple[str, Oriented]]] = {}
    indeg: dict[Oriented, int] = {}
    outdeg: dict[Oriented, int] = {}
    for link in cdbg.links:
        orientation = link.orientation or "++"
        source, target = (link.source, orientation[0]), (link.target, orientation[1])
        outgoing.setdefault(source, []).append((link.link_id, target))
        outdeg[source] = outdeg.get(source, 0) + 1
        indeg[target] = indeg.get(target, 0) + 1
    for arcs in outgoing.values():
        arcs.sort()
    found: dict[tuple, StrandBubble] = {}
    for source in sorted(outgoing):
        arcs = outgoing[source]
        if len(arcs) < 2:
            continue
        grouped: dict[Oriented, list[tuple[tuple[Oriented, ...], tuple[str, ...]]]] = {}
        for link_id, node in arcs:
            path: list[Oriented] = []
            links = [link_id]
            seen = {source}
            closed = True
            while True:
                if node in seen:
                    closed = False
                    break
                if indeg.get(node, 0) != 1 or outdeg.get(node, 0) != 1:
                    break
                if len(path) >= max_internal:
                    closed = False
                    break
                seen.add(node)
                path.append(node)
                links.append(outgoing[node][0][0])
                node = outgoing[node][0][1]
            if closed and node != source:
                grouped.setdefault(node, []).append((tuple(path), tuple(links)))
        for sink, branches in grouped.items():
            unique = list(dict.fromkeys(sorted(branches, key=lambda item: item[1])))
            if len(unique) < 2:
                continue
            internal = [node[0] for path, _ in unique for node in path]
            if len(internal) != len(set(internal)) or {source[0], sink[0]} & set(internal):
                continue
            bubble = StrandBubble(
                bubble_id=_signature(sequences, source, sink, unique),
                source=source,
                sink=sink,
                branches=tuple(path for path, _ in unique),
            )
            key = (
                frozenset((source[0], sink[0])),
                frozenset(frozenset(node[0] for node in path) for path in bubble.branches),
            )
            if key not in found or bubble.bubble_id < found[key].bubble_id:
                found[key] = bubble
    return sorted(found.values(), key=lambda item: item.bubble_id)


def bubble_rows(
    cdbg: Cdbg, node_species: dict[str, tuple[str, ...]], genome_of: dict[str, str]
) -> list[dict[str, str]]:
    """Gold rows for every bubble of ``cdbg``. ``node_species`` is keyed by CFA node id."""
    members = {unitig.unitig_id: list(unitig.members) for unitig in cdbg.unitigs}
    rows: list[dict[str, str]] = []
    for bubble in strand_bubbles(cdbg):
        tokens: list[str] = []
        labels: list[str] = []
        for path in bubble.branches:
            nodes = [member for node in path for member in members[node[0]]]
            tokens.append(";".join(nodes))
            found = {node_species.get(node, ()) for node in nodes}
            labels.append(next(iter(found))[0] if len(found) == 1 and len(next(iter(found))) == 1 else "")
        if any(label == "" for label in labels):
            kind = "unresolved"
        elif len(set(labels)) == 1:
            kind = "intra_species"
        else:
            kind = "inter_species"
        rows.append(
            {
                "bubble_id": bubble.bubble_id,
                "type": "strain" if kind == "inter_species" else kind,
                "strain_ids": ",".join(sorted({label for label in labels if label}, key=int)),
                "source_id": ";".join(members[bubble.source[0]]),
                "sink_id": ";".join(members[bubble.sink[0]]),
                "branch_ids": ",".join(tokens),
                "expected_action": "retain" if kind == "inter_species" else "",
                "class": kind,
                "branch_genomes": ",".join(genome_of.get(label, "") for label in labels),
                "branch_species": ",".join(labels),
            }
        )
    return rows


def write_phage_bubbles(outdir: Path) -> Path:
    """Write ``ground_truth/bubbles.tsv`` of a built phage bench from its CDBG and gold taxa."""
    from metametro.bench.ground_truth import load_bench_gold
    from metametro.formats.cdbg.io import load_cdbg

    outdir = Path(outdir)
    gold = load_bench_gold(outdir)
    _, values = read_tsv(outdir / "ground_truth" / "taxon_values.tsv")
    genome_of = {row["taxon"]: row["genome"] for row in values}
    rows = bubble_rows(load_cdbg(outdir / "cdbg"), gold.nodes, genome_of)
    path = outdir / "ground_truth" / "bubbles.tsv"
    write_tsv(path, list(BUBBLE_COLUMNS), rows)
    return path
