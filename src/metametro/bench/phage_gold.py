"""Gold ``taxon`` labels of the phage benches, derived from the pinned genomes.

A CFA node belongs to a pinned genome when exactly one genome has the largest
count of the node's k-mer positions found in it (on either strand) and that
count is more than half of the node's ``len - k + 1`` k-mer positions; a tie, a
count of half or less, or a node shorter than ``k`` leaves it without a label. Its gold value
is that genome's NCBI species taxid from ``genomes.tsv``. An edge carries the
label its two ends share. ``materialize`` moves the ``taxon`` colours to
``ground_truth/taxon_nodes.tsv`` and ``taxon_edges.tsv``, so they are used only
for scoring and never reach a graph a model reads. ``ground_truth/taxon_values.tsv``
names each value: rank, species name, genome, accession, and record taxid.

``ground_truth/taxon_membership.tsv`` is the origin matrix the single labels are
derived from: per node and pinned genome, whether any of the node's k-mers occur
in that genome, and the fraction of its k-mer positions that do. No alignment is
used; it is the same k-mer index as the majority rule.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from metametro.bench.genomes import pinned_genomes
from metametro.contracts.assembly import load_genomes
from metametro.contracts.colouring import paint_namespace
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import read_tsv, write_tsv

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
VALUE_COLUMNS = ("taxon", "rank", "name", "genome", "accession", "record_taxid")


@dataclass(frozen=True)
class Membership:
    """Which pinned genomes each CFA node comes from.

    ``genomes`` is the genome order (sorted names). ``counts[node][i]`` is the
    number of the node's k-mer positions found in genome ``i`` (either strand),
    out of ``positions[node]`` = ``len - k + 1`` (0 for a node shorter than k).
    """

    genomes: tuple[str, ...]
    positions: dict[str, int]
    counts: dict[str, tuple[int, ...]]

    def members(self, node_id: str) -> tuple[str, ...]:
        """Genomes that hold at least one k-mer of ``node_id``, in genome order."""
        counts = self.counts.get(node_id, ())
        return tuple(name for name, count in zip(self.genomes, counts) if count > 0)

    def owner(self, node_id: str) -> str:
        """The single genome with the most k-mer positions, if more than half, else ``""``."""
        counts = self.counts.get(node_id, ())
        best = max(counts, default=0)
        winners = [position for position, count in enumerate(counts) if count == best and count > 0]
        if len(winners) == 1 and best * 2 > self.positions.get(node_id, 0):
            return self.genomes[winners[0]]
        return ""


def node_membership(sequences: Mapping[str, str], genomes: Mapping[str, str], k: int) -> Membership:
    """Count, per node and genome, the node's k-mer positions found in the genome."""
    names = tuple(sorted(genomes))
    genome_kmers: dict[str, set[str]] = {}
    for name in names:
        sequence = genomes[name]
        kmers: set[str] = set()
        for strand in (sequence, sequence.translate(_COMPLEMENT)[::-1]):
            if len(strand) >= k:
                kmers.update(strand[offset : offset + k] for offset in range(len(strand) - k + 1))
        genome_kmers[name] = kmers
    positions: dict[str, int] = {}
    counts: dict[str, tuple[int, ...]] = {}
    for node_id in sorted(sequences):
        sequence = sequences[node_id]
        kmers = [sequence[offset : offset + k] for offset in range(len(sequence) - k + 1)] if len(sequence) >= k else []
        positions[node_id] = len(kmers)
        counts[node_id] = tuple(sum(kmer in genome_kmers[name] for kmer in kmers) for name in names)
    return Membership(genomes=names, positions=positions, counts=counts)


def node_genomes(sequences: Mapping[str, str], genomes: Mapping[str, str], k: int) -> dict[str, str]:
    """Map node id to the single genome with most of its k-mer positions, if more than half, else ``""``."""
    membership = node_membership(sequences, genomes, k)
    return {node_id: membership.owner(node_id) for node_id in sequences}


def phage_membership(graph: CfaGraph, bench: str, genome_dir: Path, k: int) -> Membership:
    """The genome membership of every node of a phage CFA, from the pinned genomes."""
    pinned = {genome.name for genome in pinned_genomes(bench)}
    genomes = {name: sequence for name, sequence in load_genomes(genome_dir).items() if name in pinned}
    return node_membership(graph.sequences, genomes, k)


def paint_phage_gold(
    graph: CfaGraph, bench: str, genome_dir: Path, k: int, *, membership: Membership | None = None
) -> CfaGraph:
    """Add the gold ``taxon`` namespace (species taxids) to a phage CFA.

    The single label of a node is its strict-majority genome (``Membership.owner``).
    """
    pinned = {genome.name: genome for genome in pinned_genomes(bench)}
    if membership is None:
        membership = phage_membership(graph, bench, genome_dir, k)
    owners = {node_id: membership.owner(node_id) for node_id in graph.sequences}
    node_values = {
        row["node_id"]: [str(pinned[owners[row["node_id"]]].species_taxid)] if owners.get(row["node_id"]) else []
        for row in graph.nodes
    }
    edge_values = {
        row["edge_id"]: sorted(set(node_values.get(row["source"], [])) & set(node_values.get(row["target"], [])))
        for row in graph.edges
    }
    values = sorted({str(genome.species_taxid) for genome in pinned.values()}, key=int)
    return paint_namespace(graph, node_values, edge_values, namespace="taxon", operation="merge", values=values)


def write_taxon_membership(outdir: Path, membership: Membership) -> Path:
    """Write ``ground_truth/taxon_membership.tsv``.

    Columns: ``node_id``, ``kmers`` (the node's k-mer positions), then per genome
    ``<genome>`` (1 when the node has a k-mer of it, else 0) and
    ``<genome>_fraction`` (found positions / ``kmers``, six decimals).
    """
    header = ["node_id", "kmers"]
    for name in membership.genomes:
        header.extend((name, f"{name}_fraction"))
    rows = []
    for node_id in sorted(membership.counts):
        total = membership.positions[node_id]
        row = {"node_id": node_id, "kmers": str(total)}
        for name, count in zip(membership.genomes, membership.counts[node_id]):
            row[name] = "1" if count > 0 else "0"
            row[f"{name}_fraction"] = f"{count / total:.6f}" if total else "0.000000"
        rows.append(row)
    path = Path(outdir) / "ground_truth" / "taxon_membership.tsv"
    write_tsv(path, header, rows)
    return path


def read_taxon_membership(outdir: Path) -> dict[str, tuple[str, ...]] | None:
    """Map node id to its member genomes from ``taxon_membership.tsv``, or ``None`` without the file."""
    path = Path(outdir) / "ground_truth" / "taxon_membership.tsv"
    if not path.is_file():
        return None
    header, rows = read_tsv(path)
    genomes = [column for column in header[2:] if not column.endswith("_fraction")]
    return {row["node_id"]: tuple(name for name in genomes if row[name] == "1") for row in rows}


def write_taxon_values(outdir: Path, bench: str) -> Path:
    """Write ``ground_truth/taxon_values.tsv``: what each gold species taxid names."""
    path = Path(outdir) / "ground_truth" / "taxon_values.tsv"
    rows = [
        {
            "taxon": str(genome.species_taxid),
            "rank": "species",
            "name": genome.species,
            "genome": genome.name,
            "accession": genome.accession,
            "record_taxid": str(genome.taxid),
        }
        for genome in sorted(pinned_genomes(bench), key=lambda item: item.species_taxid)
    ]
    write_tsv(path, list(VALUE_COLUMNS), rows)
    return path
