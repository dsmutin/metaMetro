"""Gold ``taxon`` labels of the phage benches, derived from the pinned genomes.

A CFA node belongs to the pinned genome that contains a strict majority of
the node's k-mers (both strands); any other node has no label. Its gold value
is that genome's NCBI species taxid from ``genomes.tsv``. An edge carries the
label its two ends share. ``materialize`` moves the ``taxon`` colours to
``ground_truth/taxon_nodes.tsv`` and ``taxon_edges.tsv``, so they are used only
for scoring and never reach a graph a model reads. ``ground_truth/taxon_values.tsv``
names each value: rank, species name, genome, accession, and record taxid.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from metametro.bench.genomes import pinned_genomes
from metametro.contracts.assembly import load_genomes
from metametro.contracts.colouring import paint_namespace
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import write_tsv

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
VALUE_COLUMNS = ("taxon", "rank", "name", "genome", "accession", "record_taxid")


def node_genomes(sequences: Mapping[str, str], genomes: Mapping[str, str], k: int) -> dict[str, str]:
    """Map node id to the genome owning a strict majority of its k-mers, or ``""``."""
    names = sorted(genomes)
    genome_kmers: dict[str, set[str]] = {}
    for name in names:
        sequence = genomes[name]
        kmers: set[str] = set()
        for strand in (sequence, sequence.translate(_COMPLEMENT)[::-1]):
            if len(strand) >= k:
                kmers.update(strand[offset : offset + k] for offset in range(len(strand) - k + 1))
        genome_kmers[name] = kmers
    owners: dict[str, str] = {}
    for node_id, sequence in sequences.items():
        owners[node_id] = ""
        if len(sequence) < k:
            continue
        kmers = [sequence[offset : offset + k] for offset in range(len(sequence) - k + 1)]
        counts = [sum(kmer in genome_kmers[name] for kmer in kmers) for name in names]
        best = max(counts) if counts else 0
        winners = [position for position, count in enumerate(counts) if count == best and count > 0]
        if len(winners) == 1 and best * 2 > len(kmers):
            owners[node_id] = names[winners[0]]
    return owners


def paint_phage_gold(graph: CfaGraph, bench: str, genome_dir: Path, k: int) -> CfaGraph:
    """Add the gold ``taxon`` namespace (species taxids) to a phage CFA."""
    pinned = {genome.name: genome for genome in pinned_genomes(bench)}
    genomes = {name: sequence for name, sequence in load_genomes(genome_dir).items() if name in pinned}
    owners = node_genomes(graph.sequences, genomes, k)
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
