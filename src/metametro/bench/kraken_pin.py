"""Kraken2 outputs pinned in the repository, so a build needs neither kraken2 nor its index.

``data/pins/<bench>/kraken2/nodes.kraken`` is Kraken2's per-sequence output for
the CFA node FASTA of that bench (headers are node ids), and ``kraken2.yaml``
records the Kraken2 version, the index name, date and SHA-256, the command,
and the SHA-256 of the input FASTA and of the output. ``taxdump/`` holds a
trimmed NCBI ``nodes.dmp`` and ``names.dmp`` with the lineage to root of every
taxid in the output and in the gold labels; ``pinned_taxdump(name)`` returns it,
so a resolver can read the taxa without the index. A build that was given
no ``kraken_output`` or ``kraken_db`` paints the ``kraken2`` layer from the pin
after checking that its own node FASTA has the pinned SHA-256.
"""

from __future__ import annotations

from pathlib import Path

from metametro.bench.genomes import PINS, sha256_file
from metametro.contracts.classifiers import write_node_fasta
from metametro.errors import ContractError
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import read_yaml


def kraken_pin(bench: str, *, pins: Path | None = None) -> tuple[Path, dict] | None:
    """Return ``(nodes.kraken, metadata)`` pinned for ``bench``, or ``None``."""
    from metametro.bench.registry import resolve

    root = (PINS if pins is None else Path(pins)) / resolve(bench).name / "kraken2"
    output = root / "nodes.kraken"
    if not output.is_file():
        return None
    meta = read_yaml(root / "kraken2.yaml")
    if sha256_file(output) != meta["output_sha256"]:
        raise ContractError([f"{output} does not have the SHA-256 pinned in {root / 'kraken2.yaml'}"])
    taxdump = meta.get("taxdump")
    if taxdump:
        folder = root / taxdump["path"]
        for name in ("nodes", "names"):
            path = folder / f"{name}.dmp"
            if not path.is_file() or sha256_file(path) != taxdump[f"{name}_sha256"]:
                raise ContractError([f"{path} is missing or does not have the SHA-256 pinned in kraken2.yaml"])
        taxdump["dir"] = folder
    return output, meta


def pinned_taxdump(bench: str, *, pins: Path | None = None) -> Path | None:
    """Return the trimmed taxdump directory (``nodes.dmp``, ``names.dmp``) pinned for ``bench``, or ``None``."""
    found = kraken_pin(bench, pins=pins)
    if found is None or not found[1].get("taxdump"):
        return None
    return found[1]["taxdump"]["dir"]


def checked_kraken_pin(bench: str, graph: CfaGraph, fasta: Path, *, pins: Path | None = None) -> Path | None:
    """Write ``graph``'s node FASTA to ``fasta`` and return the pinned output for it.

    Returns ``None`` when ``bench`` has no pin. Raises ``ContractError`` when the
    node FASTA differs from the one the pinned output was made from.
    """
    found = kraken_pin(bench, pins=pins)
    if found is None:
        return None
    output, meta = found
    digest = sha256_file(write_node_fasta(graph, Path(fasta)))
    if digest != meta["input_sha256"]:
        raise ContractError(
            [
                f"the {meta['bench']} node FASTA has SHA-256 {digest}, but the pinned Kraken2 output "
                f"{output} was made from {meta['input_sha256']}; pass --kraken-db or --kraken-output "
                "to classify this graph"
            ]
        )
    return output
