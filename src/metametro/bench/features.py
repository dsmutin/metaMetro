"""Feature layer of a benchmark CGT: k-mer composition and coverage.

These columns are model inputs in ``node_features``. They are not colours, so they
add no namespace and do not change ``node_colors``, the colour dictionary, or
the bench identity. Evaluation targets never enter this layer.
"""

from __future__ import annotations

import numpy as np

from metametro.contracts.composition import canonical_kmer_index, kmer_vector
from metametro.errors import ContractError
from metametro.formats.cdbg.model import Cdbg
from metametro.formats.cfa.model import CfaGraph

DEFAULT_K = 4


def _node_coverage(cfa: CfaGraph) -> dict[str, float] | None:
    if not cfa.nodes or "coverage" not in cfa.node_header:
        return None
    values: dict[str, float] = {}
    for row in cfa.nodes:
        try:
            number = float(row["coverage"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ContractError([f"node {row.get('node_id', '')} has a non-numeric coverage"]) from exc
        if not np.isfinite(number):
            raise ContractError([f"node {row['node_id']} has a non-finite coverage"])
        values[row["node_id"]] = number
    return values


def node_feature_layer(cdbg: Cdbg, cfa: CfaGraph, *, k: int = DEFAULT_K) -> tuple[dict[str, np.ndarray], list[str]]:
    """Return per-unitig features and their names.

    ``kmer:k{k}:<kmer>`` columns are L2-normalised canonical k-mer counts of
    the unitig sequence (136 columns for ``k=4``). They are written only when
    every unitig has a sequence. ``coverage`` is the mean CFA node coverage of
    the unitig members, written only when the CFA has a ``coverage`` column.
    A graph with neither gets no columns.
    """
    use_kmers = bool(cdbg.unitigs) and all(unitig.sequence for unitig in cdbg.unitigs)
    index = canonical_kmer_index(k) if use_kmers else {}
    coverage = _node_coverage(cfa)
    names = [f"kmer:k{k}:{kmer}" for kmer in index]
    if coverage is not None:
        names.append("coverage")
    rows: dict[str, np.ndarray] = {}
    for unitig in cdbg.unitigs:
        parts = [kmer_vector(unitig.sequence, k, index)] if use_kmers else []
        if coverage is not None:
            parts.append(np.asarray([np.mean([coverage[member] for member in unitig.members])]))
        rows[unitig.unitig_id] = (
            np.concatenate(parts).astype(np.float32) if parts else np.zeros(0, dtype=np.float32)
        )
    return rows, names
