"""Decaying colour leakage along CFA neighbours.

Each step mixes a node with its neighbours, then restarts on that node's own
evidence. The restart term is what makes leaked mass decay with hop distance.
CFA colours stay integer presence masks. A colour is kept when its mass is at
least ``threshold`` after leakage; otherwise the argmax is kept when it is
strictly positive.

Decaying reads one source layer and never a gold namespace
(``GOLD_NAMESPACES``, i.e. ``taxon``). Without an explicit ``source`` it takes
the first of ``DECAY_SOURCES`` on the graph; when none is there it takes every
namespace except gold, ``decaying``, and ``accession``.
"""

from __future__ import annotations

from typing import Sequence

from metametro.contracts.colouring import GOLD_NAMESPACES, paint_namespace
from metametro.errors import ContractError
from metametro.formats.cfa.model import CfaGraph
from metametro.formats.cfa.validator import parse_color_set

DECAY_SOURCES = ("sample", "kraken2", "kaiju", "composition")
"""Default source layers for decaying, most preferred first."""

_NEVER_SOURCE = frozenset({"decaying", "accession"})


def decay_sources(cfa: CfaGraph, source: str | Sequence[str] | None = None) -> tuple[str, ...]:
    """Return the colour namespaces decaying reads on ``cfa``.

    An explicit ``source`` naming a gold namespace, ``decaying``, or
    ``accession`` raises ``ContractError``. An empty tuple means no allowed
    layer is on the graph.
    """
    present = []
    for row in cfa.colors or []:
        namespace = str(row.get("namespace", ""))
        if namespace not in present:
            present.append(namespace)
    if source is not None:
        wanted = (source,) if isinstance(source, str) else tuple(source)
        bad = [name for name in wanted if name in GOLD_NAMESPACES or name in _NEVER_SOURCE]
        if bad:
            raise ContractError([f"decaying cannot read gold or derived colours: {', '.join(bad)}"])
        return tuple(name for name in wanted if name in present)
    for name in DECAY_SOURCES:
        if name in present:
            return (name,)
    return tuple(name for name in present if name not in GOLD_NAMESPACES and name not in _NEVER_SOURCE)


def normalize(weights: dict[int, float]) -> dict[int, float]:
    """Return a probability vector over positive masses."""
    total = sum(value for value in weights.values() if value > 0)
    if total <= 0:
        return {}
    return {key: value / total for key, value in weights.items() if value > 0}


def decaying_distributions(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    *,
    decay: float = 0.5,
    iterations: int = 4,
) -> dict[str, dict[int, float]]:
    """Spread label mass along neighbours without replacing a node's own evidence.

    ``decay`` must lie in ``[0, 1]`` and ``iterations`` must be at least 1.
    Input maps are not modified.
    """
    if not 0.0 <= decay <= 1.0:
        raise ContractError(["decay must be in [0, 1]"])
    if iterations < 1:
        raise ContractError(["iterations must be >= 1"])
    evidence = {node: normalize(weights) for node, weights in distributions.items()}
    current = {node: dict(weights) for node, weights in evidence.items()}
    neighbours = _neighbours(edges)
    for _ in range(iterations):
        updated: dict[str, dict[int, float]] = {}
        for node, evidence_weights in evidence.items():
            if node not in neighbours:
                updated[node] = dict(evidence_weights)
                continue
            vote = _neighbour_vote(neighbours[node], current)
            if not vote:
                updated[node] = dict(evidence_weights)
                continue
            mixed = {
                taxon_id: (1.0 - decay) * evidence_weights.get(taxon_id, 0.0) + decay * vote.get(taxon_id, 0.0)
                for taxon_id in set(evidence_weights) | set(vote)
            }
            updated[node] = normalize(mixed)
        current = updated
    return current


def _neighbours(edges: list[dict]) -> dict[str, list[tuple[str, float]]]:
    neighbours: dict[str, list[tuple[str, float]]] = {}
    for edge in edges:
        weight = float(edge.get("weight", 1.0))
        if weight <= 0:
            continue
        neighbours.setdefault(edge["source"], []).append((edge["target"], weight))
    return neighbours


def _neighbour_vote(
    links: list[tuple[str, float]],
    current: dict[str, dict[int, float]],
) -> dict[int, float]:
    total = 0.0
    accumulated: dict[int, float] = {}
    for other, weight in links:
        other_weights = current.get(other)
        if not other_weights:
            continue
        total += weight
        for taxon_id, probability in other_weights.items():
            accumulated[taxon_id] = accumulated.get(taxon_id, 0.0) + weight * probability
    if total <= 0:
        return {}
    return {taxon_id: value / total for taxon_id, value in accumulated.items() if value > 0}


def _keep(mass: dict[int, float], lookup: dict[int, str], threshold: float) -> list[str]:
    kept = [lookup[color_id] for color_id, value in mass.items() if value >= threshold and color_id in lookup]
    if kept:
        return kept
    if not mass:
        return []
    color_id, value = max(mass.items(), key=lambda item: item[1])
    if value > 0 and color_id in lookup:
        return [lookup[color_id]]
    return []


def colour_decaying(
    cfa: CfaGraph,
    *,
    decay: float = 0.5,
    iterations: int = 4,
    threshold: float = 0.5,
    operation: str = "merge",
    source: str | Sequence[str] | None = None,
) -> CfaGraph:
    """Paint a ``decaying`` namespace by leaking one allowed colour layer along edges.

    ``source`` picks the layer (see ``decay_sources``). Gold colours are never read.
    """
    namespaces = decay_sources(cfa, source)
    if not namespaces:
        raise ContractError(["decaying needs a non-gold colour layer on the graph"])
    lookup = {
        int(row["color_id"]): str(row.get("value") or row["color_id"])
        for row in cfa.colors or []
        if str(row.get("namespace", "")) in namespaces
    }
    distributions = {
        row["node_id"]: {
            color_id: 1.0 for color_id in parse_color_set(row.get("color_set", "")) if color_id in lookup
        }
        for row in cfa.nodes
    }
    edges = [{"source": row["source"], "target": row["target"], "weight": 1.0} for row in cfa.edges]
    leaked = decaying_distributions(distributions, edges, decay=decay, iterations=iterations)
    node_values = {node_id: _keep(mass, lookup, threshold) for node_id, mass in leaked.items()}
    by_node = {node_id: set(values) for node_id, values in node_values.items()}
    edge_values: dict[str, list[str]] = {}
    for row in cfa.edges:
        shared = sorted(by_node.get(row["source"], set()) & by_node.get(row["target"], set()))
        edge_values[row["edge_id"]] = shared
    return paint_namespace(cfa, node_values, edge_values, namespace="decaying", operation=operation)
