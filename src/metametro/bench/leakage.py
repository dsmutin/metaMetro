"""Refuse evaluation targets that were written into a graph a model can read.

Ground truth belongs under ``ground_truth/``. It must not be a feature column
or a colour namespace on CFA, CDBG, or CGT. Each bench declares its gold
namespaces (``BenchSpec.gold_namespaces``; ``GOLD_NAMESPACES``, i.e. ``taxon``,
for a graph of no known bench). A colour namespace, a column, or a label
dictionary namespace in that set is rejected like ``truth`` or ``target``.
"""

from __future__ import annotations

from collections.abc import Iterable

from metametro.contracts.colouring import GOLD_NAMESPACES
from metametro.errors import ContractError
from metametro.formats.cfa.model import CfaGraph
from metametro.formats.cgt.model import Cgt

_FORBIDDEN = (
    "truth_taxon",
    "truth_taxon_id",
    "expected_action",
    "genome_abundance",
    "target_metric",
    "gold_standard",
    "gold_taxon",
)
_FORBIDDEN_NAMESPACES = frozenset({"truth", "gold", "target", "expected_action"})


def _names(graph: CfaGraph) -> list[str]:
    found: list[str] = []
    features = graph.metadata.get("features") or {}
    if isinstance(features, dict):
        for kind in ("node", "edge"):
            declared = features.get(kind) or {}
            if isinstance(declared, dict):
                found.extend(str(name) for name in declared)
    found.extend(graph.node_header)
    found.extend(graph.edge_header)
    return found


def assert_no_target_leak(
    graph: CfaGraph, *, where: str = "CFA", gold_namespaces: Iterable[str] = GOLD_NAMESPACES
) -> None:
    """Raise ``ContractError`` when a target column, namespace, or label is on the graph.

    ``gold_namespaces`` is the bench's declared gold set.
    """
    gold = frozenset(str(name).lower() for name in gold_namespaces)
    forbidden = _FORBIDDEN_NAMESPACES | gold
    bad = []
    for name in dict.fromkeys(_names(graph)):
        folded = name.lower()
        if any(token in folded for token in _FORBIDDEN) or folded in gold:
            bad.append(f"{where} column {name}")
    for row in graph.colors or []:
        namespace = str(row.get("namespace", "")).lower()
        if namespace in forbidden:
            bad.append(f"{where} colour namespace {namespace}")
    for namespace in sorted({str(row.get("namespace", "")).lower() for row in graph.labels or []} & forbidden):
        bad.append(f"{where} label namespace {namespace}")
    if bad:
        raise ContractError(["target leaked into the graph: " + "; ".join(bad)])


def assert_cgt_no_target_leak(graph: Cgt) -> None:
    """Raise ``ContractError`` when a CGT feature name is an evaluation target."""
    names: list[str] = []
    for key in ("node_feature_names", "edge_feature_names"):
        raw = graph.metadata.get(key) or []
        if isinstance(raw, list):
            names.extend(str(item) for item in raw)
    bad = [name for name in names if any(token in name.lower() for token in _FORBIDDEN)]
    if bad:
        raise ContractError(["target leaked into CGT features: " + ", ".join(bad)])
