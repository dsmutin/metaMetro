"""Gold labels of a bench, kept out of the graphs a model reads.

A builder may hand ``materialize`` a CFA whose colours include a gold
namespace (``GOLD_NAMESPACES``, i.e. ``taxon``). Those colours are removed from
the graph before any colouring runs and written under ``ground_truth/``:

- ``ground_truth/<namespace>_nodes.tsv``: ``node_id``, ``<namespace>``
- ``ground_truth/<namespace>_edges.tsv``: ``edge_id``, ``source``, ``target``, ``<namespace>``

One row per CFA node or edge. The label cell is a comma-separated list of gold
values, empty when the node or edge has none. Evaluation reads them with
``load_bench_gold`` (per CFA id) or ``load_bench_gold_cgt`` (aligned with the
nodes and CSR edges of ``load_bench_cgt``).
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from metametro.contracts.colouring import GOLD_NAMESPACES, paint_namespace
from metametro.converters.cdbg_to_cgt import cdbg_to_cgt
from metametro.converters.cfa_to_cdbg import cfa_to_cdbg
from metametro.errors import ContractError
from metametro.formats.cfa.io import load_cfa
from metametro.formats.cfa.model import CfaGraph
from metametro.formats.cfa.validator import parse_color_set
from metametro.formats.cgt.model import Cgt
from metametro.tables import read_tsv, write_tsv


@dataclass(frozen=True)
class GoldLabels:
    """Gold values of one namespace by CFA node id and CFA edge id."""

    namespace: str
    values: tuple[str, ...]
    nodes: dict[str, tuple[str, ...]]
    edges: dict[str, tuple[str, ...]]


def split_gold(graph: CfaGraph) -> tuple[CfaGraph, dict[str, GoldLabels]]:
    """Return ``graph`` without gold colours, and the gold labels by namespace."""
    gold_ids: dict[int, tuple[str, str]] = {}
    for row in graph.colors or []:
        namespace = str(row.get("namespace", ""))
        if namespace in GOLD_NAMESPACES:
            gold_ids[int(row["color_id"])] = (namespace, str(row.get("value") or row["color_id"]))
    if not gold_ids:
        return graph, {}
    stripped = deepcopy(graph)
    namespaces = sorted({namespace for namespace, _ in gold_ids.values()})
    nodes: dict[str, dict[str, tuple[str, ...]]] = {name: {} for name in namespaces}
    edges: dict[str, dict[str, tuple[str, ...]]] = {name: {} for name in namespaces}
    for rows, key, out in ((stripped.nodes, "node_id", nodes), (stripped.edges, "edge_id", edges)):
        for row in rows:
            column = "color_set" if "color_set" in row else ("colors" if "colors" in row else None)
            ids = parse_color_set(row.get(column, "")) if column else []
            for name in namespaces:
                out[name][row[key]] = tuple(gold_ids[i][1] for i in sorted(ids) if i in gold_ids and gold_ids[i][0] == name)
            if column:
                row[column] = ",".join(str(i) for i in ids if i not in gold_ids)
    stripped.colors = [row for row in graph.colors or [] if int(row["color_id"]) not in gold_ids]
    gold = {
        name: GoldLabels(
            namespace=name,
            values=tuple(sorted({value for namespace, value in gold_ids.values() if namespace == name})),
            nodes=nodes[name],
            edges=edges[name],
        )
        for name in namespaces
    }
    return stripped, gold


def write_gold(outdir: Path, gold: dict[str, GoldLabels], graph: CfaGraph) -> list[Path]:
    """Write ``ground_truth/<namespace>_nodes.tsv`` and ``_edges.tsv``. Return the paths."""
    written: list[Path] = []
    root = Path(outdir) / "ground_truth"
    for name, labels in sorted(gold.items()):
        node_path = root / f"{name}_nodes.tsv"
        write_tsv(
            node_path,
            ["node_id", name],
            [{"node_id": row["node_id"], name: ",".join(labels.nodes.get(row["node_id"], ()))} for row in graph.nodes],
        )
        edge_path = root / f"{name}_edges.tsv"
        write_tsv(
            edge_path,
            ["edge_id", "source", "target", name],
            [
                {
                    "edge_id": row["edge_id"],
                    "source": row["source"],
                    "target": row["target"],
                    name: ",".join(labels.edges.get(row["edge_id"], ())),
                }
                for row in graph.edges
            ],
        )
        written.extend([node_path, edge_path])
    return written


def gold_identity_lines(gold: dict[str, GoldLabels]) -> list[str]:
    """Lines that put gold labels into a build identity."""
    lines: list[str] = []
    for name, labels in sorted(gold.items()):
        for node_id, values in labels.nodes.items():
            lines.append(f"gold\t{name}\tnode\t{node_id}\t{','.join(values)}")
        for edge_id, values in labels.edges.items():
            lines.append(f"gold\t{name}\tedge\t{edge_id}\t{','.join(values)}")
    return lines


def _split(cell: str) -> tuple[str, ...]:
    return tuple(part for part in cell.split(",") if part)


def load_bench_gold(outdir: Path, namespace: str = "taxon") -> GoldLabels:
    """Read the gold labels of ``namespace`` from ``<outdir>/ground_truth/``."""
    root = Path(outdir) / "ground_truth"
    node_path = root / f"{namespace}_nodes.tsv"
    edge_path = root / f"{namespace}_edges.tsv"
    if not node_path.is_file() or not edge_path.is_file():
        raise ContractError([f"no gold {namespace} labels under {root}"])
    _, node_rows = read_tsv(node_path)
    _, edge_rows = read_tsv(edge_path)
    nodes = {row["node_id"]: _split(row.get(namespace, "")) for row in node_rows}
    edges = {row["edge_id"]: _split(row.get(namespace, "")) for row in edge_rows}
    values = sorted({value for found in (*nodes.values(), *edges.values()) for value in found})
    return GoldLabels(namespace=namespace, values=tuple(values), nodes=nodes, edges=edges)


def load_bench_gold_cgt(outdir: Path, namespace: str = "taxon") -> Cgt:
    """Return a CGT of the bench whose only colours are the gold labels.

    Nodes and CSR edges are in the order of ``load_bench_cgt(outdir)``, so
    ``node_colors[:, j]`` is the multi-hot gold mask of ``GoldLabels.values[j]``
    for each CGT node (union over the unitig's CFA members). It has no node
    features. Use it for evaluation only, never as a model input.
    """
    labels = load_bench_gold(outdir, namespace)
    graph = load_cfa(Path(outdir) / "cfa")
    painted = paint_namespace(
        graph,
        {row["node_id"]: list(labels.nodes.get(row["node_id"], ())) for row in graph.nodes},
        {row["edge_id"]: list(labels.edges.get(row["edge_id"], ())) for row in graph.edges},
        namespace=namespace,
        operation="replace",
        values=list(labels.values),
    )
    return cdbg_to_cgt(cfa_to_cdbg(painted), node_feature_names=[], edge_feature_names=[])
