"""Gold labels of a bench, kept out of the graphs a model reads.

A builder may hand ``materialize`` a CFA that carries a gold namespace of its
bench (``BenchSpec.gold_namespaces``; ``GOLD_NAMESPACES``, i.e. ``taxon``, by
default) as colours, as a column, or as the ``label`` column. That gold is
removed from the graph before any colouring runs and written under
``ground_truth/``:

- ``ground_truth/<namespace>_nodes.tsv``: ``node_id``, ``<namespace>``
- ``ground_truth/<namespace>_edges.tsv``: ``edge_id``, ``source``, ``target``, ``<namespace>``
- ``ground_truth/<namespace>_values.tsv``: ``<namespace>``, every value of the namespace
  in colour-dictionary order, also one that labels no node or edge (it may have more
  columns, as on the phage benches)

One row per CFA node or edge. The label cell is a comma-separated list of gold
values, empty when the node or edge has none. Evaluation reads them with
``load_bench_gold`` (per CFA id) or ``load_bench_gold_cgt`` (aligned with the
nodes and CSR edges of ``load_bench_cgt``).
"""

from __future__ import annotations

from copy import deepcopy
from collections.abc import Iterable
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


def split_gold(
    graph: CfaGraph, gold_namespaces: Iterable[str] = GOLD_NAMESPACES
) -> tuple[CfaGraph, dict[str, GoldLabels]]:
    """Return ``graph`` without gold, and the gold labels by namespace.

    ``gold_namespaces`` is the bench's declared gold set (default
    ``GOLD_NAMESPACES``). Three forms of gold leave the graph: colours in a gold
    namespace, a node or edge column named after a gold namespace, and the
    ``label`` column when its label dictionary namespace is gold. Colours and
    labels become ``GoldLabels``; a named column is dropped (its values are the
    colours of the same namespace).
    """
    declared = frozenset(str(name) for name in gold_namespaces)
    gold_ids: dict[int, tuple[str, str]] = {}
    for row in graph.colors or []:
        namespace = str(row.get("namespace", ""))
        if namespace in declared:
            gold_ids[int(row["color_id"])] = (namespace, str(row.get("value") or row["color_id"]))
    label_of: dict[str, tuple[str, str]] = {}
    for row in graph.labels or []:
        namespace = str(row.get("namespace", ""))
        if namespace in declared:
            label_of[str(row["label_id"])] = (namespace, str(row.get("value") or row["label_id"]))
    columns = declared & (set(graph.node_header) | set(graph.edge_header))
    if not gold_ids and not label_of and not columns:
        return graph, {}
    stripped = deepcopy(graph)
    namespaces = sorted({namespace for namespace, _ in (*gold_ids.values(), *label_of.values())})
    nodes: dict[str, dict[str, tuple[str, ...]]] = {name: {} for name in namespaces}
    edges: dict[str, dict[str, tuple[str, ...]]] = {name: {} for name in namespaces}
    drop_label = bool(label_of)
    for rows, key, out in ((stripped.nodes, "node_id", nodes), (stripped.edges, "edge_id", edges)):
        for row in rows:
            column = "color_set" if "color_set" in row else ("colors" if "colors" in row else None)
            ids = parse_color_set(row.get(column, "")) if column else []
            label = label_of.get(str(row.get("label", "")).strip()) if drop_label else None
            for name in namespaces:
                found = [gold_ids[i][1] for i in sorted(ids) if i in gold_ids and gold_ids[i][0] == name]
                if label is not None and label[0] == name and label[1] not in found:
                    found.append(label[1])
                out[name][row[key]] = tuple(found)
            if column:
                row[column] = ",".join(str(i) for i in ids if i not in gold_ids)
            for name in columns:
                row.pop(name, None)
            if drop_label:
                row.pop("label", None)
    dropped = set(columns) | ({"label"} if drop_label else set())
    stripped.node_header = [name for name in stripped.node_header if name not in dropped]
    stripped.edge_header = [name for name in stripped.edge_header if name not in dropped]
    features = stripped.metadata.get("features")
    if isinstance(features, dict):
        for kind in ("node", "edge"):
            declared_kind = features.get(kind)
            if isinstance(declared_kind, dict):
                features[kind] = {name: value for name, value in declared_kind.items() if name not in dropped}
    stripped.colors = [row for row in graph.colors or [] if int(row["color_id"]) not in gold_ids]
    if drop_label:
        stripped.labels = None
    gold = {
        name: GoldLabels(
            namespace=name,
            values=tuple(
                dict.fromkeys(
                    value for namespace, value in (*gold_ids.values(), *label_of.values()) if namespace == name
                )
            ),
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
        value_path = root / f"{name}_values.tsv"
        write_tsv(value_path, [name], [{name: value} for value in labels.values])
        written.extend([node_path, edge_path, value_path])
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
    value_path = root / f"{namespace}_values.tsv"
    if value_path.is_file():
        _, value_rows = read_tsv(value_path)
        listed = list(dict.fromkeys(row[namespace] for row in value_rows if row.get(namespace)))
        values = listed + [value for value in values if value not in set(listed)]
    return GoldLabels(namespace=namespace, values=tuple(values), nodes=nodes, edges=edges)


def load_bench_gold_all(outdir: Path) -> dict[str, GoldLabels]:
    """Every gold namespace of a bench that has labels under ``ground_truth/``.

    The namespaces are those declared in ``manifest.yaml`` (``gold_namespaces``),
    or ``taxon`` for an older manifest. A declared namespace without files is
    left out, so a bench with no gold gives an empty mapping.
    """
    from metametro.bench.select import bench_gold_namespaces

    found: dict[str, GoldLabels] = {}
    root = Path(outdir) / "ground_truth"
    for name in sorted(bench_gold_namespaces(outdir)):
        if (root / f"{name}_nodes.tsv").is_file() and (root / f"{name}_edges.tsv").is_file():
            found[name] = load_bench_gold(outdir, name)
    return found


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


@dataclass(frozen=True)
class GoldBubble:
    """One row of ``ground_truth/bubbles.tsv``.

    ``branches`` holds the CFA node ids of each branch. ``bubble_class`` and
    ``branch_species`` are empty on benches whose file has no such columns.
    """

    bubble_id: str
    type: str
    strain_ids: tuple[str, ...]
    source_id: str
    sink_id: str
    branches: tuple[tuple[str, ...], ...]
    expected_action: str
    bubble_class: str
    branch_genomes: tuple[str, ...]
    branch_species: tuple[str, ...]


def load_bench_gold_bubbles(outdir: Path) -> tuple[GoldBubble, ...]:
    """Read ``<outdir>/ground_truth/bubbles.tsv``. Raise ``ContractError`` when it is missing."""
    path = Path(outdir) / "ground_truth" / "bubbles.tsv"
    if not path.is_file():
        raise ContractError([f"no gold bubbles under {path.parent}"])
    _, rows = read_tsv(path)
    found = []
    for row in rows:
        branches = tuple(tuple(part for part in token.split(";") if part) for token in row["branch_ids"].split(","))
        found.append(
            GoldBubble(
                bubble_id=row["bubble_id"],
                type=row.get("type", ""),
                strain_ids=_split(row.get("strain_ids", "")),
                source_id=row.get("source_id", ""),
                sink_id=row.get("sink_id", ""),
                branches=branches,
                expected_action=row.get("expected_action", ""),
                bubble_class=row.get("class", ""),
                branch_genomes=tuple(row["branch_genomes"].split(",")) if "branch_genomes" in row else (),
                branch_species=tuple(row["branch_species"].split(",")) if "branch_species" in row else (),
            )
        )
    return tuple(found)


def gold_digest_lines(outdir: Path) -> list[str]:
    """``sha256sum`` lines for every file under ``outdir/ground_truth/``, sorted by path."""
    import hashlib

    outdir = Path(outdir)
    root = outdir / "ground_truth"
    files = sorted(path for path in root.rglob("*") if path.is_file()) if root.is_dir() else []
    return [f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(outdir).as_posix()}" for path in files]


def gold_digest(outdir: Path) -> str:
    """SHA-256 of the ``gold.sha256`` text for the current ``ground_truth/`` files."""
    import hashlib

    text = "".join(f"{line}\n" for line in gold_digest_lines(outdir))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_gold_digest(outdir: Path) -> str:
    """Write ``gold.sha256`` (one ``sha256sum`` line per ground_truth file) and return ``gold_digest``.

    The gold digest is separate from ``identity.sha256``: changing a gold file
    changes it and leaves the identity as it is.
    """
    outdir = Path(outdir)
    (outdir / "gold.sha256").write_text("".join(f"{line}\n" for line in gold_digest_lines(outdir)), encoding="utf-8")
    return gold_digest(outdir)
