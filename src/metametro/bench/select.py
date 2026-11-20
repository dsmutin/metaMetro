"""Load a benchbuild ToCUMG and keep the colour namespaces a tool asked for.

Downstream packages do not invent colourings. They name the MetaMetro
namespaces they consume.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from metametro.bench.colourings import colouring_namespace
from metametro.bench.features import DEFAULT_K, node_feature_layer
from metametro.contracts.colour_filter import filter_colours
from metametro.contracts.colouring import GOLD_NAMESPACES
from metametro.converters.cdbg_to_cgt import cdbg_to_cgt
from metametro.errors import ContractError
from metametro.formats.cdbg.io import load_cdbg
from metametro.formats.cdbg.model import Cdbg
from metametro.formats.cfa.io import load_cfa
from metametro.formats.cfa.model import CfaGraph
from metametro.formats.cgt.model import Cgt
from metametro.tables import read_yaml


def namespaces_for(colourings: Sequence[str]) -> list[str]:
    """Map colouring names to colour-dictionary namespaces.

    ``as_built`` has no single namespace and is skipped here. Select those
    colours by naming the namespaces that the builder already wrote
    (``sample``, ``composition``, ...). A bench's gold namespaces (``taxon``;
    ``transport_type`` and ``route`` on spb) are not colour namespaces of a
    built graph; read them with ``load_bench_gold``.
    """
    found: list[str] = []
    missing: list[str] = []
    for name in colourings:
        namespace = colouring_namespace(name)
        if namespace is None:
            if name != "as_built":
                missing.append(str(name))
            continue
        if namespace not in found:
            found.append(namespace)
    if missing:
        raise ContractError([f"unknown colouring: {', '.join(missing)}"])
    return found


def bench_gold_namespaces(outdir: Path) -> frozenset[str]:
    """Gold namespaces declared in ``manifest.yaml``, plus ``GOLD_NAMESPACES``.

    A directory without a manifest, or a manifest written before benches
    declared their gold, gives ``GOLD_NAMESPACES``.
    """
    path = Path(outdir) / "manifest.yaml"
    declared = read_yaml(path).get("gold_namespaces") if path.is_file() else None
    return GOLD_NAMESPACES | frozenset(str(name) for name in declared or [])


def _reject_gold(outdir: Path, namespaces: Sequence[str] | None) -> frozenset[str]:
    """Fail loudly instead of returning an uncoloured graph for a gold namespace."""
    declared = bench_gold_namespaces(outdir)
    gold = sorted({str(name) for name in namespaces or []} & declared)
    if gold:
        raise ContractError(
            [f"{name} is gold ground truth, not a colour namespace; read it with load_bench_gold" for name in gold]
        )
    return declared


def load_bench_cfa(outdir: Path, *, namespaces: Sequence[str] | None = None) -> CfaGraph:
    """Load the CFA from a ``benchbuild`` directory and optionally filter colours."""
    declared = _reject_gold(outdir, namespaces)
    graph = load_cfa(Path(outdir) / "cfa")
    if namespaces:
        graph = filter_colours(graph, namespaces=namespaces, gold_namespaces=declared)
    return graph


def load_bench_cdbg(outdir: Path, *, namespaces: Sequence[str] | None = None) -> Cdbg:
    """Load the CDBG (ToCUMG) and keep the named colour namespaces."""
    declared = _reject_gold(outdir, namespaces)
    graph = load_cdbg(Path(outdir) / "cdbg")
    if namespaces:
        graph = filter_colours(graph, namespaces=namespaces, gold_namespaces=declared)
    return graph


def load_bench_cgt(
    outdir: Path,
    *,
    namespaces: Sequence[str] | None = None,
    features: bool = True,
) -> Cgt:
    """Load a CGT. Namespace selection is applied on the CDBG, then retensorised.

    A stored CGT has no colour dictionary, so a namespace filter cannot be
    applied to the ``.npy`` masks directly. ``features=True`` rebuilds the
    feature layer (k-mer composition and coverage) in ``node_features`` with the
    ``k`` recorded in ``manifest.yaml``. ``features=False`` gives an empty
    ``node_features`` as before. Colours do not depend on this flag.
    """
    cdbg = load_bench_cdbg(outdir, namespaces=namespaces)
    if not features:
        return cdbg_to_cgt(cdbg, node_feature_names=[], edge_feature_names=[])
    manifest_path = Path(outdir) / "manifest.yaml"
    manifest = read_yaml(manifest_path) if manifest_path.is_file() else {}
    k = int((manifest.get("features") or {}).get("kmer_k", DEFAULT_K))
    node_features, node_names = node_feature_layer(cdbg, load_cfa(Path(outdir) / "cfa"), k=k)
    return cdbg_to_cgt(
        cdbg,
        node_features=node_features if node_names else None,
        node_feature_names=node_names,
        edge_feature_names=[],
    )


def load_bench_states(outdir: Path, *, node_namespace: str, edge_namespaces: tuple[str, ...]):
    """``build_geometric_states`` on a bench CFA, with targets read from its gold.

    The node classes (``node_labels``, ``node_class_names``) and edge targets of
    a gold namespace come from ``ground_truth/``; the returned CDBG and CGT carry
    only the bench's colours. Requires longitude and latitude columns (spb, roxel).
    """
    from metametro.bench.ground_truth import load_bench_gold_all
    from metametro.states import build_geometric_states

    return build_geometric_states(
        load_cfa(Path(outdir) / "cfa"),
        node_namespace=node_namespace,
        edge_namespaces=tuple(edge_namespaces),
        gold=load_bench_gold_all(outdir),
    )
