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
    (``sample``, ``transport_type``, ``route``, ...). Gold ``taxon`` is not a
    colour namespace of a built graph; read it with ``load_bench_gold``.
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


def load_bench_cfa(outdir: Path, *, namespaces: Sequence[str] | None = None) -> CfaGraph:
    """Load the CFA from a ``benchbuild`` directory and optionally filter colours."""
    graph = load_cfa(Path(outdir) / "cfa")
    if namespaces:
        graph = filter_colours(graph, namespaces=namespaces)
    return graph


def load_bench_cdbg(outdir: Path, *, namespaces: Sequence[str] | None = None) -> Cdbg:
    """Load the CDBG (ToCUMG) and keep the named colour namespaces."""
    graph = load_cdbg(Path(outdir) / "cdbg")
    if namespaces:
        graph = filter_colours(graph, namespaces=namespaces)
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
