"""Write CFA, CDBG, and CGT after every applicable colouring.

The CGT carries the feature layer from ``metametro.bench.features``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from metametro.bench.colourings import BuildContext, applicable
from metametro.bench.features import node_feature_layer
from metametro.bench.identity import graph_identity, write_identity
from metametro.bench.leakage import assert_cgt_no_target_leak, assert_no_target_leak
from metametro.converters.cdbg_to_cgt import cdbg_to_cgt
from metametro.converters.cfa_to_cdbg import cfa_to_cdbg
from metametro.formats.cdbg.io import dump_cdbg
from metametro.formats.cfa.io import dump_cfa
from metametro.formats.cgt.io import dump_cgt
from metametro.identity import assert_cgt_matches_cdbg
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import write_yaml


@dataclass(frozen=True)
class Materialized:
    """Paths and the identity of one coloured build."""

    outdir: Path
    colourings: tuple[str, ...]
    identity: str
    n_nodes: int
    n_edges: int
    n_unitigs: int


def materialize(graph: CfaGraph, outdir: Path, ctx: BuildContext) -> Materialized:
    """Apply colourings, reject target leakage, and write the three graphs."""
    coloured = graph
    names: list[str] = []
    for colouring in applicable(graph, ctx):
        if ctx.selected is not None and not colouring.available(coloured, ctx):
            from metametro.errors import ContractError

            raise ContractError([f"selected colouring cannot run on this graph: {colouring.name}"])
        coloured = colouring.apply(coloured, ctx)
        names.append(colouring.name)
    if not names:
        from metametro.errors import ContractError

        raise ContractError(["no colouring could be applied"])
    assert_no_target_leak(coloured)
    outdir.mkdir(parents=True, exist_ok=True)
    dump_cfa(coloured, outdir / "cfa")
    cdbg = cfa_to_cdbg(coloured)
    dump_cdbg(cdbg, outdir / "cdbg")
    node_features, node_names = node_feature_layer(cdbg, coloured, k=ctx.composition_k)
    cgt = cdbg_to_cgt(
        cdbg,
        node_features=node_features if node_names else None,
        node_feature_names=node_names,
        edge_feature_names=[],
    )
    assert_cgt_no_target_leak(cgt)
    assert_cgt_matches_cdbg(cgt, cdbg)
    dump_cgt(cgt, outdir / "cgt")
    applied = tuple(names)
    digest = graph_identity(coloured, colourings=applied)
    write_identity(outdir / "identity.sha256", digest)
    write_yaml(
        outdir / "manifest.yaml",
        {
            "status": "built",
            "colourings": list(applied),
            "structures": ["cfa", "cdbg", "cgt"],
            "identity": digest,
            "n_nodes": len(coloured.nodes),
            "n_edges": len(coloured.edges),
            "n_unitigs": len(cdbg.unitigs),
            "features": {"kmer_k": ctx.composition_k, "node": len(node_names)},
        },
    )
    return Materialized(
        outdir=outdir,
        colourings=applied,
        identity=digest,
        n_nodes=len(coloured.nodes),
        n_edges=len(coloured.edges),
        n_unitigs=len(cdbg.unitigs),
    )
