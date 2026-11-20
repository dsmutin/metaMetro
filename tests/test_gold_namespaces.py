"""Mandatory: each bench declares its gold namespaces, and none of them is readable as a colour."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metametro.bench import (
    bench_gold_namespaces,
    build,
    load_bench_cdbg,
    load_bench_cfa,
    load_bench_gold,
    load_bench_gold_all,
    load_bench_states,
    list_specs,
    resolve,
)
from metametro.bench.colourings import BuildContext
from metametro.bench.leakage import assert_no_target_leak
from metametro.bench.materialize import materialize
from metametro.digraph_cfa import coloured_digraph_cfa
from metametro.errors import ContractError
from metametro.formats.cdbg.io import load_cdbg
from metametro.states import build_geometric_states
from metametro.transit_cfa import TransitRoute, transit_cfa

pytestmark = pytest.mark.mandatory


def _transit():
    stops = {"a": "Alpha", "b": "Beta", "c": "Gamma", "d": "Delta"}
    routes = [
        TransitRoute("bus2", "2", "bus"),
        TransitRoute("tram1", "1", "tram"),
        TransitRoute("bus9", "9", "bus"),
    ]
    trips = [("bus2", ["a", "b", "c"]), ("tram1", ["b", "c", "d"])]
    coordinates = {"a": (30.0, 59.0), "b": (30.1, 59.1), "c": (30.2, 59.2), "d": (30.3, 59.3)}
    return transit_cfa(stops, routes, trips, graph_id="toy_transit", k=5, coordinates=coordinates)


def _gold_columns(graph, gold) -> set[str]:
    names = {row["namespace"] for row in graph.colors or []}
    names |= {row["namespace"] for row in getattr(graph, "labels", None) or []}
    return names & set(gold)


def test_declared_gold_sets() -> None:
    assert resolve("spb_ground_transit").gold_namespaces == ("transport_type", "route")
    assert resolve("roxel").gold_namespaces == ("type", "name")
    assert resolve("phage_10").gold_namespaces == ("taxon",)
    assert resolve("bubble_strain_2").gold_namespaces == ("taxon",)
    assert resolve("mock_bubble").gold_namespaces == ("taxon", "genome")
    assert resolve("bubble_reads_2").gold_namespaces == ("taxon", "genome", "accession")


@pytest.mark.parametrize("name", [spec.name for spec in list_specs() if spec.kind == "inprocess"])
def test_no_gold_namespace_is_a_colour_on_any_offline_bench(name: str, tmp_path: Path) -> None:
    result = build(name, outdir=tmp_path / name)
    gold = resolve(name).gold_namespaces
    assert set(gold) <= bench_gold_namespaces(result.outdir)
    cfa = load_bench_cfa(result.outdir)
    assert_no_target_leak(cfa, gold_namespaces=gold)
    assert not _gold_columns(cfa, gold)
    assert not set(gold) & (set(cfa.node_header) | set(cfa.edge_header) | {"label"})
    assert not _gold_columns(load_cdbg(result.outdir / "cdbg"), gold)
    for namespace in gold:
        with pytest.raises(ContractError, match="load_bench_gold"):
            load_bench_cdbg(result.outdir, namespaces=[namespace])


def test_mock_label_column_moves_to_ground_truth(tmp_path: Path) -> None:
    result = build("mock_bubble", outdir=tmp_path / "mock")
    cfa = load_bench_cfa(result.outdir)
    assert "label" not in cfa.node_header and cfa.labels is None
    genome = load_bench_gold(result.outdir, "genome")
    assert genome.nodes["n000003"] == ("genome_002",)
    assert set(load_bench_gold_all(result.outdir)) == {"genome"}


def test_transit_gold_leaves_the_colours_and_keeps_the_state_labels(tmp_path: Path) -> None:
    graph = _transit()
    targets = ("transport_type", "route")
    before = build_geometric_states(graph, node_namespace="transport_type", edge_namespaces=targets)
    materialize(graph, tmp_path, BuildContext(), gold_namespaces=("transport_type", "route"))
    cfa = load_bench_cfa(tmp_path)
    assert not {row["namespace"] for row in cfa.colors or []} & {"transport_type", "route"}
    assert "coverage" in cfa.node_header
    gold = load_bench_gold_all(tmp_path)
    assert set(gold) == {"transport_type", "route"}
    assert "bus:9" in gold["route"].values
    after = load_bench_states(tmp_path, node_namespace="transport_type", edge_namespaces=targets)
    assert after.cgt.node_labels.size == before.cgt.node_labels.size > 0
    assert np.array_equal(after.cgt.node_labels, before.cgt.node_labels)
    assert after.node_class_names == before.node_class_names
    assert after.edge_label_names == before.edge_label_names
    assert np.array_equal(after.edge_targets, before.edge_targets)
    assert all(row["namespace"] not in {"transport_type", "route"} for row in after.cdbg.colors or [])
    with pytest.raises(ContractError, match="load_bench_gold"):
        load_bench_cfa(tmp_path, namespaces=["route"])


def test_digraph_attribute_columns_and_colours_leave_the_graph(tmp_path: Path) -> None:
    nodes = [
        {"node_id": "n1", "longitude": "0.1", "latitude": "1.2"},
        {"node_id": "n2", "longitude": "0.2", "latitude": "1.3"},
        {"node_id": "n3", "longitude": "0.3", "latitude": "1.4"},
    ]
    edges = [
        {"source": "n1", "target": "n2", "type": "residential", "name": "Roxelweg"},
        {"source": "n1", "target": "n3", "type": "path", "name": "Pfad"},
    ]
    graph = coloured_digraph_cfa(
        nodes,
        edges,
        graph_id="toy_roxel",
        k=5,
        colour_columns=("type", "name"),
        keep_columns=("longitude", "latitude"),
    )
    before = build_geometric_states(graph, node_namespace="type", edge_namespaces=("name",))
    materialize(graph, tmp_path, BuildContext(), gold_namespaces=("type", "name"))
    cfa = load_bench_cfa(tmp_path)
    assert "type" not in cfa.edge_header and "name" not in cfa.edge_header
    assert load_bench_gold(tmp_path, "name").edges == {"e000001": ("Roxelweg",), "e000002": ("Pfad",)}
    states = load_bench_states(tmp_path, node_namespace="type", edge_namespaces=("name",))
    assert states.cgt.node_labels.size == 3
    assert np.array_equal(states.cgt.node_labels, before.cgt.node_labels)
    assert states.node_class_names == before.node_class_names


def test_guard_rejects_a_declared_gold_column_or_label() -> None:
    graph = _transit()
    with pytest.raises(ContractError, match="colour namespace route"):
        assert_no_target_leak(graph, gold_namespaces=("transport_type", "route"))
    assert_no_target_leak(graph)


def test_a_bench_built_before_gold_namespaces_must_be_rebuilt(tmp_path: Path) -> None:
    from metametro.bench import bench_gold_namespaces, build, load_bench_cdbg, load_bench_gold_all
    from metametro.errors import ContractError
    from metametro.tables import read_yaml, write_yaml

    build("mock_bubble", outdir=tmp_path / "mock")
    manifest = read_yaml(tmp_path / "mock" / "manifest.yaml")
    assert "genome" in bench_gold_namespaces(tmp_path / "mock")
    del manifest["gold_namespaces"]
    write_yaml(tmp_path / "mock" / "manifest.yaml", manifest)
    for load in (bench_gold_namespaces, load_bench_cdbg, load_bench_gold_all):
        with pytest.raises(ContractError, match="rebuild the bench"):
            load(tmp_path / "mock")
    assert bench_gold_namespaces(tmp_path / "no_manifest") == frozenset({"taxon"})
