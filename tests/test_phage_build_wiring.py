"""Mandatory: what the phage branch of ``build`` hands to its steps, with every tool replaced. No network."""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from metametro.bench.build import build
from metametro.formats.cfa.io import dump_cfa
from metametro.formats.cfa.model import CfaGraph

pytestmark = pytest.mark.mandatory

build_module = importlib.import_module("metametro.bench.build")


class Stop(Exception):
    """Raised by the fake ``materialize`` to end the build after the wiring under test."""


def _graph() -> CfaGraph:
    return CfaGraph(
        metadata={"schema_version": "1.0", "graph_id": "toy_phage", "graph_type": "de_bruijn", "k": 3},
        sequences={"n1": "ACGTACGT", "n2": "TTTTGGGG"},
        nodes=[{"node_id": "n1"}, {"node_id": "n2"}],
        edges=[{"edge_id": "e1", "source": "n1", "target": "n2"}],
        node_header=["node_id"],
        edge_header=["edge_id", "source", "target"],
    )


@pytest.fixture
def phage_build(monkeypatch: pytest.MonkeyPatch) -> dict[str, list]:
    """Replace the genomes check, the tools, the baseline script, and the gold painter."""
    calls: dict[str, list] = {"pin": [], "ctx": []}

    def run(command: list[str], check: bool) -> object:
        dump_cfa(_graph(), Path(command[command.index("--work") + 1]) / "coloured_cfa")
        return type("Completed", (), {"returncode": 0})()

    def pin(bench: str, graph: CfaGraph, fasta: Path) -> Path:
        calls["pin"].append(bench)
        return Path("/pinned/nodes.kraken")

    def materialize(graph: CfaGraph, outdir: Path, ctx: object, **_kwargs: object) -> object:
        calls["ctx"].append(ctx)
        raise Stop

    monkeypatch.setattr("metametro.bench.genomes.check_genomes", lambda *_args: None)
    monkeypatch.setattr(build_module, "require_programs", lambda *_args: None)
    monkeypatch.setattr(build_module.subprocess, "run", run)
    monkeypatch.setattr("metametro.bench.phage_gold.paint_phage_gold", lambda graph, *_args: graph)
    monkeypatch.setattr("metametro.bench.kraken_pin.checked_kraken_pin", pin)
    monkeypatch.setattr(build_module, "materialize", materialize)
    return calls


@pytest.mark.parametrize(
    ("colourings", "pinned"),
    [(None, True), (("as_built", "kraken2"), True), (("as_built", "composition_kmeans"), False)],
)
def test_phage_build_uses_the_kraken_pin_only_when_kraken2_can_run(
    tmp_path: Path, phage_build: dict[str, list], colourings: tuple[str, ...] | None, pinned: bool
) -> None:
    with pytest.raises(Stop):
        build("phage_10", outdir=tmp_path / "phage_10", execute=True, colourings=colourings)
    assert phage_build["pin"] == (["phage_species_5_x10"] if pinned else [])
    expected = Path("/pinned/nodes.kraken") if pinned else None
    assert phage_build["ctx"][0].kraken_output == expected


def test_a_given_kraken_output_wins_over_the_pin(tmp_path: Path, phage_build: dict[str, list]) -> None:
    given = tmp_path / "given.kraken"
    given.write_text("U\tn1\t0\t8\t0:6\n", encoding="utf-8")
    with pytest.raises(Stop):
        build("phage_10", outdir=tmp_path / "phage_10", execute=True, kraken_output=given)
    assert phage_build["pin"] == []
    assert phage_build["ctx"][0].kraken_output == given
