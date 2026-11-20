"""Mandatory: what the phage branch of ``build`` hands to its steps, with every tool replaced. No network."""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from metametro.bench.build import build
from metametro.bench.phage_gold import Membership
from metametro.bench.ground_truth import gold_digest
from metametro.formats.cfa.io import dump_cfa
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import read_yaml

pytestmark = pytest.mark.mandatory

build_module = importlib.import_module("metametro.bench.build")


class Stop(Exception):
    """Raised by the fake ``materialize`` to end the build after the wiring under test."""


def _graph() -> CfaGraph:
    return CfaGraph(
        metadata={"schema_version": "1.0", "graph_id": "toy_phage", "graph_type": "de_bruijn", "k": 3},
        sequences={"n1": "ACGTACGT", "n2": "GTTTGGGG"},
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

    monkeypatch.setattr("metametro.bench.genomes.check_genomes", lambda *_args: None)
    monkeypatch.setattr(build_module, "require_programs", lambda *_args: None)
    monkeypatch.setattr(build_module.subprocess, "run", run)
    monkeypatch.setattr("metametro.bench.phage_gold.paint_phage_gold", lambda graph, *_args, **_kwargs: graph)
    monkeypatch.setattr(
        "metametro.bench.phage_gold.phage_membership",
        lambda graph, *_args: Membership(
            genomes=("T4",), positions={"n1": 6, "n2": 6}, counts={"n1": (6,), "n2": (0,)}
        ),
    )
    monkeypatch.setattr("metametro.bench.kraken_pin.checked_kraken_pin", pin)
    return calls


@pytest.fixture
def stop_at_materialize(phage_build: dict[str, list], monkeypatch: pytest.MonkeyPatch) -> dict[str, list]:
    """End the build at ``materialize``, recording the context it was given."""

    def materialize(graph: CfaGraph, outdir: Path, ctx: object, **_kwargs: object) -> object:
        phage_build["ctx"].append(ctx)
        raise Stop

    monkeypatch.setattr(build_module, "materialize", materialize)
    return phage_build


@pytest.mark.parametrize(
    ("colourings", "pinned"),
    [(None, True), (("as_built", "kraken2"), True), (("as_built", "composition_kmeans"), False)],
)
def test_phage_build_uses_the_kraken_pin_only_when_kraken2_can_run(
    tmp_path: Path, stop_at_materialize: dict[str, list], colourings: tuple[str, ...] | None, pinned: bool
) -> None:
    with pytest.raises(Stop):
        build("phage_10", outdir=tmp_path / "phage_10", execute=True, colourings=colourings)
    assert stop_at_materialize["pin"] == (["phage_species_5_x10"] if pinned else [])
    expected = Path("/pinned/nodes.kraken") if pinned else None
    assert stop_at_materialize["ctx"][0].kraken_output == expected


def test_a_given_kraken_output_wins_over_the_pin(tmp_path: Path, stop_at_materialize: dict[str, list]) -> None:
    given = tmp_path / "given.kraken"
    given.write_text("U\tn1\t0\t8\t0:6\n", encoding="utf-8")
    with pytest.raises(Stop):
        build("phage_10", outdir=tmp_path / "phage_10", execute=True, kraken_output=given)
    assert stop_at_materialize["pin"] == []
    assert stop_at_materialize["ctx"][0].kraken_output == given


def test_phage_build_writes_taxon_values_and_bubbles_before_identity_and_manifest(
    tmp_path: Path, phage_build: dict[str, list], monkeypatch: pytest.MonkeyPatch
) -> None:
    outdir = tmp_path / "phage_10"
    seen: list[tuple[bool, bool, bool, bool]] = []

    def bubbles(path: Path) -> Path:
        seen.append(
            (
                (path / "ground_truth" / "taxon_values.tsv").is_file(),
                (path / "ground_truth" / "taxon_membership.tsv").is_file(),
                (path / "identity.sha256").exists(),
                (path / "manifest.yaml").exists(),
            )
        )
        return path / "ground_truth" / "bubbles.tsv"

    monkeypatch.setattr("metametro.bench.bubble_gold.write_phage_bubbles", bubbles)
    result = build("phage_10", outdir=outdir, execute=True, colourings=("as_built",))
    assert seen == [(True, True, False, False)]
    assert result.status == "built" and (outdir / "identity.sha256").is_file()
    assert (outdir / "ground_truth" / "taxon_membership.tsv").read_text().splitlines() == [
        "node_id\tkmers\tT4\tT4_fraction",
        "n1\t6\t1\t1.000000",
        "n2\t6\t0\t0.000000",
    ]
    manifest = read_yaml(outdir / "manifest.yaml")
    assert manifest["gold_sha256"] == gold_digest(outdir)
    assert "ground_truth/taxon_membership.tsv" in (outdir / "gold.sha256").read_text()


def test_failed_phage_ground_truth_leaves_no_identity_or_manifest(
    tmp_path: Path, phage_build: dict[str, list], monkeypatch: pytest.MonkeyPatch
) -> None:
    outdir = tmp_path / "phage_10"

    def bubbles(path: Path) -> Path:
        raise RuntimeError("bubble gold failed")

    monkeypatch.setattr("metametro.bench.bubble_gold.write_phage_bubbles", bubbles)
    with pytest.raises(RuntimeError, match="bubble gold failed"):
        build("phage_10", outdir=outdir, execute=True, colourings=("as_built",))
    assert (outdir / "cdbg").is_dir()
    assert not (outdir / "identity.sha256").exists()
    assert not (outdir / "manifest.yaml").exists()
    assert not (outdir / "gold.sha256").exists()
