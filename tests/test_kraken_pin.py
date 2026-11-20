"""Mandatory: the pinned Kraken2 output of phage_10 and its input check. No kraken2, no index."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from metametro.bench.kraken_pin import checked_kraken_pin, kraken_pin
from metametro.contracts.classifiers import colour_from_kraken2_file, parse_kraken2_output
from metametro.errors import ContractError
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import write_yaml

pytestmark = pytest.mark.mandatory


def _graph(sequence: str = "ACGTACGT") -> CfaGraph:
    return CfaGraph(
        metadata={"schema_version": "1.0", "graph_id": "toy", "graph_type": "de_bruijn", "k": 3},
        sequences={"n1": sequence, "n2": "TTTTGGGG"},
        nodes=[{"node_id": "n1"}, {"node_id": "n2"}],
        edges=[{"edge_id": "e1", "source": "n1", "target": "n2"}],
        node_header=["node_id"],
        edge_header=["edge_id", "source", "target"],
    )


def _pin(tmp_path: Path) -> Path:
    root = tmp_path / "pins"
    folder = root / "phage_species_5_x10" / "kraken2"
    folder.mkdir(parents=True)
    output = folder / "nodes.kraken"
    output.write_text("C\tn1\t10665\t8\t10665:6\nU\tn2\t0\t8\t0:6\n", encoding="utf-8")
    fasta = ">n1\nACGTACGT\n>n2\nTTTTGGGG\n"
    write_yaml(
        folder / "kraken2.yaml",
        {
            "bench": "phage_species_5_x10",
            "input_sha256": hashlib.sha256(fasta.encode()).hexdigest(),
            "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        },
    )
    return root


def test_phage_10_pin_records_tool_index_and_hashes() -> None:
    found = kraken_pin("phage_10")
    assert found is not None
    output, meta = found
    assert meta["bench"] == "phage_species_5_x10"
    assert meta["tool"] == "kraken2" and meta["version"] == "2.17.2"
    assert meta["index"]["name"] == "k2_viral_20261022"
    assert len(meta["index"]["tar_sha256"]) == 64
    calls = parse_kraken2_output(output)
    assert len(calls) == meta["input_sequences"] == 983
    assert sum(1 for line in output.read_text().splitlines() if line.startswith("C\t")) == meta["classified"]
    assert kraken_pin("phage_species_5") is None


def test_pin_is_used_only_for_the_graph_it_was_made_from(tmp_path: Path) -> None:
    pins = _pin(tmp_path)
    output = checked_kraken_pin("phage_10", _graph(), tmp_path / "work" / "nodes.fna", pins=pins)
    assert output == pins / "phage_species_5_x10" / "kraken2" / "nodes.kraken"
    coloured = colour_from_kraken2_file(_graph(), output)
    assert [(row["namespace"], row["value"]) for row in coloured.colors or []] == [("kraken2", "10665")]
    with pytest.raises(ContractError, match="--kraken-db or --kraken-output"):
        checked_kraken_pin("phage_10", _graph("ACGTACGA"), tmp_path / "other.fna", pins=pins)
    assert checked_kraken_pin("phage_species_5", _graph(), tmp_path / "x.fna", pins=pins) is None
