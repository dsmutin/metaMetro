"""Mandatory: the pinned Kraken2 output of phage_10 and its input check. No kraken2, no index."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from metametro.bench.genomes import pinned_genomes
from metametro.bench.kraken_pin import checked_kraken_pin, kraken_pin, pinned_taxdump
from metametro.contracts.classifiers import colour_from_kraken2_file, parse_kraken2_output
from metametro.errors import ContractError
from metametro.formats.cfa.model import CfaGraph
from metametro.tables import read_yaml, write_yaml

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


def _dmp(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").rstrip("|").rstrip("\t").split("\t|\t") for line in path.read_text().splitlines(True)]


def test_pinned_taxdump_has_every_lineage_to_root() -> None:
    folder = pinned_taxdump("phage_10")
    assert folder is not None and folder == kraken_pin("phage_10")[1]["taxdump"]["dir"]
    parent = {int(row[0]): int(row[1]) for row in _dmp(folder / "nodes.dmp")}
    named = {int(row[0]) for row in _dmp(folder / "names.dmp") if row[3] == "scientific name"}
    output = kraken_pin("phage_10")[0]
    wanted: set[int] = set()
    for line in output.read_text().splitlines():
        fields = line.split("\t")
        wanted.add(int(fields[2]))
        wanted.update(int(token.split(":")[0]) for token in fields[4].split() if token.split(":")[0].isdigit())
    for genome in pinned_genomes("phage_10"):
        wanted.update((genome.taxid, genome.species_taxid))
    wanted.discard(0)
    for taxon in wanted:
        seen = set()
        while taxon != 1:
            assert taxon in parent and taxon in named, taxon
            assert taxon not in seen
            seen.add(taxon)
            taxon = parent[taxon]
    assert parent[1] == 1 and 1 in named
    assert len(parent) == kraken_pin("phage_10")[1]["taxdump"]["taxids"]
    assert pinned_taxdump("phage_species_5") is None


def test_incomplete_pin_is_an_error(tmp_path: Path) -> None:
    pins = _pin(tmp_path)
    folder = pins / "phage_species_5_x10" / "kraken2"
    (folder / "taxdump").mkdir()
    with pytest.raises(ContractError, match="no taxdump block"):
        kraken_pin("phage_10", pins=pins)
    with pytest.raises(ContractError, match="no taxdump block"):
        pinned_taxdump("phage_10", pins=pins)
    (folder / "taxdump").rmdir()
    meta = folder / "kraken2.yaml"
    record = meta.read_text()
    taxdump = {"path": "taxdump", "nodes_sha256": "0" * 64, "names_sha256": "0" * 64}
    write_yaml(meta, {**read_yaml(meta), "taxdump": taxdump})
    with pytest.raises(ContractError, match="nodes.dmp is missing"):
        pinned_taxdump("phage_10", pins=pins)
    meta.write_text(record)
    (folder / "nodes.kraken").unlink()
    with pytest.raises(ContractError, match="nodes.kraken is missing"):
        kraken_pin("phage_10", pins=pins)
    meta.unlink()
    assert kraken_pin("phage_10", pins=pins) is None
    _pin(tmp_path / "again")
    (tmp_path / "again" / "pins" / "phage_species_5_x10" / "kraken2" / "kraken2.yaml").unlink()
    with pytest.raises(ContractError, match="kraken2.yaml is missing"):
        kraken_pin("phage_10", pins=tmp_path / "again" / "pins")


def test_package_data_ships_every_pin_and_schema() -> None:
    import tomllib

    root = Path(__file__).resolve().parents[1]
    package = root / "src" / "metametro"
    globs = tomllib.loads((root / "pyproject.toml").read_text())["tool"]["setuptools"]["package-data"]["metametro"]
    shipped = {path for pattern in globs for path in package.glob(pattern) if path.is_file()}
    data = {
        path
        for path in package.rglob("*")
        if path.is_file() and path.suffix not in {".py", ".pyc"} and "__pycache__" not in path.parts
    }
    assert data and data <= shipped
    assert package / "bench/data/pins/phage_species_5_x10/kraken2/taxdump/nodes.dmp" in shipped


def test_tampered_pinned_output_is_rejected(tmp_path: Path) -> None:
    pins = _pin(tmp_path)
    output = pins / "phage_species_5_x10" / "kraken2" / "nodes.kraken"
    output.write_text(output.read_text().replace("10665", "10666"), encoding="utf-8")
    with pytest.raises(ContractError, match="does not have the SHA-256 pinned"):
        kraken_pin("phage_10", pins=pins)
    with pytest.raises(ContractError, match="does not have the SHA-256 pinned"):
        checked_kraken_pin("phage_10", _graph(), tmp_path / "nodes.fna", pins=pins)
