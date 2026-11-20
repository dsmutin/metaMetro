"""Mandatory: pinned phage genomes, the missing-genome report, and the hash check. No network."""

from __future__ import annotations

import hashlib
import importlib
from pathlib import Path

import pytest

import metametro.bench.genomes as genomes
from metametro.bench import build
from metametro.bench.genomes import check_genomes, fetch_genomes, pinned_genomes, pipeline_fasta
from metametro.cli import main
from metametro.errors import ContractError
from metametro.tables import read_yaml

pytestmark = pytest.mark.mandatory

# ``metametro.bench.build`` the attribute is the function; this is the module.
build_module = importlib.import_module("metametro.bench.build")

_FAKE = {"T1": "ACGTACGTAA", "T3": "CCGGTTAAGG", "T4": "ACGMTTSAGR", "T5": "GGGGCCCCAA", "T7": "TTTTAAAACC"}


def _ncbi(accession: str, sequence: str) -> str:
    return f">{accession} fake phage, complete genome\n{sequence[:6]}\n{sequence[6:].lower()}\n\n"


@pytest.fixture
def fake_pins(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, str]:
    """Point the pin at five tiny fake genomes and return the NCBI text by URL."""
    pin = tmp_path / "pins" / "phage_species_5"
    pin.mkdir(parents=True)
    lines = ["\t".join(genomes.COLUMNS)]
    responses: dict[str, str] = {}
    for index, (name, sequence) in enumerate(_FAKE.items()):
        accession = f"NC_99999{index}.1"
        raw = _ncbi(accession, sequence)
        fasta = pipeline_fasta(name, raw)
        responses[genomes.EFETCH.format(accession=accession)] = raw
        lines.append(
            "\t".join(
                [
                    name,
                    accession,
                    str(100 + index),
                    f"fake phage {name}",
                    str(200 + index),
                    f"fake species {name}",
                    str(len(sequence)),
                    hashlib.sha256(raw.encode()).hexdigest(),
                    hashlib.sha256(fasta.encode()).hexdigest(),
                    "",
                ]
            )
        )
    (pin / "genomes.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    monkeypatch.setattr(genomes, "PINS", tmp_path / "pins")
    return responses


def test_phage_pin_lists_five_versioned_refseq_genomes() -> None:
    """Both phage benches share one pin, and the contract repeats its accessions."""
    pinned = pinned_genomes("phage_10")
    assert [item.name for item in pinned] == ["T1", "T3", "T4", "T5", "T7"]
    assert [item.accession for item in pinned] == [
        "NC_005833.1",
        "NC_003298.1",
        "NC_000866.4",
        "NC_005859.1",
        "NC_001604.1",
    ]
    assert pinned == pinned_genomes("phage_species_5")
    assert all(len(item.fasta_sha256) == 64 and item.taxid > 0 for item in pinned)
    with pytest.raises(ContractError, match="has no pinned genomes"):
        pinned_genomes("bubble_strain_2")


def test_contract_parameters_come_from_the_pin(tmp_path: Path) -> None:
    """``parameters.yaml`` lists the pinned names and accessions in pin order."""
    result = build("phage_species_5", outdir=tmp_path / "phage", execute=False)
    payload = read_yaml(result.outdir / "contract" / "parameters.yaml")
    pinned = pinned_genomes("phage_species_5")
    assert payload["genomes"] == [item.name for item in pinned]
    assert payload["accessions"] == [item.accession for item in pinned]


def test_pipeline_fasta_renames_and_masks_iupac() -> None:
    """The record id becomes the genome name; non-ACGTN bases become N on one line."""
    assert pipeline_fasta("T4", ">NC_1.1 x\nACGM\ntsr\n") == ">T4\nACGNTNN\n"
    with pytest.raises(ContractError, match="one FASTA record"):
        pipeline_fasta("T4", ">a\nAC\n>b\nGT\n")


def test_missing_genomes_report_names_each_one_and_the_fetch_command(tmp_path: Path) -> None:
    """An empty genome directory lists all five genomes with accessions and the command."""
    with pytest.raises(ContractError) as caught:
        check_genomes("phage_10", tmp_path)
    message = str(caught.value)
    for name, accession in (("T1", "NC_005833.1"), ("T4", "NC_000866.4"), ("T7", "NC_001604.1")):
        assert f"{name} ({accession})" in message
    assert "metametro benchfetch phage_species_5_x10" in message


def test_changed_genome_is_reported_with_force(tmp_path: Path, fake_pins: dict[str, str]) -> None:
    """Only the missing and the altered genomes are named; a changed file asks for --force."""
    fetch_genomes("phage_species_5", tmp_path, download=lambda url: fake_pins[url].encode(), pause=0)
    check_genomes("phage_species_5", tmp_path)
    (tmp_path / "T3.fna").unlink()
    (tmp_path / "T5.fna").write_text(">T5\nAAAA\n", encoding="utf-8")
    with pytest.raises(ContractError) as caught:
        check_genomes("phage_species_5", tmp_path)
    errors = caught.value.errors
    assert errors[0].endswith(": T3 (NC_999991.1)")
    assert errors[1].endswith(": T5 (NC_999993.1)")
    assert errors[2] == "fetch them with: metametro benchfetch phage_species_5 --force"


def test_fetch_is_idempotent_and_checks_the_hash(tmp_path: Path, fake_pins: dict[str, str]) -> None:
    """A second fetch downloads nothing; a changed file needs force; a bad download writes nothing."""
    calls: list[str] = []

    def download(url: str) -> bytes:
        calls.append(url)
        return fake_pins[url].encode()

    first = fetch_genomes("phage_species_5", tmp_path, download=download, pause=0)
    assert set(first.values()) == {"fetched"}
    assert (tmp_path / "T4.fna").read_text(encoding="utf-8") == ">T4\nACGNTTNAGN\n"
    assert len(calls) == 5
    again = fetch_genomes("phage_species_5", tmp_path, download=download, pause=0)
    assert set(again.values()) == {"present"}
    assert len(calls) == 5
    (tmp_path / "T1.fna").write_text(">T1\nAAAA\n", encoding="utf-8")
    with pytest.raises(ContractError, match="--force"):
        fetch_genomes("phage_species_5", tmp_path, download=download, pause=0)
    forced = fetch_genomes("phage_species_5", tmp_path, download=download, force=True, pause=0)
    assert forced["T1"] == "replaced"
    (tmp_path / "T7.fna").unlink()

    def corrupt(url: str) -> bytes:
        return fake_pins[url].replace("TTTTAA", "TTTTAT").encode()

    with pytest.raises(ContractError, match="SHA-256"):
        fetch_genomes("phage_species_5", tmp_path, download=corrupt, pause=0)
    assert not (tmp_path / "T7.fna").exists()
    assert not list(tmp_path.glob("*.part"))


def test_phage_build_without_genomes_stops_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``build(..., execute=True)`` raises the missing-genome report and creates no output."""
    monkeypatch.setattr(build_module, "repo_root", lambda: tmp_path / "checkout")
    outdir = tmp_path / "out"
    with pytest.raises(ContractError, match=r"missing pinned genome FASTA.*benchfetch phage_species_5_x10"):
        build("phage_10", outdir=outdir, execute=True)
    assert not outdir.exists()


def test_benchfetch_list_prints_the_pin(capsys: pytest.CaptureFixture[str]) -> None:
    """``metametro benchfetch NAME --list`` prints the pin without downloading."""
    assert main(["benchfetch", "phage_10", "--list"]) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0] == "name\taccession\ttaxid\torganism\tlength"
    assert lines[1].startswith("T1\tNC_005833.1\t2492962\t")
    assert len(lines) == 6

