"""Pinned reference genomes for the phage benchmarks, and their download.

``data/pins/phage_species_5/genomes.tsv`` pins one RefSeq record per genome:
the versioned accession, the NCBI taxid, and the SHA-256 of the FASTA file the
build reads. ``fetch_genomes`` downloads each record from NCBI E-utilities and
writes ``<genome_dir>/<name>.fna`` with the record id set to the genome name and
every base that is not ``ACGTN`` replaced by ``N``. ``check_genomes`` raises
``ContractError`` that names each missing or changed genome and the command that
fetches it.
"""

from __future__ import annotations

import hashlib
import re
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlopen

from metametro.bench.registry import resolve
from metametro.errors import ContractError
from metametro.tables import read_tsv

PINS = Path(__file__).resolve().parent / "data" / "pins"
EFETCH = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={accession}&rettype=fasta&retmode=text&tool=metametro"
)
COLUMNS = (
    "name",
    "accession",
    "taxid",
    "organism",
    "species_taxid",
    "species",
    "length",
    "source_sha256",
    "fasta_sha256",
    "note",
)
_NOT_ACGTN = re.compile(r"[^ACGTN]")


@dataclass(frozen=True)
class PinnedGenome:
    """One row of a ``genomes.tsv`` pin."""

    name: str
    accession: str
    taxid: int
    organism: str
    species_taxid: int
    species: str
    length: int
    source_sha256: str
    fasta_sha256: str
    note: str


def pin_name(bench: str) -> str:
    """Return the pin directory for ``bench``. Both phage benches share one pin."""
    spec = resolve(bench)
    if spec.name.startswith("phage_species_5"):
        return "phage_species_5"
    raise ContractError([f"{spec.name} has no pinned genomes; only phage_species_5 and phage_species_5_x10 do"])


def pinned_genomes(bench: str = "phage_species_5", *, pins: Path | None = None) -> tuple[PinnedGenome, ...]:
    """Read the pinned genome table for ``bench`` (from ``PINS`` unless ``pins`` is given)."""
    path = (PINS if pins is None else Path(pins)) / pin_name(bench) / "genomes.tsv"
    header, rows = read_tsv(path)
    if tuple(header) != COLUMNS:
        raise ContractError([f"{path} columns are {header}, expected {list(COLUMNS)}"])
    return tuple(
        PinnedGenome(
            name=row["name"],
            accession=row["accession"],
            taxid=int(row["taxid"]),
            organism=row["organism"],
            species_taxid=int(row["species_taxid"]),
            species=row["species"],
            length=int(row["length"]),
            source_sha256=row["source_sha256"],
            fasta_sha256=row["fasta_sha256"],
            note=row["note"],
        )
        for row in rows
    )


def sha256_file(path: Path) -> str:
    """Return the hex SHA-256 of ``path``."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pipeline_fasta(name: str, ncbi_fasta: str) -> str:
    """Turn one NCBI FASTA record into the file the build reads.

    The record id becomes ``name``, the sequence is one upper-case line, and
    every base that is not ``ACGTN`` becomes ``N``.
    """
    lines = [line.strip() for line in ncbi_fasta.splitlines() if line.strip()]
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) != 1 or not lines[0].startswith(">"):
        raise ContractError([f"{name}: expected one FASTA record, found {len(headers)}"])
    sequence = _NOT_ACGTN.sub("N", "".join(lines[1:]).upper())
    return f">{name}\n{sequence}\n"


def genome_problems(bench: str, genome_dir: Path) -> tuple[list[PinnedGenome], list[PinnedGenome]]:
    """Return the pinned genomes that are missing and those whose hash differs."""
    missing: list[PinnedGenome] = []
    changed: list[PinnedGenome] = []
    for genome in pinned_genomes(bench):
        path = Path(genome_dir) / f"{genome.name}.fna"
        if not path.is_file():
            missing.append(genome)
        elif sha256_file(path) != genome.fasta_sha256:
            changed.append(genome)
    return missing, changed


def check_genomes(bench: str, genome_dir: Path) -> None:
    """Raise ``ContractError`` unless every pinned genome is in ``genome_dir`` with its pinned hash."""
    missing, changed = genome_problems(bench, genome_dir)
    if not missing and not changed:
        return
    spec = resolve(bench)
    errors: list[str] = []
    if missing:
        listed = ", ".join(f"{item.name} ({item.accession})" for item in missing)
        errors.append(f"{spec.name} is missing pinned genome FASTA under {genome_dir}: {listed}")
    if changed:
        listed = ", ".join(f"{item.name} ({item.accession})" for item in changed)
        errors.append(f"{spec.name} genome FASTA under {genome_dir} differs from the pinned SHA-256: {listed}")
    command = f"metametro benchfetch {spec.name}"
    if changed:
        command += " --force"
    errors.append(f"fetch them with: {command}")
    raise ContractError(errors)


def _download(url: str) -> bytes:
    with urlopen(url, timeout=120) as response:  # noqa: S310 - fixed https endpoint
        return response.read()


def fetch_genomes(
    bench: str,
    genome_dir: Path,
    *,
    force: bool = False,
    download: Callable[[str], bytes] = _download,
    pause: float = 0.4,
) -> dict[str, str]:
    """Download every pinned genome for ``bench`` into ``genome_dir``.

    A file that already has the pinned hash is kept (status ``present``). A
    file with another hash raises ``ContractError`` unless ``force`` is true.
    A download whose result does not match the pinned length and hash raises
    ``ContractError`` and leaves no file behind. Returns status by genome name.
    """
    genome_dir = Path(genome_dir)
    status: dict[str, str] = {}
    _, changed = genome_problems(bench, genome_dir)
    if changed and not force:
        listed = ", ".join(f"{item.name}.fna" for item in changed)
        raise ContractError(
            [f"{genome_dir} has genome FASTA with another SHA-256: {listed}; rerun with --force to replace them"]
        )
    first = True
    for genome in pinned_genomes(bench):
        path = genome_dir / f"{genome.name}.fna"
        if path.is_file() and sha256_file(path) == genome.fasta_sha256:
            status[genome.name] = "present"
            continue
        if not first and pause > 0:
            time.sleep(pause)
        first = False
        url = EFETCH.format(accession=genome.accession)
        raw = download(url).decode("utf-8")
        header = raw.lstrip().split("\n", 1)[0]
        if not header.startswith(f">{genome.accession}"):
            raise ContractError([f"{genome.name}: NCBI returned {header[:80]!r}, expected >{genome.accession}"])
        text = pipeline_fasta(genome.name, raw)
        length = len(text.split("\n")[1])
        if length != genome.length:
            raise ContractError([f"{genome.name} ({genome.accession}): length {length}, pinned {genome.length}"])
        digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        if digest != genome.fasta_sha256:
            raise ContractError(
                [f"{genome.name} ({genome.accession}): SHA-256 {digest}, pinned {genome.fasta_sha256}"]
            )
        genome_dir.mkdir(parents=True, exist_ok=True)
        partial = path.with_suffix(".fna.part")
        partial.write_text(text, encoding="utf-8")
        partial.replace(path)
        status[genome.name] = "replaced" if genome in changed else "fetched"
    return status
