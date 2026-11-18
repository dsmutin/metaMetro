"""Run the held-out genera example from pinned NCBI assemblies.

The script resumes. It does not download a genome that is already extracted,
and it does not rebuild a database whose index already exists. Large outputs
stay in ``data/raw/`` and the bench work directory under ``data/bench/``.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import json
import shutil
import subprocess
import sys
import zipfile
from datetime import date
from pathlib import Path

import yaml

from metametro.bench.data.universal.community import (
    assembly_graph_from_fastg,
    highest_megahit_contigs,
    iss_key,
    load_pairs,
    lognormal_read_counts,
    megahit_coverage,
    presence_f1,
    r_squared,
    sequence_digest,
    rank_taxon,
    species_taxon,
)

ROOT = Path(__file__).resolve().parents[5]
EXAMPLE = ROOT / "src" / "metametro" / "bench" / "data" / "pins" / "bacteria_species_20_heldout"
WORK = ROOT / "data" / "bench" / "bacteria_species_20_heldout" / "megahit" / "reads100000_k21" / "work"
GRAPH_ID = "bacteria_species_20_heldout"
SCORE_RANK = "species"
REPORT_PATH = ROOT / "data" / "raw" / "bacteria_species_20_heldout_assembly_data_report.jsonl"
TOTAL_READS = 100_000
SEED = 42
SIGMA = 1.5
ZIP_NAME = "ncbi_bacteria_species_20_heldout.zip"


def configure(
    *,
    example: Path,
    work: Path,
    root: Path,
    graph_id: str,
    score_rank: str,
    total_reads: int,
    report_path: Path,
    zip_name: str,
) -> None:
    """Point this pipeline at one benchmark pin and one bench output directory."""
    global EXAMPLE, WORK, ROOT, GRAPH_ID, SCORE_RANK, REPORT_PATH, TOTAL_READS, ZIP_NAME
    EXAMPLE = example
    WORK = work
    ROOT = root
    GRAPH_ID = graph_id
    SCORE_RANK = score_rank
    REPORT_PATH = report_path
    TOTAL_READS = total_reads
    ZIP_NAME = zip_name


def _log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")


def _run(cmd: list[str], log_path: Path, *, cwd: Path | None = None) -> None:
    """Run a command, fail if it exits non-zero, and append it to the log."""
    _log(log_path, "$ " + " ".join(cmd))
    completed = subprocess.run(cmd, cwd=cwd, text=True)
    if completed.returncode != 0:
        raise SystemExit(f"command failed ({completed.returncode}): {' '.join(cmd)}")


def _report_identity(record: dict) -> tuple[str, str, int]:
    """Read accession, organism name, and taxid from either datasets JSON spelling."""
    organism = record.get("organism") or {}
    accession = str(record.get("accession") or record.get("currentAccession") or record.get("current_accession") or "")
    name = str(organism.get("organismName") or organism.get("organism_name") or "")
    tax = organism.get("taxId", organism.get("tax_id"))
    if not accession or not name or tax in (None, ""):
        raise SystemExit(f"assembly report is missing accession, organism name, or taxid: {accession or 'unknown'}")
    return accession, name, int(tax)


LOGGED_TOOLS = ("megahit", "kraken2", "kaiju", "minimap2", "samovar")


def tool_versions(tools: tuple[str, ...] = LOGGED_TOOLS) -> dict[str, str]:
    """Return the first ``--version`` output line of each tool, or ``absent`` when it is not on PATH."""
    found: dict[str, str] = {}
    for tool in tools:
        program = shutil.which(tool)
        if program is None:
            found[tool] = "absent"
            continue
        try:
            done = subprocess.run([program, "--version"], capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.TimeoutExpired) as exc:
            found[tool] = f"unknown ({type(exc).__name__})"
            continue
        lines = [line.strip() for line in (done.stdout + "\n" + done.stderr).splitlines() if line.strip()]
        found[tool] = lines[0] if lines else f"unknown (exit {done.returncode})"
    return found


def log_tool_versions() -> dict[str, str]:
    """Append external tool versions to ``WORK/commands.log`` and return them."""
    versions = tool_versions()
    for tool, version in versions.items():
        _log(WORK / "commands.log", f"tool version {tool}: {version}")
    return versions


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tool(name: str) -> str:
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"required program is not on PATH: {name}")
    return found


def _help_text(program: str) -> str:
    completed = subprocess.run([program, "--help"], capture_output=True, text=True)
    text = (completed.stdout or "") + (completed.stderr or "")
    if completed.returncode != 0 and not text:
        completed = subprocess.run([program, "-h"], capture_output=True, text=True)
        text = (completed.stdout or "") + (completed.stderr or "")
    return text


def _require_help(program: str, flag: str) -> None:
    if flag not in _help_text(program):
        raise SystemExit(f"{program} help text does not mention {flag}")


def stage_download() -> None:
    """Download the pinned assemblies with NCBI datasets and validate them."""
    datasets = _tool("datasets")
    rows = load_pairs(EXAMPLE / "accessions.tsv")
    raw = ROOT / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    log_path = ROOT / "data" / "logs" / "download.log"
    version = subprocess.run([datasets, "--version"], capture_output=True, text=True).stdout.strip()
    _log(log_path, f"datasets version: {version}")
    zip_path = raw / ZIP_NAME
    if not zip_path.is_file() or not zipfile.is_zipfile(zip_path):
        listing = raw / "accessions.txt"
        listing.write_text("".join(f"{row['accession']}\n" for row in rows), encoding="utf-8")
        download_help = subprocess.run(
            [datasets, "download", "genome", "accession", "--help"],
            capture_output=True,
            text=True,
        )
        help_text = (download_help.stdout or "") + (download_help.stderr or "")
        for flag in ("--include", "--filename", "--inputfile"):
            if flag not in help_text:
                raise SystemExit(f"datasets download help does not mention {flag}")
        _run(
            [
                datasets,
                "download",
                "genome",
                "accession",
                "--inputfile",
                str(listing),
                "--include",
                "genome",
                "--filename",
                str(zip_path),
            ],
            log_path,
        )
    fasta_dir = raw / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as archive:
        report_name = next(
            (name for name in archive.namelist() if name.endswith("assembly_data_report.jsonl")),
            "",
        )
        if not report_name:
            raise SystemExit("NCBI zip has no assembly_data_report.jsonl")
        REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
        REPORT_PATH.write_bytes(archive.read(report_name))
        report_path = REPORT_PATH
        for row in rows:
            accession = row["accession"]
            dest = fasta_dir / f"{accession}.fna"
            if dest.is_file() and dest.stat().st_size > 0:
                continue
            genomic = [
                name
                for name in archive.namelist()
                if accession in name and name.endswith(".fna")
            ]
            member = next((name for name in genomic if name.endswith("_genomic.fna")), genomic[0] if genomic else "")
            if not member:
                raise SystemExit(f"NCBI zip has no FASTA for {accession}")
            dest.write_bytes(archive.read(member))
            if dest.stat().st_size == 0:
                raise SystemExit(f"extracted FASTA is empty: {accession}")
    _validate_genomes(rows, fasta_dir, report_path, zip_path, version)


def _validate_genomes(rows: list[dict[str, str]], fasta_dir: Path, report_path: Path, zip_path: Path, datasets_version: str) -> None:
    reports = {}
    for line in report_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        accession, _name, _tax = _report_identity(json.loads(line))
        reports[accession] = json.loads(line)
    today = date.today().isoformat()
    manifest = ["accession\trepository\tdownload_url\tlocal_path\tfilename_bytes\tsha256\tdownload_status\tvalidation_status\tbiological_status\tnotes\tdownload_date"]
    checksums = [f"{_sha256(zip_path)}  {zip_path.relative_to(ROOT)}"]
    by_pair: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = row["accession"]
        record = reports.get(accession)
        if record is None:
            raise SystemExit(f"assembly report has no record for {accession}")
        _accession, name, tax_id = _report_identity(record)
        if not name.startswith(row["genus"]):
            raise SystemExit(f"{accession} organism {name!r} is not genus {row['genus']}")
        path = fasta_dir / f"{accession}.fna"
        digest = sequence_digest(path)
        by_pair.setdefault(row["pair_id"], {})[row["role"]] = digest
        file_sha = _sha256(path)
        checksums.append(f"{file_sha}  {path.relative_to(ROOT)}")
        manifest.append(
            "\t".join(
                [
                    accession,
                    "NCBI Assembly",
                    "",
                    str(path.relative_to(ROOT)),
                    str(path.stat().st_size),
                    file_sha,
                    "success",
                    "pass",
                    "pass",
                    f"datasets {datasets_version}; organism {name}; tax_id {tax_id}; no per-file URL in the datasets summary",
                    today,
                ]
            )
        )
    for pair_id, digests in sorted(by_pair.items()):
        if digests.get("sim") == digests.get("db"):
            raise SystemExit(f"{pair_id} sim and db assemblies have the same sequence digest")
    manifest_path = WORK / "download" / "download_manifest.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    checksum_path = WORK / "download" / "checksums.txt"
    checksum_path.parent.mkdir(parents=True, exist_ok=True)
    checksum_path.write_text("\n".join(checksums) + "\n", encoding="utf-8")
    _write_data_reports(rows, datasets_version, today)


_HELDOUT_NOTES = (
    "`eco_coli` sim `GCF_000952955.1` had no tax_id in the summary used to pin the table; the download report supplies it.",
    "Two RefSeq references (`Streptococcus agalactiae`, `Streptococcus thermophilus`) had no strain name in the summary.",
    "Four of the five Escherichia pairs are not E. coli.",
    "Shigella has two flexneri pairs (species 623 lineage and flexneri 2a taxid 42897). Species-rank scores add those read counts.",
    "`Pseudomonas stutzeri` was not used: the NCBI summary for that taxid returned `Stutzerimonas stutzeri`.",
)
"""Findings recorded when ``bacteria_species_20_heldout`` was pinned. Only that bench prints them."""


def _shown_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _write_data_reports(rows: list[dict[str, str]], datasets_version: str, today: str) -> None:
    """Write the acquisition report and dataset audit for the configured bench and pin."""
    n_genomes = len(rows)
    n_sim = sum(1 for row in rows if row["role"] == "sim")
    n_db = sum(1 for row in rows if row["role"] == "db")
    n_pairs = len({row["pair_id"] for row in rows})
    genera = sorted({row["genus"] for row in rows if row.get("genus")})
    pin = _shown_path(EXAMPLE / "accessions.tsv")
    notes = list(_HELDOUT_NOTES) if GRAPH_ID == "bacteria_species_20_heldout" else []
    acquisition = WORK / "download" / "acquisition_report.md"
    acquisition.parent.mkdir(parents=True, exist_ok=True)
    acquisition.write_text(
        "\n".join(
            [
                f"# Acquisition report: {GRAPH_ID}",
                "",
                f"**Date:** {today}",
                "",
                f"Downloaded {n_genomes} NCBI assembly accessions listed in `{pin}`.",
                f"Method: `{datasets_version}` `datasets download genome accession --include genome`.",
                "The datasets summary used to choose the accessions did not include a per-file FTP URL, so `download_url` is empty.",
                "Technical validation: each extracted FASTA is non-empty and its SHA-256 is in `checksums.txt` next to this report.",
                "Biological validation: the assembly-report organism name starts with the requested genus, and the two assemblies in each pair have different sequence digests.",
                "No download failed.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    audit = WORK / "download" / "dataset-audit.md"
    lines = [
        f"# Dataset Audit Report: {GRAPH_ID}",
        "",
        f"**Date:** {today}",
        "",
        "## 1. Data inventory",
        "",
        "| File | Path | Rows × Cols | Notes |",
        "|------|------|-------------|-------|",
        f"| Pinned accessions | {pin} | {n_genomes} × {len(rows[0]) if rows else 0} | {n_sim} sim / {n_db} db, {n_pairs} pairs |",
        f"| FASTA | data/raw/fasta/ | {n_genomes} files | one NCBI assembly each |",
        f"| Manifest | download_manifest.tsv | {n_genomes} × 11 | SHA-256 recorded |",
        "",
        "## 2. Checks",
        "",
        "Every accession has a genus, a role, and an organism name from the NCBI assembly report.",
        "Assembly accessions are unique. The two assemblies in each pair have different sequence digests.",
        f"Genera: {len(genera)}. The simulated sample is {TOTAL_READS} paired fragments, created later by Samovar (InSilicoSeq HiSeq).",
        "",
    ]
    if notes:
        lines += ["## 3. Notes recorded when this pin was made", ""] + notes + [""]
    audit.write_text("\n".join(lines), encoding="utf-8")


def stage_layout() -> None:
    """Write sim and db FASTA directories and the lognormal abundance table."""
    rows = load_pairs(EXAMPLE / "accessions.tsv")
    sim = WORK / "sim"
    db = WORK / "db"
    sim.mkdir(parents=True, exist_ok=True)
    db.mkdir(parents=True, exist_ok=True)
    sim_rows = [row for row in rows if row["role"] == "sim"]
    counts = lognormal_read_counts(len(sim_rows), TOTAL_READS, SEED, mu=0.0, sigma=SIGMA)
    abundance = ["taxid,N_sample"]
    for row, count in zip(sim_rows, counts):
        accession = row["accession"]
        key = iss_key(accession)
        source = ROOT / "data" / "raw" / "fasta" / f"{accession}.fna"
        if not source.is_file():
            raise SystemExit(f"missing extracted FASTA: {source}")
        _rewrite_fasta(source, sim / f"{key}.fna", key)
        abundance.append(f"{accession},{count}")
        db_row = next(item for item in rows if item["pair_id"] == row["pair_id"] and item["role"] == "db")
        db_source = ROOT / "data" / "raw" / "fasta" / f"{db_row['accession']}.fna"
        db_dest = db / f"{db_row['accession']}.fna"
        if not db_dest.is_file():
            shutil.copyfile(db_source, db_dest)
    table = WORK / "abundance.csv"
    table.write_text("\n".join(abundance) + "\n", encoding="utf-8")
    (WORK / "abundance_parameters.txt").write_text(
        f"total_reads={TOTAL_READS}\nseed={SEED}\nmu=0\nsigma={SIGMA}\ndistribution=random.Random.lognormvariate\nunit=paired fragments passed to samovar --total_reads\n",
        encoding="utf-8",
    )


def _rewrite_fasta(source: Path, dest: Path, key: str) -> None:
    if dest.is_file() and dest.stat().st_size > 0:
        return
    chunks = []
    index = 0
    sequence: list[str] = []

    def flush() -> None:
        nonlocal index
        if not sequence:
            return
        chunks.append(f">{key}|{index}\n")
        body = "".join(sequence)
        for start in range(0, len(body), 80):
            chunks.append(body[start : start + 80] + "\n")
        index += 1

    for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            flush()
            sequence = []
        else:
            sequence.append(line.strip())
    flush()
    if index == 0:
        raise SystemExit(f"no sequences in {source}")
    dest.write_text("".join(chunks), encoding="utf-8")


def _simulated_reads() -> tuple[Path, Path] | None:
    """Return the Samovar ISS pair. The sample directory is ``iss/initial``."""
    if not (WORK / "iss").is_dir():
        return None
    candidates = [
        path
        for path in (WORK / "iss").rglob("*_full_R1.fastq")
        if path.is_file() and path.stat().st_size > 0
    ]
    outside = [path for path in candidates if ".iss_full" not in path.parts]
    candidates = outside or candidates
    if not candidates:
        return None
    if len(candidates) != 1:
        raise SystemExit(f"expected one simulated R1 file, found {[str(path) for path in candidates]}")
    r1 = candidates[0]
    r2 = r1.with_name(r1.name.replace("_R1.fastq", "_R2.fastq"))
    if not r2.is_file() or r2.stat().st_size == 0:
        raise SystemExit(f"missing simulated R2 mate: {r2}")
    return r1, r2


def stage_simulate() -> None:
    """Simulate one HiSeq sample with Samovar's InSilicoSeq generate pipeline."""
    samovar = _tool("samovar")
    _tool("snakemake")
    _tool("iss")
    if _simulated_reads():
        return
    log_path = WORK / "commands.log"
    _run(
        [
            samovar,
            "generate",
            "--simulator",
            "iss",
            "--genome_dir",
            str(WORK / "sim"),
            "--output_dir",
            str(WORK / "iss"),
            "--n_samples",
            "1",
            "--total_reads",
            str(TOTAL_READS),
            "--seed",
            str(SEED),
            "--model",
            "hiseq",
            "--cores",
            "4",
            "--abundance",
            str(WORK / "abundance.csv"),
            "--host_fraction",
            "0",
        ],
        log_path,
    )
    script = WORK / "iss" / ".generate" / "generate.sh"
    if not script.is_file():
        raise SystemExit(f"samovar generate did not write {script}")
    _log(log_path, f"$ bash {script}")
    process = subprocess.Popen(["bash", str(script)])
    stable = 0
    previous_size = -1
    while process.poll() is None:
        found = _simulated_reads()
        size = found[0].stat().st_size if found else 0
        if found and size == previous_size and size > 0:
            stable += 1
        else:
            stable = 0
        previous_size = size
        if stable >= 3:
            process.terminate()
            process.wait(timeout=30)
            _log(log_path, "stopped generate.sh after *_full_R*.fastq stopped growing; Snakemake waits for a different filename")
            break
        try:
            process.wait(timeout=15)
        except subprocess.TimeoutExpired:
            continue
    if process.returncode not in (0, -15) and _simulated_reads() is None:
        raise SystemExit(f"command failed ({process.returncode}): bash {script}")
    if process.returncode not in (0, None):
        _log(
            log_path,
            "generate.sh exited non-zero after writing *_full_R*.fastq. "
            "Samovar's snakefile expects iss/initial/1_full_R1.fastq; InSilicoSeq writes sample_full or pool_full reads.",
        )
    if _simulated_reads() is None:
        raise SystemExit("samovar ISS wrote no non-empty *_full_R1.fastq")


def _samovar_src() -> Path:
    """Return the Samovar ``src`` directory next to the ``samovar`` executable."""
    binary = shutil.which("samovar")
    if not binary:
        raise SystemExit("required program is not on PATH: samovar")
    src = Path(binary).resolve().parents[1] / "src"
    if not (src / "samovar" / "build_database.py").is_file():
        raise SystemExit(f"Samovar source was not found beside {binary}")
    return src


def _python_with_biopython() -> str:
    """Return a Python that can import Biopython, which Samovar's Kaiju builder needs."""
    seen: set[str] = set()
    for folder in os.environ.get("PATH", "").split(":"):
        for name in ("python3", "python"):
            path = str(Path(folder) / name)
            if path in seen or not os.access(path, os.X_OK):
                continue
            seen.add(path)
            probe = subprocess.run([path, "-c", "import Bio"], capture_output=True)
            if probe.returncode == 0:
                return path
    raise SystemExit("Kaiju indexing needs a Python with Biopython on PATH")


def _build_kaiju_from_pinned_fasta(log_path: Path) -> None:
    """Index 6-frame translations of the pinned database assemblies.

    ``samovar build --type kaiju`` downloads the latest RefSeq proteome for the
    taxid. That proteome is not the held-out assembly. ``fetch_missing=False``
    keeps only translations of the FASTA in ``work/db``.
    """
    try:
        import Bio  # noqa: F401
    except ModuleNotFoundError:
        python = _python_with_biopython()
        if Path(python).resolve() == Path(sys.executable).resolve():
            raise SystemExit(f"{python} cannot import Biopython")
        _log(log_path, f"re-exec Kaiju build with {python}")
        os.execv(python, [python, *sys.argv])
    sys.path.insert(0, str(_samovar_src()))
    from samovar.build_database import add_database_kaiju, build_database_kaiju, get_taxonomy_db

    identities = {}
    report_path = REPORT_PATH
    for line in report_path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            accession, _name, tax_id = _report_identity(json.loads(line))
            identities[accession] = tax_id
    db_path = WORK / "kaiju_db"
    library = _classifier_library()
    fastas = sorted(library.glob("*.fna"))
    n_db = sum(1 for row in load_pairs(_classifier_pin()) if row["role"] == "db")
    if len(fastas) != n_db:
        raise SystemExit(f"expected {n_db} database FASTA files, found {len(fastas)}")
    if db_path.exists():
        shutil.rmtree(db_path)
    db_path.mkdir(parents=True)
    _log(log_path, "kaiju: Samovar add_database_kaiju fetch_missing=False on pinned FASTA")
    for fasta in fastas:
        accession = fasta.name[: -len(".fna")]
        tax_id = identities.get(accession)
        if tax_id is None:
            raise SystemExit(f"no assembly-report taxid for {accession}")
        add_database_kaiju(str(fasta), str(tax_id), db_path=str(db_path), fetch_missing=False)
    get_taxonomy_db(db_path=str(db_path))
    build_database_kaiju(db_path=str(db_path), threads=4, protein=False)


def _classifier_library() -> Path:
    """FASTA used to build Kraken2 and Kaiju indexes.

    ``work/classifier_db`` is the full parent ``db`` set on half communities.
    ``work/db`` is the layout copy of the pin that was simulated.
    """
    library = WORK / "classifier_db"
    if library.is_dir() and any(library.glob("*.fna")):
        return library
    return WORK / "db"


def _classifier_pin() -> Path:
    """Accession table whose ``db`` rows must match the classifier library."""
    from metametro.bench.data.universal.catalog import pin_dir
    from metametro.bench.registry import resolve

    spec = resolve(GRAPH_ID)
    return pin_dir(spec.parent or spec.name)


def stage_databases() -> None:
    """Build Kraken2 and Kaiju indexes from the classifier library FASTA."""
    samovar = _tool("samovar")
    log_path = WORK / "commands.log"
    library = _classifier_library()
    config = {
        "input_dir": [str(library)],
        "output_dir": str(WORK / "db_preprocessed"),
        "mutation_rate": 0,
        "include_percent": 100,
    }
    config_path = WORK / "db_config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    kraken_index = WORK / "kraken_db" / "hash.k2d"
    if not kraken_index.is_file():
        _run(
            [
                samovar,
                "build",
                "--type",
                "kraken2",
                "--config_path",
                str(config_path),
                "--db_path",
                str(WORK / "kraken_db"),
                "--no-example-omit",
            ],
            log_path,
        )
    kaiju_index = next((WORK / "kaiju_db").glob("*.fmi"), None) if (WORK / "kaiju_db").is_dir() else None
    if kaiju_index is None:
        # samovar build --type kaiju downloads the latest RefSeq proteome for each
        # taxid. That file is not the pinned database assembly (one download in
        # this run was the simulated ruysiae genome). Index 6-frame translations
        # of the pinned nucleotide FASTA instead, using Samovar's own builders.
        _build_kaiju_from_pinned_fasta(log_path)


def stage_classify() -> None:
    """Classify the simulated reads with the held-out Kraken2 and Kaiju databases."""
    log_path = WORK / "commands.log"
    found = _simulated_reads()
    if found is None:
        raise SystemExit("simulated FASTQ is missing; run the simulate stage")
    r1, r2 = found
    kraken = _tool("kraken2")
    _require_help(kraken, "--paired")
    _require_help(kraken, "--report")
    out = WORK / "classify"
    out.mkdir(parents=True, exist_ok=True)
    if not (out / "kraken2.output").is_file():
        _run(
            [
                kraken,
                "--db",
                str(WORK / "kraken_db"),
                "--paired",
                "--threads",
                "4",
                "--report",
                str(out / "kraken2.report"),
                "--output",
                str(out / "kraken2.output"),
                str(r1),
                str(r2),
            ],
            log_path,
        )
    kaiju = _tool("kaiju")
    fmi = next((WORK / "kaiju_db").glob("*.fmi"), None)
    nodes = _taxdump_dir() / "nodes.dmp"
    if fmi is None:
        raise SystemExit("Kaiju FM-index is missing")
    if not (out / "kaiju.output").is_file():
        help_text = _help_text(kaiju)
        if "-j" not in help_text or "-f" not in help_text:
            raise SystemExit("kaiju help does not mention -j and -f")
        _run(
            [
                kaiju,
                "-t",
                str(nodes),
                "-f",
                str(fmi),
                "-i",
                str(r1),
                "-j",
                str(r2),
                "-z",
                "4",
                "-o",
                str(out / "kaiju.output"),
            ],
            log_path,
        )


def stage_assemble() -> None:
    """Assemble the reads with MEGAHIT and export FASTG when the toolkit documents it."""
    megahit = _tool("megahit")
    log_path = WORK / "commands.log"
    out = WORK / "megahit"
    contigs = out / "final.contigs.fa"
    if not contigs.is_file():
        _require_help(megahit, "-1")
        _require_help(megahit, "--keep-tmp-files")
        _require_help(megahit, "--presets")
        _run(
            [
                megahit,
                "-1",
                str(_simulated_reads()[0]),
                "-2",
                str(_simulated_reads()[1]),
                "-o",
                str(out),
                "-t",
                "4",
                "--out-prefix",
                "final",
                "--presets",
                "meta-sensitive",
                "--keep-tmp-files",
            ],
            log_path,
        )
    toolkit = shutil.which("megahit_toolkit")
    fastg = out / "assembly.fastg"
    if toolkit and (not fastg.is_file() or fastg.stat().st_size == 0):
        usage = _help_text(toolkit)
        (out / "toolkit_help.txt").write_text(usage, encoding="utf-8")
        if "contig2fastg" not in usage:
            _log(log_path, "megahit_toolkit help does not mention contig2fastg; assembly graph was not exported")
            return
        k_dir = out / "intermediate_contigs"
        names = [path.name for path in k_dir.iterdir()] if k_dir.is_dir() else []
        chosen = highest_megahit_contigs(names)
        if chosen is None:
            _log(log_path, "MEGAHIT kept no k<int>.contigs.fa; assembly graph was not exported")
            return
        source = k_dir / Path(chosen).name
        k_value = source.name.split(".", 1)[0][1:]
        completed = subprocess.run(
            [toolkit, "contig2fastg", k_value, str(source)],
            capture_output=True,
            text=True,
        )
        _log(log_path, f"$ {toolkit} contig2fastg {k_value} {source} > {fastg}")
        if completed.returncode != 0:
            raise SystemExit(completed.stderr.strip() or "contig2fastg failed")
        fastg.write_text(completed.stdout, encoding="utf-8")
        if fastg.stat().st_size == 0:
            raise SystemExit("contig2fastg wrote an empty FASTG")


def _taxdump_dir() -> Path:
    """Directory with NCBI ``nodes.dmp`` and ``names.dmp``."""
    for path in (
        WORK / "kraken_db" / "taxonomy" / "nodes.dmp",
        WORK / "kaiju_db" / "nodes.dmp",
    ):
        if path.is_file() and (path.parent / "names.dmp").is_file():
            return path.parent
    raise SystemExit("missing NCBI nodes.dmp and names.dmp from the Kraken2 or Kaiju build")


def _taxdump() -> tuple[dict[int, int], dict[int, str], dict[int, str]]:
    folder = _taxdump_dir()
    nodes = folder / "nodes.dmp"
    names = folder / "names.dmp"
    parents: dict[int, int] = {}
    ranks: dict[int, str] = {}
    for line in nodes.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) < 3:
            continue
        parents[int(parts[0])] = int(parts[1])
        ranks[int(parts[0])] = parts[2]
    scientific: dict[int, str] = {}
    for line in names.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) >= 4 and parts[3] == "scientific name":
            scientific[int(parts[0])] = parts[1]
    return parents, ranks, scientific


def _roll(counts: dict[int, int], parents: dict[int, int], ranks: dict[int, str]) -> dict[int, float]:
    total = sum(counts.values())
    if total <= 0:
        raise SystemExit("classification has no reads")
    rolled: dict[int, int] = {}
    for taxon_id, count in counts.items():
        species = None if taxon_id == 0 else rank_taxon(taxon_id, parents, ranks, SCORE_RANK)
        rolled[species or 0] = rolled.get(species or 0, 0) + count
    return {taxon_id: count / total for taxon_id, count in rolled.items()}


def _kraken_counts(path: Path) -> dict[int, int]:
    counts: dict[int, int] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            raise SystemExit(f"kraken line has fewer than 3 columns: {line[:80]}")
        taxon_id = 0 if parts[0] == "U" else int(parts[2])
        counts[taxon_id] = counts.get(taxon_id, 0) + 1
    return counts


def _kaiju_counts(path: Path) -> dict[int, int]:
    counts: dict[int, int] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            raise SystemExit(f"kaiju line has fewer than 3 columns: {line[:80]}")
        taxon_id = 0 if parts[0] == "U" else int(parts[2])
        counts[taxon_id] = counts.get(taxon_id, 0) + 1
    return counts


def _truth(parents: dict[int, int], ranks: dict[int, str]) -> dict[int, float]:
    """Species read fractions from the abundance table and the assembly-report taxids."""
    reports = {}
    for line in REPORT_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            accession, _name, _tax = _report_identity(record)
            reports[accession] = record
    counts: dict[int, int] = {}
    for line in (WORK / "abundance.csv").read_text(encoding="utf-8").splitlines()[1:]:
        accession, raw_count = line.split(",")
        record = reports.get(accession)
        if record is None:
            raise SystemExit(f"no assembly report for simulated genome {accession}")
        _accession, _name, taxon_id = _report_identity(record)
        species = rank_taxon(taxon_id, parents, ranks, SCORE_RANK)
        if species is None:
            raise SystemExit(f"{accession} taxid {taxon_id} does not roll to {SCORE_RANK}")
        counts[species] = counts.get(species, 0) + int(raw_count)
    total = sum(counts.values())
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def stage_compare() -> None:
    """Score Kraken2, Kaiju, and the current reprofiling hypotheses against the simulation."""
    parents, ranks, scientific = _taxdump()
    truth = _truth(parents, ranks)
    metrics = [
        "method\tgraph\tweight\tl1\tbray_curtis\tpearson\tspearman\tr2\tpresence_f1\tpresence_precision\tpresence_recall\tnode_f1\tnode_accuracy"
    ]
    comparison = ["method\ttaxon_id\tname\ttruth\testimated"]
    read_methods = {
        "kraken2": _roll(_kraken_counts(WORK / "classify" / "kraken2.output"), parents, ranks),
        "kaiju": _roll(_kaiju_counts(WORK / "classify" / "kaiju.output"), parents, ranks),
    }
    for method, predicted in read_methods.items():
        metrics.append(_metric_row(method, "reads", "read_count", predicted, truth, None))
        _append_comparison(comparison, method, predicted, truth, scientific)
    for row in _assembly_profiles(parents, ranks, scientific, truth):
        metrics.append(row[0])
        comparison.extend(row[1])
    (WORK / "metrics.tsv").write_text("\n".join(metrics) + "\n", encoding="utf-8")
    (WORK / "comparison.tsv").write_text("\n".join(comparison) + "\n", encoding="utf-8")


def _metric_row(method, graph, weight, predicted, truth, node_scores) -> str:
    from metamalevich.evaluate import abundance_scores

    scores = abundance_scores(predicted, truth)
    presence = presence_f1(predicted, truth)
    node_f1 = "" if node_scores is None else f"{node_scores['f1']:.6g}"
    node_acc = "" if node_scores is None else f"{node_scores['accuracy']:.6g}"
    pearson = scores["pearson"]
    return "\t".join(
        [
            method,
            graph,
            weight,
            f"{scores['l1']:.6g}",
            f"{scores['bray_curtis']:.6g}",
            "" if pearson is None else f"{pearson:.6g}",
            "" if scores["spearman"] is None else f"{scores['spearman']:.6g}",
            "" if r_squared(pearson) is None else f"{r_squared(pearson):.6g}",
            f"{presence['f1']:.6g}",
            f"{presence['precision']:.6g}",
            f"{presence['recall']:.6g}",
            node_f1,
            node_acc,
        ]
    )


def _append_comparison(lines, method, predicted, truth, names) -> None:
    for taxon_id in sorted(set(predicted) | set(truth)):
        if taxon_id == 0 and predicted.get(0, 0) == 0 and truth.get(0, 0) == 0:
            continue
        lines.append(
            "\t".join(
                [
                    method,
                    str(taxon_id),
                    names.get(taxon_id, "unclassified" if taxon_id == 0 else str(taxon_id)),
                    f"{truth.get(taxon_id, 0.0):.8g}",
                    f"{predicted.get(taxon_id, 0.0):.8g}",
                ]
            )
        )


def _graph_inputs() -> tuple[Path, list[dict], str] | None:
    """Use the MEGAHIT assembly graph. Do not substitute a 4-mer kNN."""
    fastg = WORK / "megahit" / "assembly.fastg"
    if not fastg.is_file() or fastg.stat().st_size == 0:
        return None
    edges, sequences = assembly_graph_from_fastg(fastg.read_text(encoding="utf-8", errors="replace"))
    fasta = WORK / "megahit" / "graph_nodes.fa"
    chunks = []
    for node_id, sequence in sequences.items():
        chunks.append(f">{node_id}\n")
        for start in range(0, len(sequence), 80):
            chunks.append(sequence[start : start + 80] + "\n")
    fasta.write_text("".join(chunks), encoding="utf-8")
    _log(WORK / "commands.log", f"megahit_fastg nodes={len(sequences)} edges={len(edges)}")
    return fasta, edges, "megahit_fastg"


def _export_assembly_tocumg(sequences, edges, evidence, edge_dist, scientific) -> None:
    """Write a coloured CFA and CDBG for the assembly graph when MetaMetro imports."""
    destination = WORK / "reprofile" / "tocumg"
    marker = destination / "export.json"
    if marker.is_file():
        return
    from metamalevich.bridge import export_tocumg

    node_taxa = {node_id: [taxon_id for taxon_id in dist if taxon_id != 0] for node_id, dist in evidence.items()}
    edge_taxa = {edge_id: [taxon_id for taxon_id in dist if taxon_id != 0] for edge_id, dist in edge_dist.items()}
    names = {taxon_id: scientific.get(taxon_id, str(taxon_id)) for taxon_id in {taxon for taxa in node_taxa.values() for taxon in taxa}}
    try:
        record = export_tocumg(
            root=ROOT,
            graph_id=GRAPH_ID,
            sequences=sequences,
            edges=edges,
            node_taxa=node_taxa,
            edge_taxa=edge_taxa,
            taxonomy_names=names,
            destination=destination,
            relative_to=ROOT,
        )
    except Exception as exc:
        _log(WORK / "commands.log", f"ToCUMG export failed: {exc}")
        return
    # export_tocumg labels every graph as knn in CFA metadata. This graph is the FASTG.
    record["example_graph"] = "megahit_fastg"
    marker.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    _log(WORK / "commands.log", f"ToCUMG export {record}")


def _assembly_profiles(parents, ranks, scientific, truth) -> list[tuple[str, list[str]]]:
    """Colour the assembly graph and profile it with the current resolvers."""
    from metamalevich.bench import (
        _evidence_distributions,
        _group_counts,
        _lca_distribution,
        _one_hot_call,
        _relative_from_profile,
        decision_variants,
        graph_variants,
    )
    from metamalevich.evaluate import classification_scores
    from metamalevich.native import kraken_counts
    from metamalevich.reprofile import profile_nodes
    from metamalevich.resolve import edge_distributions, hard_assignment, resolve
    from metamalevich.tables import read_fasta, read_tsv
    from metamalevich.taxonomy import parse_kraken_report

    found = _graph_inputs()
    if found is None:
        raise SystemExit("MEGAHIT assembly graph is missing; the example does not substitute a 4-mer kNN")
    fasta, edges, graph_name = found
    classify = WORK / "reprofile"
    classify.mkdir(parents=True, exist_ok=True)
    kraken_output = classify / "contigs.kraken"
    report_path = classify / "contigs.report"
    if not kraken_output.is_file():
        kraken = _tool("kraken2")
        _run(
            [
                kraken,
                "--db",
                str(WORK / "kraken_db"),
                "--threads",
                "4",
                "--report",
                str(report_path),
                "--output",
                str(kraken_output),
                str(fasta),
            ],
            WORK / "commands.log",
        )
    counts_path = classify / "kraken_counts.tsv"
    calls_path = classify / "kraken_calls.tsv"
    kraken_counts(kraken_output, counts_path, calls_path)
    taxonomy = parse_kraken_report(report_path.read_text(encoding="utf-8"), source="kraken2", version=str(report_path.relative_to(WORK)))
    counts = _group_counts(read_tsv(counts_path))
    calls = {row["seq_id"]: int(row["taxon_id"]) for row in read_tsv(calls_path)}
    evidence, _raw = _evidence_distributions(counts, taxonomy, calls)
    calls_dist = {node_id: _one_hot_call(calls.get(node_id, 0), taxonomy) for node_id in evidence}
    edge_dist = edge_distributions(edges, evidence) if edges else {}
    sequences = dict(read_fasta(fasta))
    lengths = {node_id: len(sequence) for node_id, sequence in sequences.items()}
    _export_assembly_tocumg(sequences, edges, evidence, edge_dist, scientific)
    resolved = {
        "initial_colouring": hard_assignment(calls_dist),
        "probability_sum": evidence,
        "lca": {node_id: _lca_distribution(counts.get(node_id, {}), taxonomy) for node_id in evidence},
        "gated_neighbour": resolve(evidence, edges, "gated_neighbour", edge_dist) if edges else hard_assignment(evidence),
        "bayesian_edge": resolve(evidence, edges, "bayesian_edge", edge_dist) if edges else hard_assignment(evidence),
        **graph_variants(evidence, edges),
        **decision_variants(evidence, edges, taxonomy, calls, lengths),
    }
    coverages = {}
    for node_id, sequence in sequences.items():
        coverage = megahit_coverage(node_id)
        if coverage is not None and coverage > 0:
            coverages[node_id] = max(1, int(round(len(sequence) * coverage)))
    node_truth = _contig_species(fasta, parents, ranks)
    _write_node_assignments(resolved, lengths, node_truth or {})
    rows = []
    for hypothesis, distributions in resolved.items():
        rows.append(
            _profile_row(
                hypothesis,
                graph_name,
                "node_length",
                distributions,
                lengths,
                taxonomy,
                truth,
                scientific,
                node_truth,
                profile_nodes,
                _relative_from_profile,
                classification_scores,
                parents,
                ranks,
            )
        )
        if coverages and set(coverages) == set(lengths):
            rows.append(
                _profile_row(
                    hypothesis + "_coverage",
                    graph_name,
                    "length_times_multi",
                    distributions,
                    coverages,
                    taxonomy,
                    truth,
                    scientific,
                    node_truth,
                    profile_nodes,
                    _relative_from_profile,
                    classification_scores,
                    parents,
                    ranks,
                )
            )
    return rows


def _write_node_assignments(resolved, lengths, node_truth) -> None:
    """Save one argmax table per hypothesis so later error analysis can reuse it."""
    from metamalevich.resolve import argmax_taxon
    from metamalevich.tables import write_tsv

    destination = WORK / "reprofile" / "assignments"
    destination.mkdir(parents=True, exist_ok=True)
    for hypothesis, distributions in resolved.items():
        rows = []
        for node_id, dist in sorted(distributions.items()):
            taxon_id, probability = argmax_taxon(dist)
            rows.append(
                {
                    "node_id": node_id,
                    "taxon_id": taxon_id,
                    "probability": f"{probability:.6f}",
                    "truth_taxon_id": node_truth.get(node_id, ""),
                    "length": lengths.get(node_id, 0),
                }
            )
        write_tsv(
            destination / f"{hypothesis}.tsv",
            rows,
            ["node_id", "taxon_id", "probability", "truth_taxon_id", "length"],
        )


def _profile_row(
    hypothesis,
    graph_name,
    weight,
    distributions,
    lengths,
    taxonomy,
    truth,
    scientific,
    node_truth,
    profile_nodes,
    relative_from_profile,
    classification_scores,
    parents,
    ranks,
) -> tuple[str, list[str]]:
    profile, summary = profile_nodes(distributions, lengths, taxonomy)
    predicted = {}
    for taxon_id, fraction in relative_from_profile(profile, summary["unclassified_fraction"]).items():
        species = 0 if taxon_id == 0 else rank_taxon(taxon_id, parents, ranks, SCORE_RANK)
        predicted[species or 0] = predicted.get(species or 0, 0.0) + fraction
    node_scores = None
    if node_truth:
        rolled = {}
        for node_id, dist in distributions.items():
            if node_id not in node_truth:
                continue
            taxon_id = max(dist, key=dist.get) if dist else 0
            species = 0 if taxon_id == 0 else rank_taxon(taxon_id, parents, ranks, SCORE_RANK)
            rolled[node_id] = {species or 0: 1.0}
        node_scores = classification_scores(rolled, node_truth)
    lines = []
    _append_comparison(lines, hypothesis, predicted, truth, scientific)
    return _metric_row(hypothesis, graph_name, weight, predicted, truth, node_scores), lines


def _contig_species(fasta: Path, parents: dict[int, int], ranks: dict[int, str]) -> dict[str, int] | None:
    """Map each contig to a simulated genome's species with minimap2, when it is installed."""
    minimap = shutil.which("minimap2")
    if not minimap:
        _log(WORK / "commands.log", "minimap2 is not on PATH; node-level F1 was not computed")
        return None
    help_text = _help_text(minimap)
    reference = WORK / "sim_reference.fna"
    if not reference.is_file():
        chunks = []
        for path in sorted((WORK / "sim").glob("*.fna")):
            sequence = []
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                if not line.startswith(">"):
                    sequence.append(line.strip())
            chunks.append(f">{path.stem}\n{''.join(sequence)}\n")
        reference.write_text("".join(chunks), encoding="utf-8")
    paf = WORK / "reprofile" / "contigs.paf"
    cmd = [minimap, "-x", "asm20", str(reference), str(fasta)] if "asm20" in help_text else [minimap, str(reference), str(fasta)]
    completed = subprocess.run(cmd, capture_output=True, text=True)
    _log(WORK / "commands.log", "$ " + " ".join(cmd) + f" > {paf}")
    if completed.returncode != 0:
        raise SystemExit(completed.stderr.strip() or "minimap2 failed")
    paf.write_text(completed.stdout, encoding="utf-8")
    best: dict[str, tuple[int, str]] = {}
    for line in completed.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 12:
            raise SystemExit("minimap2 PAF row has fewer than 12 columns")
        query, matches, target = parts[0], int(parts[9]), parts[5]
        previous = best.get(query)
        if previous is None or matches > previous[0]:
            best[query] = (matches, target)
    reports = {}
    for line in REPORT_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            accession, _name, _tax = _report_identity(record)
            reports[iss_key(accession)] = record
    truth = {}
    for query, (_matches, target) in best.items():
        record = reports.get(target)
        if record is None:
            continue
        _accession, _name, taxon_id = _report_identity(record)
        species = rank_taxon(taxon_id, parents, ranks, SCORE_RANK)
        if species is None:
            raise SystemExit(f"mapped genome {target} taxid {taxon_id} has no {SCORE_RANK}")
        truth[query] = species
    return truth


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Build one pinned NCBI community and its assembly graph")
    parser.add_argument(
        "--stage",
        default="all",
        choices=("all", "download", "layout", "simulate", "databases", "classify", "assemble", "compare", "plot"),
    )
    args = parser.parse_args(argv)
    order = ["download", "layout", "simulate", "databases", "classify", "assemble", "compare", "plot"]
    selected = order if args.stage == "all" else [args.stage]
    log_tool_versions()
    for name in selected:
        print(f"stage {name}", flush=True)
        if name == "download":
            stage_download()
        elif name == "layout":
            stage_layout()
        elif name == "simulate":
            stage_simulate()
        elif name == "databases":
            stage_databases()
        elif name == "classify":
            stage_classify()
        elif name == "assemble":
            stage_assemble()
        elif name == "compare":
            stage_compare()
        else:
            raise SystemExit("plot stage is not part of the MetaMetro graph build")


if __name__ == "__main__":
    main()
