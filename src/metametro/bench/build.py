"""Build one benchmark into a contract and, when tools allow, three graphs."""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from metametro.bench.colourings import BuildContext
from metametro.bench.data.universal.catalog import (
    classifier_pin,
    require_programs,
    write_community_contract,
    write_external_contract,
)
from metametro.bench.data.universal.community import highest_megahit_contigs, load_pairs
from metametro.bench.identity import contract_identity, write_identity
from metametro.bench.materialize import Materialized, materialize
from metametro.bench.paths import default_outdir, repo_root
from metametro.bench.registry import PREPARE, resolve
from metametro.bench.spec import BenchSpec
from metametro.contracts.assembly import fastg_to_cfa
from metametro.errors import ContractError
from metametro.tables import read_yaml, write_yaml


@dataclass(frozen=True)
class BuildResult:
    """One finished ``benchbuild`` run."""

    spec: BenchSpec
    outdir: Path
    identity: str
    colourings: tuple[str, ...]
    status: str


@dataclass(frozen=True)
class ClassifierInputs:
    """Kraken2 and Kaiju paths given to ``build``. ``None`` means not given."""

    kraken_db: Path | None = None
    kraken_output: Path | None = None
    kaiju_db: Path | None = None
    kaiju_output: Path | None = None

    def given(self) -> dict[str, Path]:
        """Return the paths that were given, by ``BuildContext`` field name."""
        return {name: Path(value) for name, value in vars(self).items() if value is not None}

    def check(self) -> None:
        """Raise ``ContractError`` for every given path that does not exist."""
        missing = [
            f"{name} (--{name.replace('_', '-')}) does not exist: {path}"
            for name, path in self.given().items()
            if not path.exists()
        ]
        if missing:
            raise ContractError(missing)

    def context(self, colourings: tuple[str, ...] | None, **found: object) -> BuildContext:
        """Return a ``BuildContext``. A given path wins over one in ``found``."""
        return BuildContext(selected=colourings, **{**found, **self.given()})


def build(
    name: str,
    *,
    outdir: Path | None = None,
    colourings: tuple[str, ...] | None = None,
    root: Path | None = None,
    execute: bool = False,
    gtfs: Path | None = None,
    kraken_db: Path | None = None,
    kraken_output: Path | None = None,
    kaiju_db: Path | None = None,
    kaiju_output: Path | None = None,
) -> BuildResult:
    """Build ``name``.

    In-process benchmarks always write CFA, CDBG, and CGT. Community and
    external benchmarks always write their contract. They download and assemble
    only when ``execute`` is true, and they stop if a required program is
    missing. ``colourings`` restricts the registry. ``None`` applies every
    colouring that can run, including colourings registered later.
    ``kraken_db``, ``kraken_output``, ``kaiju_db``, and ``kaiju_output`` go
    to the colourings of every benchmark kind. A given path wins over one found
    under ``work/``. A given path that does not exist raises ``ContractError``
    before anything is written.
    """
    spec = resolve(name)
    classifiers = ClassifierInputs(kraken_db, kraken_output, kaiju_db, kaiju_output)
    classifiers.check()
    destination = Path(outdir) if outdir is not None else default_outdir(spec, root=root)
    if destination.exists() and any(destination.iterdir()):
        raise ContractError([f"output directory is not empty: {destination}"])
    if spec.kind == "inprocess":
        return _inprocess(spec, destination, colourings, classifiers, requested=name)
    if spec.kind == "community":
        return _community(
            spec, destination, colourings, classifiers, execute=execute, requested=name
        )
    if spec.kind == "external":
        return _external(
            spec, destination, colourings, classifiers, execute=execute, gtfs=gtfs, requested=name
        )
    raise ContractError([f"unknown benchmark kind {spec.kind}"])


def build_all(
    *,
    execute: bool = False,
    gtfs: Path | None = None,
    root: Path | None = None,
    kraken_db: Path | None = None,
    kraken_output: Path | None = None,
    kaiju_db: Path | None = None,
    kaiju_output: Path | None = None,
) -> list[BuildResult]:
    """Build every registered benchmark.

    In-process graphs are always materialised. Community and external benches
    write their contract. They download and assemble only when ``execute`` is
    true. A directory that already has ``manifest.yaml`` and ``identity.sha256``
    is left as it is. Classifier paths go to every benchmark that is built.
    """
    from metametro.bench.registry import list_specs

    results: list[BuildResult] = []
    for spec in list_specs():
        destination = default_outdir(spec, root=root)
        if (destination / "manifest.yaml").is_file() and (destination / "identity.sha256").is_file():
            from metametro.bench.ground_truth import check_gold_digest

            check_gold_digest(destination)
            results.append(
                BuildResult(
                    spec=spec,
                    outdir=destination,
                    identity=(destination / "identity.sha256").read_text(encoding="utf-8").split()[0],
                    colourings=(),
                    status="present",
                )
            )
            continue
        results.append(
            build(
                spec.name,
                execute=execute if spec.kind != "inprocess" else False,
                gtfs=gtfs,
                root=root,
                kraken_db=kraken_db,
                kraken_output=kraken_output,
                kaiju_db=kaiju_db,
                kaiju_output=kaiju_output,
            )
        )
    return results


def _finish_manifest(
    spec: BenchSpec,
    destination: Path,
    *,
    status: str,
    identity: str,
    colourings: tuple[str, ...],
    requested: str,
    tools: dict[str, str] | None = None,
) -> None:
    manifest = {
        "bench": spec.name,
        "assembler": spec.assembler,
        "properties": spec.properties,
        "status": status,
        "colourings": list(colourings),
        "identity": identity,
        "scoring": list(spec.scoring),
        "gold_namespaces": list(spec.gold_namespaces),
    }
    if requested != spec.name:
        manifest["requested_name"] = requested
    if tools is not None:
        manifest["tools"] = dict(tools)
    previous = destination / "manifest.yaml"
    if status == "built" and previous.is_file():
        written = read_yaml(previous)
        for key in ("features", "gold_sha256"):
            if written.get(key) is not None:
                manifest[key] = written[key]
    write_yaml(destination / "manifest.yaml", manifest)


def _inprocess(
    spec: BenchSpec,
    destination: Path,
    colourings: tuple[str, ...] | None,
    classifiers: ClassifierInputs,
    *,
    requested: str,
) -> BuildResult:
    ctx = classifiers.context(colourings)
    graph = PREPARE[spec.name](destination, ctx)
    made: Materialized = materialize(graph, destination, ctx, gold_namespaces=spec.gold_namespaces)
    _finish_manifest(
        spec,
        destination,
        status="built",
        identity=made.identity,
        colourings=made.colourings,
        requested=requested,
    )
    return BuildResult(
        spec=spec,
        outdir=destination,
        identity=made.identity,
        colourings=made.colourings,
        status="built",
    )


def _community(
    spec: BenchSpec,
    destination: Path,
    colourings: tuple[str, ...] | None,
    classifiers: ClassifierInputs,
    *,
    execute: bool,
    requested: str,
) -> BuildResult:
    write_community_contract(spec, destination)
    if not execute:
        digest = contract_identity(destination)
        write_identity(destination / "identity.sha256", digest)
        _finish_manifest(
            spec,
            destination,
            status="contract",
            identity=digest,
            colourings=(),
            requested=requested,
        )
        return BuildResult(
            spec=spec, outdir=destination, identity=digest, colourings=(), status="contract"
        )
    require_programs(spec, ("datasets", "samovar", "megahit", "megahit_toolkit"))
    graph = _assemble_community(spec, destination)
    from metametro.bench.data.universal import ncbi_pipeline as pipe

    tools = pipe.log_tool_versions()
    ctx = _community_context(spec, destination, colourings, classifiers)
    made = materialize(graph, destination, ctx, gold_namespaces=spec.gold_namespaces)
    _finish_manifest(
        spec,
        destination,
        status="built",
        identity=made.identity,
        colourings=made.colourings,
        requested=requested,
        tools=tools,
    )
    return BuildResult(
        spec=spec,
        outdir=destination,
        identity=made.identity,
        colourings=made.colourings,
        status="built",
    )


def _assemble_community(spec: BenchSpec, destination: Path):
    from metametro.bench.data.universal import ncbi_pipeline as pipe

    root = repo_root()
    pipe.configure(
        example=destination / "contract",
        work=destination / "work",
        root=root,
        graph_id=spec.name,
        score_rank=spec.score_rank,
        total_reads=spec.total_reads,
        report_path=assembly_report_path(spec, root),
        zip_name=f"ncbi_{spec.name}.zip",
    )
    if spec.parent:
        pipe.stage_download = lambda: _require_parent_fasta(spec, root)
    for stage in ("download", "layout", "simulate", "assemble"):
        getattr(pipe, f"stage_{stage}")()
    fastg = destination / "work" / "megahit" / "assembly.fastg"
    if not fastg.is_file() or fastg.stat().st_size == 0:
        raise ContractError([f"{spec.name} assembled no MEGAHIT FASTG"])
    k_dir = destination / "work" / "megahit" / "intermediate_contigs"
    names = [path.name for path in k_dir.iterdir()] if k_dir.is_dir() else []
    chosen = highest_megahit_contigs(names)
    if chosen is None:
        raise ContractError([f"{spec.name} has no MEGAHIT k-mer contig graph"])
    k_value = int(Path(chosen).name.split(".", 1)[0][1:])
    _install_classifier_library(spec, destination, root)
    _try_classifier_indexes(destination)
    return fastg_to_cfa(fastg, k=k_value, graph_id=spec.name)


def _install_classifier_library(spec: BenchSpec, destination: Path, root: Path) -> None:
    """Copy the full non-synonymous ``db`` FASTA set used to build classifiers."""
    rows = load_pairs(classifier_pin(spec) / "accessions.tsv")
    library = destination / "work" / "classifier_db"
    library.mkdir(parents=True, exist_ok=True)
    missing = []
    for row in rows:
        if row["role"] != "db":
            continue
        source = root / "data" / "raw" / "fasta" / f"{row['accession']}.fna"
        if not source.is_file() or source.stat().st_size == 0:
            missing.append(row["accession"])
            continue
        target = library / source.name
        if not target.is_file():
            shutil.copyfile(source, target)
    if missing:
        shown = ", ".join(missing[:8])
        if len(missing) > 8:
            shown += ", ..."
        raise ContractError([f"{spec.name} classifier library missing FASTA: {shown}"])


def _try_classifier_indexes(destination: Path) -> None:
    """Build Kraken2 and Kaiju indexes when those tools can run. Skip if they cannot."""
    from metametro.bench.data.universal import ncbi_pipeline as pipe

    kraken = shutil.which("kraken2")
    kaiju = shutil.which("kaiju")
    samovar = shutil.which("samovar")
    if kraken is None and kaiju is None:
        return
    if samovar is None:
        return
    try:
        pipe.stage_databases()
    except SystemExit:
        return


def _community_context(
    spec: BenchSpec,
    destination: Path,
    colourings: tuple[str, ...] | None,
    classifiers: ClassifierInputs,
) -> BuildContext:
    """Point colourings at classifier databases built from the full ``db`` set.

    A path given to ``build`` wins over the one found under ``work/``.
    """
    work = destination / "work"
    kraken_db = work / "kraken_db"
    kaiju_db = work / "kaiju_db"
    kraken_out = work / "classifier" / "contigs.kraken"
    kaiju_out = work / "classifier" / "contigs.kaiju"
    return classifiers.context(
        colourings,
        work_dir=work,
        kraken_db=kraken_db if (kraken_db / "hash.k2d").is_file() else None,
        kaiju_db=kaiju_db if kaiju_db.is_dir() else None,
        kraken_output=kraken_out if kraken_out.is_file() else None,
        kaiju_output=kaiju_out if kaiju_out.is_file() else None,
    )


def assembly_report_path(spec: BenchSpec, root: Path) -> Path:
    """Return the NCBI assembly report of the community that downloaded the FASTA.

    A community with a ``parent`` reuses the parent's FASTA and never runs
    ``datasets``, so it reads the report of its first ancestor without a parent.
    """
    source = spec
    while source.parent:
        source = resolve(source.parent)
    return root / "data" / "raw" / f"{source.name}_assembly_data_report.jsonl"


def _require_parent_fasta(spec: BenchSpec, root: Path) -> None:
    rows = load_pairs(root / "src" / "metametro" / "bench" / "data" / "pins" / spec.name / "accessions.tsv")
    missing = []
    for row in rows:
        path = root / "data" / "raw" / "fasta" / f"{row['accession']}.fna"
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(row["accession"])
    if missing:
        shown = ", ".join(missing[:8])
        if len(missing) > 8:
            shown += ", ..."
        raise ContractError([f"{spec.name} reuses FASTA from {spec.parent}; missing {shown}"])
    report = assembly_report_path(spec, root)
    if not report.is_file():
        raise ContractError([f"{spec.name} needs {report.name}; build {spec.parent} first"])


def _external(
    spec: BenchSpec,
    destination: Path,
    colourings: tuple[str, ...] | None,
    classifiers: ClassifierInputs,
    *,
    execute: bool,
    gtfs: Path | None,
    requested: str,
) -> BuildResult:
    if execute and spec.name.startswith("phage_species_5"):
        from metametro.bench.genomes import check_genomes

        check_genomes(spec.name, repo_root() / "data" / "raw" / "genomes")
        require_programs(spec, ("samovar", "megahit", "megahit_toolkit"))
    write_external_contract(spec, destination)
    if not execute:
        digest = contract_identity(destination)
        write_identity(destination / "identity.sha256", digest)
        _finish_manifest(
            spec,
            destination,
            status="contract",
            identity=digest,
            colourings=(),
            requested=requested,
        )
        return BuildResult(
            spec=spec, outdir=destination, identity=digest, colourings=(), status="contract"
        )
    if spec.name.startswith("phage_species_5"):
        genomes = repo_root() / "data" / "raw" / "genomes"
        script = repo_root() / "scripts" / "phage_baseline.py"
        completed = subprocess.run(
            [
                sys.executable,
                str(script),
                "--genomes",
                str(genomes),
                "--work",
                str(destination / "work"),
                "--k",
                str(spec.k),
                "--total-reads",
                str(spec.total_reads),
                "--samovar",
                shutil.which("samovar") or "samovar",
                "--megahit",
                shutil.which("megahit") or "megahit",
            ],
            check=False,
        )
        if completed.returncode != 0:
            raise ContractError([f"phage baseline exited {completed.returncode}"])
        from metametro.formats.cfa.io import load_cfa

        coloured = destination / "work" / "coloured_cfa"
        if not coloured.is_dir():
            raise ContractError(["phage baseline did not write coloured_cfa"])
        from metametro.bench.phage_gold import (
            paint_phage_gold,
            phage_membership,
            write_taxon_membership,
            write_taxon_values,
        )

        graph = load_cfa(coloured)
        membership = phage_membership(graph, spec.name, genomes, spec.k)
        graph = paint_phage_gold(graph, spec.name, genomes, spec.k, membership=membership)
        found: dict[str, Path] = {}
        wants_kraken = colourings is None or "kraken2" in colourings
        if wants_kraken and classifiers.kraken_db is None and classifiers.kraken_output is None:
            from metametro.bench.kraken_pin import checked_kraken_pin

            pinned = checked_kraken_pin(spec.name, graph, destination / "work" / "classifier" / "nodes.fna")
            if pinned is not None:
                found["kraken_output"] = pinned
        ctx = classifiers.context(colourings, **found)
        from metametro.bench.bubble_gold import write_phage_bubbles

        def write_phage_ground_truth(outdir: Path) -> None:
            write_taxon_values(outdir, spec.name)
            write_taxon_membership(outdir, membership)
            write_phage_bubbles(outdir)

        made = materialize(
            graph,
            destination,
            ctx,
            gold_namespaces=spec.gold_namespaces,
            before_identity=write_phage_ground_truth,
        )
        _finish_manifest(
            spec,
            destination,
            status="built",
            identity=made.identity,
            colourings=made.colourings,
            requested=requested,
        )
        return BuildResult(
            spec=spec,
            outdir=destination,
            identity=made.identity,
            colourings=made.colourings,
            status="built",
        )
    if spec.name == "roxel":
        require_programs(spec, ("Rscript",))
        from metametro.bench.data.universal.roxel import build_roxel_cfa

        graph = build_roxel_cfa(k=spec.k or 21)
        ctx = classifiers.context(colourings)
        made = materialize(graph, destination, ctx, gold_namespaces=spec.gold_namespaces)
        _finish_manifest(
            spec,
            destination,
            status="built",
            identity=made.identity,
            colourings=made.colourings,
            requested=requested,
        )
        return BuildResult(
            spec=spec,
            outdir=destination,
            identity=made.identity,
            colourings=made.colourings,
            status="built",
        )
    if spec.name == "spb_ground_transit":
        if gtfs is None or not Path(gtfs).is_file():
            raise ContractError(["spb_ground_transit needs an existing GTFS zip passed as --gtfs"])
        script = repo_root() / "examples" / "spb_transit" / "build_cfa.py"
        staged = destination / "cfa_source"
        completed = subprocess.run(
            [sys.executable, str(script), "--gtfs", str(gtfs), "--out", str(staged), "--k", "21"],
            check=False,
        )
        if completed.returncode != 0:
            raise ContractError([f"GTFS CFA build exited {completed.returncode}"])
        from metametro.formats.cfa.io import load_cfa

        graph = load_cfa(staged)
        ctx = classifiers.context(colourings)
        made = materialize(graph, destination, ctx, gold_namespaces=spec.gold_namespaces)
        _finish_manifest(
            spec,
            destination,
            status="built",
            identity=made.identity,
            colourings=made.colourings,
            requested=requested,
        )
        return BuildResult(
            spec=spec,
            outdir=destination,
            identity=made.identity,
            colourings=made.colourings,
            status="built",
        )
    raise ContractError([f"no executor for {spec.name}"])
