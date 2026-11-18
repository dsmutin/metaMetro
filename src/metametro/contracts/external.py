"""External baselines: Samovar ISS (Contract 1) and MEGAHIT (Contract 2)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def samovar_generate_command(
    genome_dir: str | Path,
    output_dir: str | Path,
    *,
    n_samples: int = 2,
    total_reads: int = 200,
    host_fraction: float = 0,
    seed: int = 1,
    samovar: str = "samovar",
) -> list[str]:
    """Argv for ``samovar generate``, which drives InSilicoSeq."""
    return [
        samovar,
        "generate",
        "--genome_dir",
        str(genome_dir),
        "--output_dir",
        str(output_dir),
        "--n_samples",
        str(n_samples),
        "--total_reads",
        str(total_reads),
        "--host_fraction",
        str(host_fraction),
        "--seed",
        str(seed),
    ]


def megahit_command(
    reads: str | Path,
    output_dir: str | Path,
    *,
    k: int = 21,
    threads: int = 4,
    megahit: str = "megahit",
    min_count: int = 2,
) -> list[str]:
    """Argv for a single-k MEGAHIT assembly. ``k`` must be odd and at least 15."""
    if k < 15 or k % 2 == 0:
        raise ValueError("MEGAHIT k must be an odd integer >= 15")
    return [
        megahit,
        "-r",
        str(reads),
        "-o",
        str(output_dir),
        "--k-min",
        str(k),
        "--k-max",
        str(k),
        "--k-step",
        "2",
        "--min-count",
        str(min_count),
        "--keep-tmp-files",
        "-t",
        str(threads),
    ]


FLYE_READ_MODES = ("nano-raw", "nano-hq")
MINIMAP2_PRESETS = ("asm5", "asm10", "asm20", "map-ont", "map-hifi", "map-pb", "sr")


def flye_command(
    reads: str | Path,
    output_dir: str | Path,
    *,
    mode: str = "nano-raw",
    meta: bool = True,
    threads: int = 4,
    flye: str = "flye",
) -> list[str]:
    """Argv for a Flye assembly of ONT reads. ``mode`` is ``nano-raw`` or ``nano-hq``."""
    if mode not in FLYE_READ_MODES:
        raise ValueError(f"Flye mode must be one of {', '.join(FLYE_READ_MODES)}")
    if threads < 1:
        raise ValueError("Flye threads must be >= 1")
    argv = [flye, f"--{mode}", str(reads), "--out-dir", str(output_dir), "--threads", str(threads)]
    if meta:
        argv.append("--meta")
    return argv


def minimap2_command(
    target: str | Path,
    query: str | Path,
    *,
    preset: str = "asm5",
    cigar: bool = True,
    threads: int = 4,
    minimap2: str = "minimap2",
) -> list[str]:
    """Argv for minimap2 that writes PAF to stdout; pass it to ``run_command_to_file``.

    ``cigar=True`` adds ``-c`` so each PAF line carries base-level alignment and
    identity can be read from columns 10 and 11.
    """
    if preset not in MINIMAP2_PRESETS:
        raise ValueError(f"minimap2 preset must be one of {', '.join(MINIMAP2_PRESETS)}")
    if threads < 1:
        raise ValueError("minimap2 threads must be >= 1")
    argv = [minimap2, "-x", preset]
    if cigar:
        argv.append("-c")
    return argv + ["-t", str(threads), str(target), str(query)]


def contig2fastg_command(
    contigs: str | Path,
    k: int,
    *,
    toolkit: str = "megahit_toolkit",
) -> list[str]:
    """Argv for ``megahit_toolkit contig2fastg``. FASTG is written to stdout."""
    if k < 15 or k % 2 == 0:
        raise ValueError("MEGAHIT k must be an odd integer >= 15")
    return [toolkit, "contig2fastg", str(k), str(contigs)]


def run_samovar_generate(
    genome_dir: str | Path,
    output_dir: str | Path,
    *,
    n_samples: int = 2,
    total_reads: int = 200,
    host_fraction: float = 0,
    seed: int = 1,
    samovar: str = "samovar",
) -> None:
    """Write the Samovar ISS config and execute the generated Snakemake script.

    ``samovar generate`` only writes ``.generate/generate.sh``. The reads
    appear after that script runs.
    """
    run_command(
        samovar_generate_command(
            genome_dir,
            output_dir,
            n_samples=n_samples,
            total_reads=total_reads,
            host_fraction=host_fraction,
            seed=seed,
            samovar=samovar,
        )
    )
    script = Path(output_dir) / ".generate" / "generate.sh"
    if not script.is_file():
        raise FileNotFoundError(script)
    run_command(["bash", str(script), "--directory", str(output_dir)])


def run_command(argv: list[str]) -> None:
    """Run a baseline command. Fail if the executable is missing or it exits non-zero."""
    if shutil.which(argv[0]) is None and not Path(argv[0]).is_file():
        raise FileNotFoundError(argv[0])
    subprocess.run(argv, check=True)


def run_command_to_file(argv: list[str], destination: str | Path) -> None:
    """Run a baseline command and write its stdout to ``destination``."""
    if shutil.which(argv[0]) is None and not Path(argv[0]).is_file():
        raise FileNotFoundError(argv[0])
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            subprocess.run(argv, check=True, stdout=handle)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    temporary.replace(path)
