"""Output locations for a benchmark build. Results are gitignored."""

from __future__ import annotations

from pathlib import Path

from metametro.bench.spec import BenchSpec


def checkout_root(package_file: Path) -> Path:
    """Return the checkout that contains ``package_file``, or raise ``ContractError``.

    An installed package (for example the conda package in site-packages) has no
    checkout above it, so ``data/bench`` has no default location there.
    """
    root = Path(package_file).resolve().parents[3]
    if (root / "pyproject.toml").is_file() and (root / "src" / "metametro").is_dir():
        return root
    from metametro.errors import ContractError

    raise ContractError(
        [
            "metametro is not running from a checkout, so data/bench has no default location; "
            "pass --outdir (or outdir=/root= in Python)"
        ]
    )


def repo_root() -> Path:
    """Return the MetaMetro checkout that contains this package."""
    return checkout_root(Path(__file__))


def default_outdir(spec: BenchSpec, *, root: Path | None = None) -> Path:
    """Return ``data/bench/{name}/{assembler}/{properties}`` under ``root``.

    ``root`` defaults to the repository. Callers that pass ``--outdir`` replace
    this path entirely; they do not append another segment.
    """
    base = repo_root() if root is None else Path(root)
    name, assembler, properties = spec.path_parts()
    return base / "data" / "bench" / name / assembler / properties
