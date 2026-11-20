"""``metametro benchfetch`` command."""

from __future__ import annotations

import argparse
from pathlib import Path

from metametro.bench.genomes import fetch_genomes, pinned_genomes
from metametro.bench.paths import repo_root


def main(argv: list[str] | None = None) -> int:
    """Download the pinned genomes of one benchmark from NCBI."""
    parser = argparse.ArgumentParser(
        prog="metametro benchfetch",
        description=(
            "Download the pinned reference genomes of a benchmark from NCBI E-utilities "
            "and check each against its pinned SHA-256."
        ),
    )
    parser.add_argument("name", help="benchmark name, for example phage_species_5 or phage_10")
    parser.add_argument(
        "--genome-dir",
        type=Path,
        default=None,
        help="where to write <name>.fna (default: data/raw/genomes in the checkout)",
    )
    parser.add_argument("--force", action="store_true", help="replace files whose SHA-256 differs from the pin")
    parser.add_argument("--list", action="store_true", help="print the pinned table and exit")
    args = parser.parse_args(argv)
    pinned = pinned_genomes(args.name)
    if args.list:
        print("name\taccession\ttaxid\torganism\tlength")
        for genome in pinned:
            print(f"{genome.name}\t{genome.accession}\t{genome.taxid}\t{genome.organism}\t{genome.length}")
        return 0
    genome_dir = args.genome_dir if args.genome_dir is not None else repo_root() / "data" / "raw" / "genomes"
    status = fetch_genomes(args.name, genome_dir, force=args.force)
    by_name = {genome.name: genome for genome in pinned}
    for name, state in status.items():
        print(f"{name}\t{by_name[name].accession}\t{state}\t{genome_dir / (name + '.fna')}")
    return 0
