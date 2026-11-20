"""Mandatory when a built phage_10 is available: its gold class counts. Skipped otherwise.

The CI job ``phage10-solve`` builds phage_10 and runs this file with
``METAMETRO_PHAGE10_DIR`` set; locally it reads the default benchbuild folder.
"""

from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import pytest

from metametro.bench.paths import default_outdir
from metametro.bench.registry import resolve
from metametro.tables import read_tsv

pytestmark = pytest.mark.mandatory


def _bench() -> Path:
    given = os.environ.get("METAMETRO_PHAGE10_DIR")
    folder = Path(given) if given else default_outdir(resolve("phage_10"))
    if not (folder / "ground_truth" / "bubbles.tsv").is_file():
        pytest.skip(f"no built phage_10 at {folder}; the phage10-solve CI job builds one")
    return folder


def test_phage_10_gold_class_counts() -> None:
    folder = _bench()
    _, bubbles = read_tsv(folder / "ground_truth" / "bubbles.tsv")
    assert Counter(row["class"] for row in bubbles) == {"inter_species": 92, "intra_species": 1, "error": 1}
    _, nodes = read_tsv(folder / "ground_truth" / "taxon_nodes.tsv")
    assert sum(1 for row in nodes if row["taxon"]) == 903
    header, membership = read_tsv(folder / "ground_truth" / "taxon_membership.tsv")
    genomes = [column for column in header[3:] if not column.endswith("_fraction")]
    counts = Counter(sum(row[name] == "1" for name in genomes) for row in membership)
    assert counts == {0: 12, 1: 854, 2: 117}
