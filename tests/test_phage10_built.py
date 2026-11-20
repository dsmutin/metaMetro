"""Mandatory when a built phage_10 is available: its gold class counts. Skipped when there is none.

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
from metametro.tables import read_tsv, read_yaml

pytestmark = pytest.mark.mandatory


def _bench() -> Path:
    """The built phage_10 folder; skip only when no build was given or found.

    A folder named by ``METAMETRO_PHAGE10_DIR``, or a default folder holding a
    built bench (``status: built`` in ``manifest.yaml``), fails the test when
    ``bubbles.tsv`` is missing.
    """
    given = os.environ.get("METAMETRO_PHAGE10_DIR")
    folder = Path(given) if given else default_outdir(resolve("phage_10"))
    manifest = folder / "manifest.yaml"
    if not given and (not manifest.is_file() or (read_yaml(manifest) or {}).get("status") != "built"):
        pytest.skip(f"no built phage_10 at {folder}; the phage10-solve CI job builds one")
    if not (folder / "ground_truth" / "bubbles.tsv").is_file():
        pytest.fail(f"{folder} has no ground_truth/bubbles.tsv; rebuild phage_10")
    return folder


def test_a_given_phage_10_without_bubbles_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("METAMETRO_PHAGE10_DIR", str(tmp_path))
    with pytest.raises(pytest.fail.Exception, match="no ground_truth/bubbles.tsv"):
        _bench()


def test_phage_10_gold_class_counts() -> None:
    folder = _bench()
    _, bubbles = read_tsv(folder / "ground_truth" / "bubbles.tsv")
    assert Counter(row["class"] for row in bubbles) == {"inter_species": 92, "error": 1, "unresolved": 1}
    _, nodes = read_tsv(folder / "ground_truth" / "taxon_nodes.tsv")
    assert sum(1 for row in nodes if row["taxon"]) == 903
    header, membership = read_tsv(folder / "ground_truth" / "taxon_membership.tsv")
    genomes = [column for column in header[3:] if not column.endswith("_fraction")]
    counts = Counter(sum(row[name] == "1" for name in genomes) for row in membership)
    assert counts == {0: 12, 1: 854, 2: 117}
