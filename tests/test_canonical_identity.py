"""Mandatory: canonical identities match the docs and the in-process builds. No network.

phage_10 needs MEGAHIT and Samovar, so no mandatory test rebuilds it: the CI job
``phage10-solve`` (``.github/workflows/full-tests.yml``) builds it and checks its
identity and gold digest, then runs ``tests/test_phage10_built.py`` on that build.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from metametro.bench import build, canonical_gold, canonical_identity, list_specs
from metametro.tables import read_yaml

pytestmark = pytest.mark.mandatory

DOCS = Path(__file__).resolve().parents[1] / "docs" / "benchmarks.md"


def _documented() -> dict[str, str]:
    section = DOCS.read_text(encoding="utf-8").split("### Canonical identities", 1)[1].split("\n### ", 1)[0]
    rows = re.findall(r"^\| `([a-z0-9_]+)` \| `([0-9a-f]{64})` \| `([0-9a-f]{64})` \|$", section, flags=re.MULTILINE)
    return {name: (identity, gold) for name, identity, gold in rows}


def test_stored_identities_match_the_docs() -> None:
    stored = {
        spec.name: (canonical_identity(spec.name), canonical_gold(spec.name))
        for spec in list_specs()
        if canonical_identity(spec.name)
    }
    assert stored == _documented()
    assert all(gold for _identity, gold in stored.values())
    assert canonical_gold("phage_10") == "c80f418585a7082a2f0fb99278566e7b30e77b86f57466d98e69fdf6085b7c44"
    assert canonical_identity("phage_10") == canonical_identity("phage_species_5_x10")
    assert canonical_identity("phage_10") == "78121c1b05a5814bfc17ac3b7650ed5223a0f584cfada24e9224b7ba0364564c"
    assert canonical_identity("roxel") is None and canonical_identity("spb_ground_transit") is None
    assert all(canonical_identity(spec.name) for spec in list_specs() if spec.kind == "inprocess")


@pytest.mark.parametrize("name", [spec.name for spec in list_specs() if spec.kind == "inprocess"])
def test_default_build_has_the_canonical_identity_and_gold(tmp_path: Path, name: str) -> None:
    assert build(name, outdir=tmp_path / name).identity == canonical_identity(name)
    assert read_yaml(tmp_path / name / "manifest.yaml")["gold_sha256"] == canonical_gold(name)
