"""Mandatory: canonical identities match the docs and the in-process builds. No network."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from metametro.bench import build, canonical_identity, list_specs

pytestmark = pytest.mark.mandatory

DOCS = Path(__file__).resolve().parents[1] / "docs" / "benchmarks.md"


def _documented() -> dict[str, str]:
    section = DOCS.read_text(encoding="utf-8").split("### Canonical identities", 1)[1].split("\n### ", 1)[0]
    return dict(re.findall(r"^\| `([a-z0-9_]+)` \| `([0-9a-f]{64})` \|$", section, flags=re.MULTILINE))


def test_stored_identities_match_the_docs() -> None:
    stored = {spec.name: canonical_identity(spec.name) for spec in list_specs() if canonical_identity(spec.name)}
    assert stored == _documented()
    assert canonical_identity("phage_10") == canonical_identity("phage_species_5_x10")
    assert canonical_identity("phage_10") == "78121c1b05a5814bfc17ac3b7650ed5223a0f584cfada24e9224b7ba0364564c"
    assert canonical_identity("roxel") is None and canonical_identity("spb_ground_transit") is None
    assert all(canonical_identity(spec.name) for spec in list_specs() if spec.kind == "inprocess")


@pytest.mark.parametrize("name", [spec.name for spec in list_specs() if spec.kind == "inprocess"])
def test_default_build_has_the_canonical_identity(tmp_path: Path, name: str) -> None:
    assert build(name, outdir=tmp_path / name).identity == canonical_identity(name)
