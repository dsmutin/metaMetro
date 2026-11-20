"""Optional: the pinned phage genomes download from NCBI and match their pin. Needs network."""

from __future__ import annotations

from pathlib import Path

import pytest

from metametro.bench.genomes import check_genomes, fetch_genomes

pytestmark = pytest.mark.optional


def test_fetch_from_ncbi_matches_the_pin(tmp_path: Path) -> None:
    """The real NCBI records reproduce every pinned hash."""
    status = fetch_genomes("phage_species_5", tmp_path)
    assert set(status.values()) == {"fetched"}
    check_genomes("phage_species_5", tmp_path)
