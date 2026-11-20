"""Optional checks that need extra software."""

from __future__ import annotations

from pathlib import Path

import pytest

pytestmark = pytest.mark.optional


def test_pyg_adapter_when_installed() -> None:
    """Build a PyG Data object when torch_geometric is installed."""
    pytest.importorskip("torch_geometric")
    from metametro.converters.cgt_to_pyg import to_pyg
    from metametro.fixtures import mock_cgt

    data = to_pyg(mock_cgt())
    assert data.x.shape[0] == 4
    assert data.edge_index.shape == (2, 5)
    assert data.y.shape == (4,)


def test_strain_matches_bubbleblower_when_present() -> None:
    """Node sequences and taxon colours match BubbleBlower's strain bubble."""
    sibling = Path(__file__).resolve().parents[2] / "BubbleBlower" / "src"
    if not (sibling / "bubbleblower" / "fixtures.py").is_file():
        pytest.skip("BubbleBlower checkout is not beside metametro")
    import sys

    sys.path.insert(0, str(sibling))
    from bubbleblower.fixtures import strain_bubble

    graph = strain_bubble()
    sequences = sorted((unitig.sequence, tuple(unitig.color_ids)) for unitig in graph.cdbg.unitigs)
    from metametro.bench.data.universal.bubbles import strain_bubble as our_strain
    from metametro.converters.cfa_to_cdbg import cfa_to_cdbg

    ours = cfa_to_cdbg(our_strain())
    our_sequences = sorted((unitig.sequence, tuple(unitig.color_ids)) for unitig in ours.unitigs)
    assert our_sequences == sequences
