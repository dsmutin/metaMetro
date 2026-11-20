"""Benchmark names, assemblers, and the properties that enter the output path."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BenchSpec:
    """One benchmark the ``benchbuild`` command can materialise.

    ``properties`` is the last path segment under
    ``data/bench/{name}/{assembler}/{properties}``. It records the knobs that
    change the graph (read budget, k, seed) and is not a machine path.
    ``gold_namespaces`` are the evaluation targets of this bench: colours,
    columns, or labels in these namespaces are moved to ``ground_truth/`` and
    are never readable as colours. ``canonical_identity`` is the
    ``identity.sha256`` of a default build (``colourings=None``, no classifier
    paths) when it is reproducible from pinned inputs, else ``""``;
    ``canonical_gold`` is that build's ``gold_sha256``.
    """

    name: str
    assembler: str
    properties: str
    summary: str
    scoring: tuple[str, ...]
    kind: str
    legacy_names: tuple[str, ...] = ()
    score_rank: str = ""
    total_reads: int = 0
    pair_count: int = 0
    domain_count: int = 0
    k: int = 0
    parent: str = ""
    gold_namespaces: tuple[str, ...] = ("taxon",)
    canonical_identity: str = ""
    canonical_gold: str = ""

    def path_parts(self) -> tuple[str, str, str]:
        """Return ``(name, assembler, properties)`` for the default output tree."""
        return (self.name, self.assembler, self.properties)
