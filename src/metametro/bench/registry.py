"""Catalog of benchmarks. Strong100 and other prebuilt bundles are not listed."""

from __future__ import annotations

from metametro.bench.data.universal.catalog import community_specs, external_specs
from metametro.bench.data.universal.inprocess import (
    prepare_bubble_error_1,
    prepare_bubble_nested_2,
    prepare_bubble_reads_2,
    prepare_bubble_shared_2,
    prepare_bubble_strain_2,
    prepare_bubble_strain_3_n50,
    prepare_mock_bubble,
    prepare_mock_chain,
    prepare_synthetic_reads_2,
)
from metametro.bench.spec import BenchSpec
from metametro.formats.cfa.model import CfaGraph
from metametro.bench.colourings import BuildContext
from pathlib import Path
from typing import Callable

Prepare = Callable[[Path, BuildContext], CfaGraph]

SPECS: dict[str, BenchSpec] = {
    "bubble_strain_2": BenchSpec(
        name="bubble_strain_2",
        canonical_identity="a6a805ca1435f8afbf4646107d21b1ecc7f6a433a75bda22e785ade851ef456c",
        canonical_gold="dae6f67ef1b912cde17eeef945e8e6708034277698fc2791e480701ea632319e",
        assembler="inprocess",
        properties="taxa2",
        summary="Two taxa, one bubble. BubbleBlower toy. Debubbler truth: retain.",
        scoring=("assembly",),
        kind="inprocess",
        legacy_names=("strain_bubble", "bubbleblower_toy"),
    ),
    "bubble_error_1": BenchSpec(
        name="bubble_error_1",
        canonical_identity="c4232c8cc17dc3c48b52a6f0b006eb841b83e9c2fbe09d31812a87cbee50af3a",
        canonical_gold="bbf0a479aa9bcb003de9d7a0b3226c8801035d1d56766acd67ddbcf4be03e4d5",
        assembler="inprocess",
        properties="taxon1",
        summary="One taxon, one error bubble. Debubbler truth: pop.",
        scoring=("assembly",),
        kind="inprocess",
        legacy_names=("error_bubble",),
    ),
    "bubble_nested_2": BenchSpec(
        name="bubble_nested_2",
        canonical_identity="3a8a773fa3ed27252bfb2f50d66736b9e19ee532d95832af523b415f2e8ae2d7",
        canonical_gold="a8b689d34e3049b695f2e6cfc6a8421a44e1707041dd89cbc177c9fe064dce81",
        assembler="inprocess",
        properties="taxa2",
        summary="Two taxa, one outer bubble and one inner bubble.",
        scoring=("assembly",),
        kind="inprocess",
        legacy_names=("nested_bubbles",),
    ),
    "bubble_shared_2": BenchSpec(
        name="bubble_shared_2",
        canonical_identity="254486559f42d4c7a4c5e16b7e3ef1eaf9c99e816450be9ae959862fc6e62969",
        canonical_gold="4c4e599efa9b2b1a1987c259ff72d4dfa5a860e5752ce7a276a411fc4083ce28",
        assembler="inprocess",
        properties="taxa2",
        summary="Two taxa sharing one node. No bubble.",
        scoring=("assembly",),
        kind="inprocess",
        legacy_names=("shared_duplicate_node",),
    ),
    "bubble_strain_3_n50": BenchSpec(
        name="bubble_strain_3_n50",
        canonical_identity="d9a1db5c2ad772545b39c120cb5906623d2bf5e532d1df1885ae03c3e58c835f",
        canonical_gold="a4385d9c9141c9c2f6ef7acec71a4872d3032a1f2e038473874815c1c792d4bb",
        assembler="inprocess",
        properties="seed42",
        summary="Three strains and 50 disjoint bubbles (25 strain, 25 error), seed 42.",
        scoring=("assembly",),
        kind="inprocess",
        legacy_names=("bubble_benchmark",),
    ),
    "bubble_reads_2": BenchSpec(
        name="bubble_reads_2",
        canonical_identity="949fa8a0578a4396d9df20af4c449b4d3c44ed4a7d9e41859d2bc1e3e15b2a41",
        canonical_gold="7cb08e179afae3c16653ae4d0b9a621dbd8280fe3f24046be22c9fc9837cf823",
        assembler="inprocess",
        properties="k5_len12_seed1",
        summary="Two strains, read length 12, de Bruijn k=5, 20 reads per sample.",
        scoring=("assembly",),
        kind="inprocess",
        gold_namespaces=("taxon", "genome", "accession"),
        legacy_names=("bubbleblower_reads",),
    ),
    "mock_bubble": BenchSpec(
        name="mock_bubble",
        canonical_identity="cd7a2e150dbb344035066030794be8ef5bef66fbf8e20942d8eccdeaeeced095",
        canonical_gold="bbe5f7acb070efa714d4890e92c7fb1aa51316c45142f12dd07f5f84490ab849",
        assembler="inprocess",
        properties="k3",
        summary="MetaMetro four-node bubble fixture.",
        scoring=("assembly",),
        kind="inprocess",
        gold_namespaces=("taxon", "genome"),
    ),
    "mock_chain": BenchSpec(
        name="mock_chain",
        canonical_identity="a397fd80a0c2779bc4362da5f79d17ace86527e4bb0b3aa7260c2600ee905377",
        canonical_gold="3bb81d4e5600db9238aeb4594301102f48da6eb903dd5bbbb65644754b61d89f",
        assembler="inprocess",
        properties="k3",
        summary="MetaMetro six-node chain that compacts to three unitigs.",
        scoring=("assembly",),
        kind="inprocess",
        gold_namespaces=("taxon", "genome"),
    ),
    "synthetic_reads_2": BenchSpec(
        name="synthetic_reads_2",
        canonical_identity="89c2cfbd147852473e31f1e6c46b92e2f41ac15b2d940dc05a47d85c8d245448",
        canonical_gold="078cf246083cc5ca56930dcdb2f718599217f0c10f7a6986b2433581d0be5e66",
        assembler="inprocess",
        properties="k8_len16_seed1",
        summary="Two synthetic genomes and an in-process de Bruijn graph at k=8.",
        scoring=("assembly", "profiling"),
        kind="inprocess",
        gold_namespaces=("taxon", "genome", "accession"),
    ),
}

SPECS.update(community_specs())
SPECS.update(external_specs())

PREPARE: dict[str, Prepare] = {
    "bubble_strain_2": prepare_bubble_strain_2,
    "bubble_error_1": prepare_bubble_error_1,
    "bubble_nested_2": prepare_bubble_nested_2,
    "bubble_shared_2": prepare_bubble_shared_2,
    "bubble_strain_3_n50": prepare_bubble_strain_3_n50,
    "bubble_reads_2": prepare_bubble_reads_2,
    "mock_bubble": prepare_mock_bubble,
    "mock_chain": prepare_mock_chain,
    "synthetic_reads_2": prepare_synthetic_reads_2,
}

_LEGACY = {
    alias: spec.name
    for spec in SPECS.values()
    for alias in spec.legacy_names
}


def resolve(name: str) -> BenchSpec:
    """Return the spec for a canonical name or a legacy alias."""
    from metametro.errors import ContractError

    canonical = _LEGACY.get(name, name)
    spec = SPECS.get(canonical)
    if spec is None:
        known = ", ".join(sorted(SPECS))
        raise ContractError([f"unknown benchmark {name!r}. Known benchmarks: {known}"])
    return spec


def canonical_identity(name: str) -> str | None:
    """Return the expected ``identity.sha256`` of a default build of ``name``, or ``None``.

    A default build applies every colouring that can run (``colourings=None``)
    and is given no classifier paths. ``None`` means the identity depends on
    inputs that are not pinned (a GTFS feed, the Roxel data, NCBI downloads).
    """
    return resolve(name).canonical_identity or None


def canonical_gold(name: str) -> str | None:
    """Return the expected ``gold_sha256`` (``manifest.yaml``) of a default build of ``name``, or ``None``.

    It is the digest of ``gold.sha256``, which covers every file under
    ``ground_truth/``; ``None`` when the bench has no ``canonical_identity``.
    """
    return resolve(name).canonical_gold or None


def list_specs() -> tuple[BenchSpec, ...]:
    """Return specs in a stable order."""
    return tuple(SPECS[name] for name in sorted(SPECS))
