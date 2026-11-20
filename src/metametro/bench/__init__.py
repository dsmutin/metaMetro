"""Benchmark construction and scoring."""

from metametro.bench.build import BuildResult, build
from metametro.bench.ground_truth import (
    GoldBubble,
    GoldLabels,
    load_bench_gold,
    load_bench_gold_bubbles,
    load_bench_gold_cgt,
)
from metametro.bench.registry import list_specs, resolve
from metametro.bench.select import load_bench_cdbg, load_bench_cfa, load_bench_cgt, namespaces_for

__all__ = [
    "BuildResult",
    "GoldBubble",
    "GoldLabels",
    "build",
    "list_specs",
    "load_bench_cdbg",
    "load_bench_cfa",
    "load_bench_cgt",
    "load_bench_gold",
    "load_bench_gold_bubbles",
    "load_bench_gold_cgt",
    "namespaces_for",
    "resolve",
]
