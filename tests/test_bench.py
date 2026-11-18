"""Mandatory: benchbuild writes three graphs and a stable identity."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metametro.bench import build
from metametro.bench.colourings import BuildContext, register, unregister
from metametro.bench.leakage import assert_no_target_leak
from metametro.bench.registry import list_specs
from metametro.bench.scoring.assembly.bubble_strain_2 import expected_action
from metametro.bench.scoring.assembly.universal.lengths import n50
from metametro.bench.scoring.assembly_binning.universal.contig_f1 import contig_f1
from metametro.bench.scoring.profiling.universal.abundance import l1, presence_f1
from metametro.cli import main
from metametro.contracts.colouring import paint_namespace
from metametro.errors import ContractError
from metametro.formats.cfa.io import load_cfa
from metametro.formats.cdbg.io import load_cdbg
from metametro.formats.cgt.io import load_cgt

pytestmark = pytest.mark.mandatory


def test_registry_lists_inprocess_benches() -> None:
    """The in-process benchmarks are registered under their new names."""
    names = {spec.name for spec in list_specs()}
    assert "bubble_strain_2" in names
    assert "bubble_reads_2" in names
    assert "synthetic_reads_2" in names
    assert "mock_bubble" in names
    assert "strong100" not in names


def test_benchbuild_strain_is_stable_and_coloured(tmp_path: Path) -> None:
    """Two builds of the two-taxon bubble share an identity and write three graphs."""
    first = build("bubble_strain_2", outdir=tmp_path / "a")
    second = build("strain_bubble", outdir=tmp_path / "b")
    assert first.identity == second.identity
    assert "as_built" in first.colourings
    assert "composition_kmeans" in first.colourings
    assert "decaying" in first.colourings
    for folder in (first.outdir, second.outdir):
        assert (folder / "cfa" / "metadata.yaml").is_file()
        assert (folder / "cdbg" / "metadata.yaml").is_file()
        assert (folder / "cgt" / "metadata.yaml").is_file()
        assert expected_action(folder) == "retain"
    cfa = load_cfa(first.outdir / "cfa")
    assert {row["node_id"] for row in cfa.nodes} == {"S", "A", "B", "T"}
    assert cfa.sequences["A"] == "ATATCG"
    assert_no_target_leak(cfa)
    text = (first.outdir / "cfa" / "nodes.tsv").read_text(encoding="utf-8")
    assert "expected_action" not in text
    cdbg = load_cdbg(first.outdir / "cdbg")
    cgt = load_cgt(first.outdir / "cgt")
    assert cgt.num_nodes == len(cdbg.unitigs)


def test_new_colouring_is_applied(tmp_path: Path) -> None:
    """A colouring registered after the benchmark still runs."""

    def mark(graph, _ctx: BuildContext):
        return paint_namespace(
            graph,
            {row["node_id"]: ["unit"] for row in graph.nodes},
            namespace="mark",
            operation="merge",
        )

    register("unit_test_mark", mark)
    try:
        result = build("bubble_error_1", outdir=tmp_path / "out")
    finally:
        unregister("unit_test_mark")
    assert "unit_test_mark" in result.colourings
    cfa = load_cfa(result.outdir / "cfa")
    assert any(row["namespace"] == "mark" for row in cfa.colors or [])


def test_target_column_is_rejected() -> None:
    """A truth column on the CFA is not a legal benchmark graph."""
    from metametro.bench.data.universal.bubbles import strain_bubble

    graph = strain_bubble()
    graph.metadata["features"]["node"]["truth_taxon_id"] = "int"
    graph.node_header = [*graph.node_header, "truth_taxon_id"]
    with pytest.raises(ContractError, match="target leaked"):
        assert_no_target_leak(graph)


def test_reads_benchmark_keeps_genome_ids_out_of_the_graph(tmp_path: Path) -> None:
    """Read colouring uses samples. Genome ids stay in ground truth."""
    result = build("bubble_reads_2", outdir=tmp_path / "reads")
    assert "read_depth" in result.colourings
    assert "read_accession" not in result.colourings
    cfa = load_cfa(result.outdir / "cfa")
    namespaces = {row["namespace"] for row in cfa.colors or []}
    assert "sample" in namespaces
    assert "truth" not in namespaces
    truth = (result.outdir / "ground_truth" / "read_to_genome.tsv").read_text(encoding="utf-8")
    assert "strain_A" in truth
    assert "strain_A" not in (result.outdir / "cfa" / "colors.tsv").read_text(encoding="utf-8")


@pytest.mark.parametrize("name", ["bubble_reads_2", "synthetic_reads_2"])
def test_reads_benchmark_has_node_coverage(tmp_path: Path, name: str) -> None:
    """Read benchmarks carry k-mer read coverage on every CFA node."""
    result = build(name, outdir=tmp_path / name)
    cfa = load_cfa(result.outdir / "cfa")
    assert "coverage" in cfa.node_header
    k = cfa.metadata["k"]
    reads = [line.strip() for line in (result.outdir / "reads" / "reads.fastq").read_text().splitlines()[1::4]]
    total = sum(max(len(read) - k + 1, 0) for read in reads)
    values = [float(row["coverage"]) for row in cfa.nodes]
    assert min(values) >= 1.0
    assert sum(values) == total


def test_synthetic_fifty_is_deterministic(tmp_path: Path) -> None:
    """The 50-bubble generator does not depend on the output path."""
    left = build("bubble_strain_3_n50", outdir=tmp_path / "left")
    right = build("bubble_strain_3_n50", outdir=tmp_path / "right")
    assert left.identity == right.identity
    assert left.n_nodes if hasattr(left, "n_nodes") else load_cfa(left.outdir / "cfa")
    cfa = load_cfa(left.outdir / "cfa")
    assert len(cfa.nodes) > 50
    truth = (left.outdir / "ground_truth" / "bubbles.tsv").read_text(encoding="utf-8").strip().splitlines()
    assert len(truth) == 51
    assert sum(line.endswith("\tpop") for line in truth) == 25


def test_cli_benchbuild_list(capsys: pytest.CaptureFixture[str]) -> None:
    """``metametro benchbuild --list`` prints the two-taxon bubble."""
    assert main(["benchbuild", "--list"]) == 0
    assert "bubble_strain_2" in capsys.readouterr().out


def test_cli_list_colourings(capsys: pytest.CaptureFixture[str]) -> None:
    """``metametro benchbuild --list-colourings`` names the migrated colourings."""
    assert main(["benchbuild", "--list-colourings"]) == 0
    text = capsys.readouterr().out
    assert "composition_kmeans" in text
    assert "kraken2" in text
    assert "kaiju" in text
    assert "decaying" in text


def test_default_outdir_shape(tmp_path: Path) -> None:
    """The default tree is data/bench/{name}/{assembler}/{properties}."""
    result = build("mock_bubble", root=tmp_path)
    assert result.outdir == tmp_path / "data" / "bench" / "mock_bubble" / "inprocess" / "k3"
    assert (result.outdir / "cgt" / "indptr.npy").is_file()


def test_scoring_helpers() -> None:
    """Universal scorers match closed-form cases."""
    assert n50([10, 20, 30]) == 30
    assert l1({1: 0.5, 2: 0.5}, {1: 1.0, 2: 1.0}) == 0.0
    assert presence_f1({1: 1.0}, {1: 1.0, 2: 1.0}) == pytest.approx(2 / 3)
    perfect = contig_f1(np.array([0, 0, 1, 1]), np.array([5, 5, 9, 9]))
    assert perfect == pytest.approx(1.0)


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


def test_community_contract_is_stable_and_renamed(tmp_path: Path) -> None:
    """Pinned communities keep their accession tables and lognormal abundances."""
    from metametro.bench.data.universal.community import every_other_strain, half_strains, load_pairs

    left = build("4domain_family_100", outdir=tmp_path / "a", execute=False)
    right = build("high100", outdir=tmp_path / "b", execute=False)
    assert left.status == "contract"
    assert left.identity == right.identity
    assert left.spec.pair_count == 100
    assert "strong100" not in {spec.name for spec in list_specs()}
    text = (left.outdir / "ground_truth" / "abundance.csv").read_text(encoding="utf-8").strip().splitlines()
    assert text[0] == "taxid,N_sample"
    total = sum(int(line.split(",")[1]) for line in text[1:])
    assert total == 100_000
    pins = Path(__file__).resolve().parents[1] / "src" / "metametro" / "bench" / "data" / "pins"
    held = load_pairs(pins / "bacteria_species_20_heldout" / "accessions.tsv")
    half = load_pairs(pins / "bacteria_strain_10" / "accessions.tsv")
    assert half == half_strains(held)
    genus = load_pairs(pins / "3domain_genus_75" / "accessions.tsv")
    genus_half = load_pairs(pins / "3domain_genus_75_half" / "accessions.tsv")
    assert genus_half == every_other_strain(genus)
    phage = build("phage_10", outdir=tmp_path / "phage", execute=False)
    assert phage.spec.name == "phage_species_5_x10"
    assert phage.spec.pair_count == 5
    assert phage.spec.total_reads == 4000


def test_benchbuild_all_builds_inprocess_and_contracts(tmp_path: Path) -> None:
    """``--all`` materialises in-process graphs and writes contracts without downloading."""
    from metametro.bench.build import build_all

    results = build_all(root=tmp_path, execute=False)
    by_name = {item.spec.name: item for item in results}
    assert by_name["bubble_strain_2"].status == "built"
    assert (by_name["bubble_strain_2"].outdir / "cgt" / "metadata.yaml").is_file()
    assert by_name["4domain_family_100"].status == "contract"
    assert "strong100" not in by_name
    again = build_all(root=tmp_path, execute=False)
    assert all(item.status == "present" for item in again)
    from metametro.bench.check import check_registry

    checked = dict(check_registry(tmp_path))
    assert checked["bubble_strain_2"] == "built"
    assert checked["bubble_reads_2"] == "built"
    assert checked["4domain_family_100"] == "contract"
    assert checked["phage_species_5"] == "contract"
    assert checked["roxel"] == "contract"
    assert checked["spb_ground_transit"] == "contract"


def test_roxel_execute_needs_rscript(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Roxel does not invent a street graph when R is absent."""
    monkeypatch.setattr("metametro.bench.data.universal.catalog.shutil.which", lambda _name: None)
    with pytest.raises(ContractError, match="Rscript"):
        build("roxel", outdir=tmp_path / "roxel", execute=True)


def test_community_execute_stops_without_tools(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A community build does not invent a graph when the assembler is absent."""
    monkeypatch.setattr("metametro.bench.data.universal.catalog.shutil.which", lambda _name: None)
    with pytest.raises(ContractError, match="needs these programs"):
        build("bacteria_species_20_heldout", outdir=tmp_path / "held", execute=True)



def test_download_validation_writes_only_under_the_bench_work_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The download manifest, checksums, and reports stay out of tracked repo files."""
    import json

    from metametro.bench.data.universal import ncbi_pipeline as pipe

    work = tmp_path / "data" / "bench" / "toy" / "work"
    monkeypatch.setattr(pipe, "ROOT", tmp_path)
    monkeypatch.setattr(pipe, "WORK", work)
    monkeypatch.setattr(pipe, "GRAPH_ID", "toy_community")
    monkeypatch.setattr(pipe, "EXAMPLE", tmp_path / "data" / "bench" / "toy" / "contract")
    fasta_dir = tmp_path / "data" / "raw" / "fasta"
    fasta_dir.mkdir(parents=True)
    rows = []
    records = []
    for role, accession, sequence in (("sim", "GCF_000000001.1", "ACGTACGTAA"), ("db", "GCF_000000002.1", "TTGCATGCAA")):
        (fasta_dir / f"{accession}.fna").write_text(f">{accession}\n{sequence}\n", encoding="utf-8")
        rows.append({"accession": accession, "genus": "Escherichia", "pair_id": "p1", "role": role})
        records.append({"accession": accession, "organism": {"organismName": "Escherichia coli", "taxId": 562}})
    report = tmp_path / "data" / "raw" / "report.jsonl"
    report.write_text("\n".join(json.dumps(record) for record in records) + "\n", encoding="utf-8")
    zip_path = tmp_path / "data" / "raw" / "ncbi.zip"
    zip_path.write_bytes(b"zip")
    pipe._validate_genomes(rows, fasta_dir, report, zip_path, "datasets 0")
    written = sorted(path.name for path in (work / "download").iterdir())
    assert written == ["acquisition_report.md", "checksums.txt", "dataset-audit.md", "download_manifest.tsv"]
    assert not (tmp_path / "data" / "manifests").exists()
    assert not (tmp_path / "data" / "checksums").exists()
    assert not (tmp_path / "docs").exists()
    audit = (work / "download" / "dataset-audit.md").read_text(encoding="utf-8")
    acquisition = (work / "download" / "acquisition_report.md").read_text(encoding="utf-8")
    assert "toy_community" in audit and "toy_community" in acquisition
    assert "data/bench/toy/contract/accessions.tsv" in acquisition
    assert "1 sim / 1 db, 1 pairs" in audit
    assert "| 2 files |" in audit
    for stale in ("heldout_genera", "20 sim", "Shigella", "40 files"):
        assert stale not in audit + acquisition


def test_derived_community_reads_the_downloading_ancestor_report(tmp_path: Path) -> None:
    """Half and x10 communities reuse the assembly report of the community that ran datasets."""
    from metametro.bench.build import _require_parent_fasta, assembly_report_path
    from metametro.bench.data.universal.catalog import pin_dir
    from metametro.bench.data.universal.community import load_pairs
    from metametro.bench.registry import resolve

    expected = tmp_path / "data" / "raw" / "4domain_family_100_assembly_data_report.jsonl"
    for name in ("4domain_family_100", "4domain_family_100_half", "4domain_family_100_half_x10"):
        assert assembly_report_path(resolve(name), tmp_path) == expected
    spec = resolve("4domain_family_100_half")
    pin = tmp_path / "src" / "metametro" / "bench" / "data" / "pins" / spec.name
    pin.mkdir(parents=True)
    (pin / "accessions.tsv").write_bytes((pin_dir(spec.name) / "accessions.tsv").read_bytes())
    fasta = tmp_path / "data" / "raw" / "fasta"
    fasta.mkdir(parents=True)
    for row in load_pairs(pin_dir(spec.name) / "accessions.tsv"):
        (fasta / f"{row['accession']}.fna").write_text(">x\nACGT\n", encoding="utf-8")
    with pytest.raises(ContractError, match="build 4domain_family_100 first"):
        _require_parent_fasta(spec, tmp_path)
    expected.write_text("{}\n", encoding="utf-8")
    _require_parent_fasta(spec, tmp_path)


@pytest.mark.parametrize("name", ["bubble_reads_2", "mock_chain"])
def test_bench_cgt_has_a_feature_layer_without_new_colours(tmp_path: Path, name: str) -> None:
    """The CGT carries k=4 composition plus coverage features; colours are untouched."""
    from metametro.bench import load_bench_cgt

    result = build(name, outdir=tmp_path / name)
    stored = load_cgt(result.outdir / "cgt")
    names = stored.metadata["node_feature_names"]
    assert len([item for item in names if item.startswith("kmer:k4:")]) == 136
    assert names[-1] == "coverage"
    assert stored.node_features.shape == (stored.num_nodes, 137)
    assert np.all(np.isfinite(stored.node_features))
    loaded = load_bench_cgt(result.outdir)
    plain = load_bench_cgt(result.outdir, features=False)
    assert np.array_equal(loaded.node_features, stored.node_features)
    assert plain.node_features.shape == (stored.num_nodes, 0)
    assert np.array_equal(plain.node_colors, stored.node_colors)
    assert np.array_equal(loaded.node_colors, stored.node_colors)
    assert stored.node_colors.size and stored.node_colors.any()


def test_tool_versions_are_logged_and_absent_tools_are_named(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``--version`` output goes to commands.log; a tool that is not on PATH is ``absent``."""
    from metametro.bench.build import _finish_manifest
    from metametro.bench.data.universal import ncbi_pipeline as pipe
    from metametro.bench.registry import resolve
    from metametro.tables import read_yaml

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "megahit"
    fake.write_text("#!/bin/sh\necho 'MEGAHIT v1.2.9' >&2\n", encoding="utf-8")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_dir))
    monkeypatch.setattr(pipe, "WORK", tmp_path / "work")
    versions = pipe.log_tool_versions()
    assert versions == {
        "megahit": "MEGAHIT v1.2.9",
        "kraken2": "absent",
        "kaiju": "absent",
        "minimap2": "absent",
        "samovar": "absent",
    }
    log = (tmp_path / "work" / "commands.log").read_text(encoding="utf-8")
    assert "tool version megahit: MEGAHIT v1.2.9" in log
    assert "tool version kaiju: absent" in log
    spec = resolve("bacteria_species_20_heldout")
    _finish_manifest(spec, tmp_path / "bench", status="built", identity="x", colourings=(), requested=spec.name, tools=versions)
    assert read_yaml(tmp_path / "bench" / "manifest.yaml")["tools"] == versions
