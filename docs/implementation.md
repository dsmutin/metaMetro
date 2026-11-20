# Implementation

The public entry points are:

```python
cfa = load_cfa(path)
validate_cfa(cfa)
cdbg = cfa_to_cdbg(cfa)
validate_cdbg(cdbg)
cgt = cdbg_to_cgt(
    cdbg,
    node_features=...,
    edge_features=...,
    node_labels=...,
    node_annotation=[("namespace", "feature")],
    edge_annotation=[("namespace", "feature")],
)
cdbg = transfer_annotations(cdbg, cfa, namespace=..., provenance=...)
cdbg = aggregate_annotations(cdbg, namespace=..., feature=..., policy="weighted_mean", provenance=...)
validate_cgt(cgt)
lineage = node_lineage(cgt, dense_id)       # source_id, cfa_node_ids
edge_ids = cgt_edge_cfa_ids(cdbg)           # CSR slot → CFA edge id
sequence = resolve_sequence(cgt, cdbg, dense_id)  # unitig DNA, not a CGT column
external = cgt_from_csr(...)  # external CSR; colours stay off X
csc = csc_from_cgt(cgt)          # incoming index; does not reorder X_edge
graph = to_pyg(cgt)              # metametro.to_pyg; requires torch_geometric
result = run_ds(cgt, seed=0)     # does not modify cgt
predictions = predictions_from_ds(result)  # node predictions; confidence is the softmax probability
edited = apply_edit_proposal(proposal, cdbg)  # new CDBG; cdbg is not modified
edge = edge_prediction(cgt, csr_slot, predicted_class=..., probability=..., confidence=..., model_id=..., model_version=..., cdbg=cdbg)
```

Fixtures in `metametro.fixtures`: `mock_cfa()`, `mock_cdbg()`, `mock_cgt()`, `chain_cfa()`, `chain_cdbg()`, `synthetic_genomes()`. All but `synthetic_genomes()` are also exported from `metametro`.

The complete list of public names is `metametro.__all__`. Besides the calls above it exports:

- reading and writing: `load_cdbg`, `load_cgt`, `dump_cfa`, `dump_cdbg`, `dump_cgt`;
- conversion back and edge arrays: `cdbg_to_cfa`, `edge_index_array`;
- colouring (Contract 4): `colour_cfa`, `colour_by_reads`, `colour_by_composition`, `colour_decaying`, `paint_namespace`, `filter_colours`, `kmeans_onehot`;
- CDBG annotations: `annotate_cdbg`, `get_node_annotations`, `get_edge_annotations`;
- identity checks: `assert_cgt_matches_cdbg`, `link_cfa_edge_id`;
- edits and predictions: `EditProposal`, `GraphEdit`, `validate_edit_proposal`, `validate_edit_application`, `GraphPrediction`;
- contract versions: `CONTRACTS`.

`metametro.bench` exports `build`, `BuildResult`, `list_specs`, `resolve`, `load_bench_cfa`, `load_bench_cdbg`, `load_bench_cgt`, and `namespaces_for`.

| Contract | Module | What the minimal code does |
| --- | --- | --- |
| 1 | `contracts/assembly.py` `simulate_metagenome` | Deterministic slices of the source genomes. External path: `samovar generate`, then the generated `.generate/generate.sh` Snakemake script, which runs InSilicoSeq. |
| 2 | `dbg_from_sequences`, `contigs_to_dbg`, `fastg_to_cfa`, `gfa_to_cfa` | k-mer de Bruijn graph, or a GFA repeat graph (`graph_type: repeat`) with per-link overlap. External commands: MEGAHIT and Flye. The T-phage baseline loads MEGAHIT's intermediate FASTG so branches and bubbles survive; `strand_junction_counts` and `directed_bubble_sources` report that structure. |
| 3 | `formats/cfa` | Directory loader, dumper, and validator. The de Bruijn builder already emits CFA. |
| 4 | `contracts/colouring.py` | Explicit set operations and the depth / density rule. k-mer graphs are counted in one pass over the reads; the counts match the general substring rule. |
| 5 | `converters/cfa_to_cdbg.py` | ToCUMG chain compaction described in `docs/formats.md`. `graph_type` is unchanged. |
| 6 | `converters/cdbg_to_cgt.py` | CSR, dense ids, aligned features and colours, optional feature registry. `csc_from_cgt` derives incoming adjacency without storing it. `resolve_sequence` reads unitig DNA from the CDBG. |
| 7 | `contracts/ds.py` | NumPy convolution by default; PyG `GCNConv` when requested. Node predictions carry an explicit confidence equal to the softmax class probability. Edge predictions are a separate constructor and are not emitted here. |

Schema files live next to each format (`formats/*/schema.yaml`) and pin schema 1.0. Validators raise `ContractError` with every violation they collected. The ONT + metaFlye repeat-graph run is in [repeat-graph-run.md](repeat-graph-run.md).

## What this tree deliberately does not do

- CFA does not import Bifrost or PyTorch.
- CDBG topology does not embed coverage, GC, entropy, k-mer composition, or feature matrices. Those values may be stored in the optional columnar annotation sidecar (`formats/cdbg/annotations.py`) or passed into `cdbg_to_cgt` as arrays. The sidecar is not `X_node` / `X_edge`. Compaction does not fill it.
- CGT construction does not allocate an `N×N` adjacency matrix. CSC is derived from CSR on demand and is not a second canonical topology.
- CGT rows do not store DNA. `resolve_sequence` fails when the unitig or the dense id is missing.
- Training labels are not copied into `X_node` or `X_edge`. Predictions and confidence stay off those arrays.
- The PyG adapter does not become a fourth on-disk format.
- The NumPy trainer is the tested Contract 7 implementation. PyG is an optional backend, not a second result schema.
- Edit proposals do not mutate the input CDBG or CFA, do not bisect a CFA node, and do not invent a biological id for a new unitig. There is no correction model.

## Tests

Mandatory pytest covers, for each format: a valid object, a missing field, a bad dtype, a duplicate id, a dangling endpoint, and a bad schema version. Conversion tests cover the bubble round trip, the chain compaction counts, colour union, feature alignment, CSR edge order, determinism, and the genome → DS id chain on the synthetic metagenome. `tests/test_identity.py` checks dense-id lineage, internal CFA edge ids, permuted CSR edge features, labels, and colours, and the on-disk CDBG fixtures against the compactor. `tests/test_cgt_ml.py` checks the feature registry, label separation, prediction joins, CSC, sequence lookup, and loading a CGT directory that has no registry. `tests/test_cdbg_annotations.py` checks CFA transfer, compaction survival for node, internal-edge, and link annotations, namespaces, scalar and vector values, categories, explicit aggregation, provenance, the sidecar round trip, CGT alignment, missing targets, bad dtypes, and refused overwrites.

`pytest -m optional` builds a PyG `Data` object when `torch_geometric` is installed and skips otherwise.

`examples/mock_pipeline/run.py` trains the NumPy model on the bubble for a few epochs and checks that four predictions come back.

## External tools

| Tool | Role | How it is invoked |
| --- | --- | --- |
| Samovar `generate` | Contract 1 baseline (InSilicoSeq) | `samovar_generate_command` |
| MEGAHIT | Contract 2 baseline | `megahit_command` (`k` odd and ≥ 15, `--keep-tmp-files`) |
| `megahit_toolkit contig2fastg` | Intermediate assembly graph | `contig2fastg_command`; parsed by `fastg_to_cfa` |
| Flye | Contract 2 repeat graph from ONT reads | `flye_command` (`--meta` by default); the GFA is read by `gfa_to_cfa` |
| minimap2 | Contig-to-genome mapping | `minimap2_command` builds the argv; the NCBI community pipeline calls `minimap2 -x asm20` itself and skips node-level F1 when `minimap2` is not on `PATH` |
| NCBI `datasets` | Genome download for NCBI communities | `benchbuild` on a community without `--contract-only` |
| Kraken2 | `kraken2` colouring | `run_kraken2`, or an existing output through `colour_from_kraken2_file` |
| Kaiju | `kaiju` colouring | `run_kaiju`, or an existing output through `colour_from_kaiju_file` |
| AMBER | Genome-binning score | `amber.py` through `amber_command` / `score_bins`; `AMBER_SRC` or `PATH` |
| `Rscript` with `sf`, `sfnetworks`, `tidygraph` | `roxel` street graph | `Rscript -e` in `bench/data/universal/roxel.py` |
| PyTorch Geometric | Contract 7 optional backend | `run_ds(..., backend="pyg")`; install with `environment-pyg.yml` |

Bifrost is a permitted CDBG backend and is not linked in this version. `cite/tools.bib` has entries for MEGAHIT, InSilicoSeq, AMBER, PyTorch Geometric, Kraken2, Kaiju, and minimap2. Flye, Samovar, and the R packages `sf`, `sfnetworks`, and `tidygraph` have no entry there.
