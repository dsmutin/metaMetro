# Benchmarks

`metametro benchbuild` generates every registered benchmark graph. Other tools depend on the directory it writes. They do not keep a second copy of the generator. `scripts/phage_baseline.py` is a second, older phage generator outside `benchbuild`; the registered benchmark for the same five phages is `phage_species_5`.

Strong100, the samovar10 bundles, and the prebuilt ONT and Illumina CSR dumps are not benchmarks here. They are frozen assemblies without a generator in this repository.

## Command

```bash
metametro benchbuild --list
metametro benchbuild --all
metametro benchbuild bubble_strain_2
metametro benchbuild bubble_strain_2 --outdir /path/to/empty_dir
metametro benchbuild bubble_strain_2 --colouring as_built
```

`--all` builds every in-process graph and writes a contract for every community and external benchmark. It does not download genomes or call R. `--all --execute` also downloads and assembles those. A single benchmark is different: `metametro benchbuild <name>` downloads and assembles by default, and `--contract-only` writes only the contract. `--execute` has an effect only together with `--all`. `--all` skips a directory that already has `manifest.yaml` and `identity.sha256`, leaves it unchanged, and prints its status as `present`. A single `metametro benchbuild <name>` into a directory that is not empty stops with `ContractError` (`output directory is not empty`); so does `--all` on a non-empty directory that lacks one of those two files.

The default output directory is

```text
data/bench/{bench_name}/{assembler}/{properties}/
```

It is relative to the root of the MetaMetro checkout, not to the current directory, and that tree is gitignored. An installed package with no checkout above it (for example the conda package) has no default: pass `--outdir`, or `outdir=` / `root=` in Python, otherwise the build stops with `ContractError`. A build writes:

| Path | Role |
| --- | --- |
| `cfa/`, `cdbg/`, `cgt/` | The three graph structures |
| `ground_truth/` | Evaluation labels. Not a model input. `mock_bubble` and `mock_chain` have none. On the five hand-built bubble benches and the two phage benches it also holds the gold taxon labels, `taxon_nodes.tsv` and `taxon_edges.tsv` (see below) |
| `genomes/`, `reads/` | Read benchmarks (`bubble_reads_2`, `synthetic_reads_2`): the simulated genomes and the reads with their metadata |
| `contract/` | Community and external benchmarks: the accession pin or tool parameters |
| `work/` | Community and external benchmarks after download and assembly: downloads, assembly, classifier databases and outputs |
| `manifest.yaml` | Bench name, `colourings` (the colourings applied, in order), identity; on an assembled NCBI community also `tools`, the `--version` line of megahit, kraken2, kaiju, minimap2, and samovar (`absent` when not on `PATH`), also written to `work/commands.log`. Not part of the identity |
| `identity.sha256` | Digest of sequences, topology, colour sets, and the gold labels moved to `ground_truth/` |

The CGT carries a feature layer in `node_features.npy`, separate from colours: for each unitig, L2-normalised canonical k-mer counts (`kmer:k4:<kmer>`, 136 columns at the default `k=4`) when every unitig has a sequence, and `coverage`, the mean CFA node coverage of the unitig, when the CFA has a `coverage` column. It adds no colour namespace and does not change `node_colors` or `identity.sha256`. `manifest.yaml` records `features.kmer_k`. `load_bench_cgt(outdir)` rebuilds the same layer; `load_bench_cgt(outdir, features=False)` gives an empty feature matrix.

The `composition` colour layer is computed from this feature matrix. `composition_kmeans` runs `kmeans_onehot(node_features, n_clusters=8, seed=0)` on the CGT rows, one row per unitig. Cluster `i` is the colour value `cluster_i`, and the `composition` colour ids follow cluster order. So the `composition` columns of `node_colors` are that one-hot matrix, the same matrix ParaGVAE used to compute when a CGT had no colours. Every CFA node of a unitig gets the unitig's cluster. An edge gets the cluster when both of its ends have the same one.

Every colouring registered in `metametro.bench.colourings` that can run is applied, including a colouring added after the benchmark was introduced. Each automatic colouring is checked against the graph as coloured so far and skipped when it cannot run, so `decaying` sees the layers written before it. `--colouring` restricts that set and runs the named colourings in order. It works only for a single benchmark: `--all --colouring` is a parser error. Each one is checked against the graph as coloured so far, and one that cannot run stops the build with an error instead of being skipped. `metametro benchbuild --list-colourings` prints the registry.

| Colouring | Namespace | Auto | When it runs |
| --- | --- | --- | --- |
| `as_built` | existing (sample, route, …) | yes | always; keeps colours already on the graph except gold `taxon`, which is moved to `ground_truth/` first |
| `read_depth` | `sample` | yes | reads and integer `k` are on the build |
| `composition_kmeans` | `composition` | yes | every node has a sequence and the graph compacts to at least two unitigs |
| `kraken2` | `kraken2` | yes | Kraken2 output or `kraken2` plus a database |
| `kaiju` | `kaiju` | yes | Kaiju output or `kaiju` plus an FM-index |
| `decaying` | `decaying` | yes | the graph has edges and an allowed source layer; it leaks the first of `sample`, `kraken2`, `kaiju`, `composition` on the graph, never `taxon` |
| `read_accession` | `accession` | no | named with `--colouring`; ISS-style read ids. Simulated accessions leak truth, so this is not auto |

A half community (`*_half`) still simulates the half metagenome. Its Kraken2 and Kaiju library is the parent pin's full `db` set (the non-synonymous partners), not the half set. A community with a parent (`*_half`, `*_x10`) does not run `datasets`: build the parent first, because the derived build reuses its FASTA and its `data/raw/<parent>_assembly_data_report.jsonl` (the first ancestor without a parent). Kraken2 and Kaiju indexes need `kraken2` or `kaiju` plus Samovar on `PATH`; without them the indexes and those colourings are skipped.

`--kraken-db`, `--kraken-output`, `--kaiju-db`, and `--kaiju-output` (the `build()` and `build_all()` parameters `kraken_db`, `kraken_output`, `kaiju_db`, `kaiju_output`) pass Kraken2 and Kaiju inputs to the `kraken2` and `kaiju` colourings on every benchmark kind: in-process, community, and external. An output file is read as it is; its sequence ids must be CFA node ids, and a file in which none of them is a node id stops the build with `ContractError` instead of adding an empty layer. A database is used only when `kraken2` or `kaiju` is on `PATH`. A path given here wins over the one a community build finds under `work/`. A given path that does not exist stops the build with an error that names the flag, before anything is written. For example, `metametro benchbuild bubble_strain_2 --colouring kraken2 --kraken-output contigs.kraken` paints the `kraken2` layer on the in-process bubble.

Downstream tools do not recolour. They load the ToCUMG and name namespaces:

```python
from metametro.bench import load_bench_cdbg, namespaces_for

graph = load_bench_cdbg(bench_dir, namespaces=namespaces_for(("kraken2", "decaying")))
```

`filter_colours(graph, namespace="kraken2")` is the same selection on a CFA, CDBG, or (by `color_ids`) a CGT.

A column or colour namespace that carries an evaluation target (`truth_taxon_id`, `expected_action`, `genome_abundance`, a `truth` namespace, the gold `taxon` namespace) is rejected. Those values stay under `ground_truth/`.

### Gold taxon labels

Gold `taxon` is not a colour of a built graph, so no tool can train on it. A builder that hands `materialize` a CFA with a `taxon` namespace (the bubble benches `bubble_strain_2`, `bubble_error_1`, `bubble_nested_2`, `bubble_shared_2`, `bubble_strain_3_n50`) has those colours removed before any colouring runs and written to:

| File | Columns |
| --- | --- |
| `ground_truth/taxon_nodes.tsv` | `node_id`, `taxon` (one row per CFA node) |
| `ground_truth/taxon_edges.tsv` | `edge_id`, `source`, `target`, `taxon` (one row per CFA edge) |

A `taxon` cell is a comma-separated list of gold values, empty when there is none. Evaluation (ParaGVAE, BubbleBlower, MetaMalevich) reads them like this:

```python
from metametro.bench import load_bench_gold, load_bench_gold_cgt

gold = load_bench_gold(bench_dir)          # GoldLabels(namespace="taxon", values, nodes, edges)
gold.nodes["S"]                            # ("taxon_1", "taxon_2"), keyed by CFA node id
labels = load_bench_gold_cgt(bench_dir)    # a CGT whose only colours are the gold labels
labels.node_colors                         # (num_nodes, len(gold.values)) uint8, column j = gold.values[j]
labels.edge_colors                         # CSR edge order of load_bench_cgt(bench_dir)
```

`load_bench_gold_cgt` has the node and CSR edge order of `load_bench_cgt(bench_dir)`; a unitig's gold set is the union over its CFA members. It replaces `load_bench_cgt(bench_dir, namespaces=["taxon"])`, which no longer finds a `taxon` namespace. Use it for scoring only.

On the bubble benches `decaying` therefore leaks `composition` (or `kraken2` / `kaiju` when that layer was painted), not gold taxon. On `mock_bubble`, `mock_chain`, `bubble_reads_2`, and `synthetic_reads_2` it leaks `sample`.

## Layout

Generation code:

```text
src/metametro/bench/data/universal/     shared builders
src/metametro/bench/data/pins/{bench_name}/  accession pin of one NCBI community
```

Scoring code:

```text
src/metametro/bench/scoring/assembly/universal/
src/metametro/bench/scoring/assembly_binning/universal/
src/metametro/bench/scoring/profiling/universal/
src/metametro/bench/scoring/{type}/{bench_name}.py   scorer for one benchmark, e.g. assembly/bubble_strain_2.py
```

`assembly_binning` is the contig-binning score used by ParaGVAE (majority-label contig F1) plus AMBER sequence F1. AMBER is `amber.py` from [CAMI-challenge/AMBER](https://github.com/CAMI-challenge/AMBER) tag `v2.0.17-beta` (`089ea20e83811e090ab69814e033adc9fc963892`). Set `AMBER_SRC` to that checkout, or put `amber.py` on `PATH`. The required tests do not call it. `profiling` is the abundance score (L1, Bray–Curtis, presence F1). `assembly` is a length summary such as N50. A scorer reads `ground_truth/` and a prediction table. It does not read a target out of the feature matrix.

## In-process benchmarks

| Name | What it is |
| --- | --- |
| `bubble_strain_2` | Two taxa, one bubble. The BubbleBlower toy. Truth: retain |
| `bubble_error_1` | One taxon, one error bubble. Truth: pop |
| `bubble_nested_2` | Two taxa, outer and inner bubbles |
| `bubble_shared_2` | Two taxa, one shared node, no bubble |
| `bubble_strain_3_n50` | Three strains, 50 bubbles, seed 42 |
| `bubble_reads_2` | Two strains, read length 12, de Bruijn k=5 |
| `mock_bubble` | MetaMetro four-node fixture |
| `mock_chain` | MetaMetro six-node compaction fixture |
| `synthetic_reads_2` | Two synthetic genomes, in-process reads, k=8 |

On `bubble_reads_2` and `synthetic_reads_2` the CFA node `coverage` is how often the node k-mer occurs in all reads of the bench.

Legacy names such as `strain_bubble` resolve to the canonical name.

Pinned NCBI communities and external contracts are in the same registry. `metametro benchbuild NAME` downloads and assembles when `datasets`, Samovar, and MEGAHIT are on `PATH`. `metametro benchbuild NAME --contract-only` writes the accession pin and the lognormal abundance table and does not invent a graph. A second contract build has the same `identity.sha256`. The download step writes its manifest, checksums, acquisition report, and dataset audit to `work/download/` inside the bench directory, never to tracked files under `data/` or `docs/`.

| Name | Was | Design |
| --- | --- | --- |
| `4domain_family_100` | high100 | 100 families, 4 domains, family rank, 100000 reads |
| `4domain_family_100_half` | half100half | every other family, 50 families |
| `4domain_family_100_x10` | high100_enriched | same families, 1000000 reads |
| `4domain_family_100_half_x10` | half100half_enriched | 50 families, 1000000 reads |
| `3domain_genus_75` | low75 | 75 genera, 3 domains, genus rank |
| `3domain_genus_75_half` | low75half | 38 genera |
| `3domain_genus_75_x10` | low75_enriched | 75 genera, 1000000 reads |
| `3domain_genus_75_half_x10` | low75half_enriched | 38 genera, 1000000 reads |
| `bacteria_species_20_heldout` | heldout_genera | 4 genera, 20 species pairs, species rank |
| `bacteria_strain_10` | half_strains | 10 of those pairs |
| `phage_species_5` | phage_baseline | T1, T3, T4, T5, T7, 400 reads, k=21 |
| `phage_species_5_x10` | phage_x10, phage_10 | same five phages, 4000 reads |
| `roxel` | roxel | sfnetworks street graph |
| `spb_ground_transit` | spb transit | ORGP GTFS, 60 m stop merge; pass `--gtfs` |

There is no ten-species phage benchmark. `phage_10` was the ten-times read budget. Strong100, samovar10, and the prebuilt ONT and Illumina bundles are not in this registry.

The phage build runs MEGAHIT with one thread. With several threads MEGAHIT writes the same contigs in a different order on each run, so the graph identity would change between two builds of the same reads.

### Phage genomes

`phage_species_5` and `phage_species_5_x10` read five RefSeq genomes from `data/raw/genomes/<name>.fna` in the checkout. That directory is gitignored. The pin is `src/metametro/bench/data/pins/phage_species_5/genomes.tsv`. It gives one versioned accession per genome, the NCBI taxid of the record, its species, the length, the SHA-256 of the NCBI FASTA, and the SHA-256 of the file the build reads.

| Name | Accession | Taxid | Organism | Species taxid | Species | Length (bp) |
| --- | --- | --- | --- | --- | --- | --- |
| T1 | NC_005833.1 | 2492962 | Escherichia phage T1 | 1921008 | Tunavirus T1 | 48836 |
| T3 | NC_003298.1 | 10759 | Enterobacteria phage T3 | 2732706 | Teetrevirus T3 | 38208 |
| T4 | NC_000866.4 | 2681598 | Escherichia phage T4 | 10665 | Tequatrovirus T4 | 168903 |
| T5 | NC_005859.1 | 2695836 | Escherichia phage T5 | 10726 | Tequintavirus T5 | 121750 |
| T7 | NC_001604.1 | 10760 | Escherichia phage T7 | 1985738 | Teseptimavirus T7 | 39937 |

Taxids and names are from NCBI nuccore `esummary` and taxonomy `efetch` (2026-11-20). The record taxid has no rank in NCBI taxonomy. The species taxid is its ICTV species.

```bash
metametro benchfetch phage_10            # download, check SHA-256, write data/raw/genomes/T*.fna
metametro benchfetch phage_10 --list     # print the pin, no download
metametro benchfetch phage_10 --genome-dir DIR
```

`benchfetch` calls NCBI E-utilities `efetch` (`db=nuccore`, `rettype=fasta`) over HTTPS, one request per accession. It writes each record as `>T1` and one upper-case sequence line. Bases other than `ACGTN` become `N`; only T4 has any (IUPAC `M`, `S`, `R` at 1-based 56143, 70346, 135175). A file is written only after its length and SHA-256 match the pin. A file that already matches is kept and reported as `present`, so a second run downloads nothing. A file with another hash stops `benchfetch` until you pass `--force`. In Python the same step is `metametro.bench.genomes.fetch_genomes(name, genome_dir)`.

`metametro benchbuild phage_10` checks the five files before it writes anything. A missing or changed genome raises `ContractError`, for example:

```text
phage_species_5_x10 is missing pinned genome FASTA under <checkout>/data/raw/genomes: T1 (NC_005833.1), T3 (NC_003298.1), T4 (NC_000866.4), T5 (NC_005859.1), T7 (NC_001604.1); fetch them with: metametro benchfetch phage_species_5_x10
```

A changed file is listed as `differs from the pinned SHA-256`, and the command gets `--force`. Then the build checks for `samovar`, `megahit`, and `megahit_toolkit` on `PATH`. `phage_species_5_x10` built from these genomes has identity `78121c1b05a5814bfc17ac3b7650ed5223a0f584cfada24e9224b7ba0364564c` (983 CFA nodes, 1270 edges; SAMOVAR 0.12.150, MEGAHIT 1.2.9; colourings `as_built`, `composition_kmeans`, `kraken2`, `decaying`).

### Phage gold labels

A phage build writes gold `taxon` labels to `ground_truth/taxon_nodes.tsv` and `taxon_edges.tsv` (`metametro.bench.phage_gold`). A CFA node gets the pinned genome that contains a strict majority of its k-mers (k = 21, both strands), the same rule as the `node_labels` of `scripts/phage_baseline.py`; any other node is empty. The value is that genome's species taxid from `genomes.tsv`. An edge gets the label its two ends share. `ground_truth/taxon_values.tsv` names each value: `taxon`, `rank` (`species`), `name`, `genome`, `accession`, and `record_taxid`. The labels are for scoring only and are not a colour of any graph. On `phage_species_5_x10`, 903 of 983 nodes and 440 of 1270 edges have a label (T4 389, T5 177, T7 163, T3 156, T1 18 nodes). The gold labels are part of the identity: they changed it from `3d1ce381…` to `af194d92…` without changing the CFA, CDBG, or CGT files, and the pinned `kraken2` layer below changed it to the value above.

### Phage Kraken2 layer

`phage_species_5_x10` paints the `kraken2` colour layer from a Kraken2 output pinned in the repository, so the build needs neither `kraken2` nor its index. The pin is `src/metametro/bench/data/pins/phage_species_5_x10/kraken2/`: `nodes.kraken` (Kraken2's per-sequence output for the CFA node FASTA; sequence ids are node ids) and `kraken2.yaml`.

| Field | Value |
| --- | --- |
| Kraken2 | 2.17.2 (bioconda `kraken2 2.17.2 pl5321h3be2455_0`) |
| Index | `k2_viral_20261022` (Ben Langmead genome-idx, RefSeq viral, 2026-10-22) |
| Index archive SHA-256 | `1532b484559347bed98c686cb6b0671910f971d358c3873ba184d7ae6778cfa4` (589941425 bytes) |
| Command | `kraken2 --db k2_viral_20261022 --threads 1 --output nodes.kraken nodes.fna` |
| Input `nodes.fna` SHA-256 | `3b6d88efb06652ed6b167e7d587715cff9fb9449593d64a842d53f6c9a32b6e2` (983 sequences) |
| Output `nodes.kraken` SHA-256 | `4596d7687ef4fce1b49ed614521c4c364a28071ae92868ff374d035c2208b1e0` |
| Classified | 889 of 983 nodes (90.4%) |

The build writes the node FASTA to `work/classifier/nodes.fna` and uses the pin only when that file has the pinned SHA-256; otherwise it stops with `ContractError` and asks for `--kraken-db` or `--kraken-output`. A given `--kraken-db` or `--kraken-output` wins over the pin. The layer keeps every taxid with k-mer hits on a node (`parse_kraken2_output`), so 894 of 983 nodes and 548 of 1270 edges carry a `kraken2` colour, 44 taxids in all. Its values are NCBI taxids at the rank Kraken2 called (often the genus or the record taxid), not the gold species taxids; `metametro.bench.kraken_pin.kraken_pin(name)` returns the pinned path and its metadata. `phage_species_5` has no pin.

