# Contributing to metametro

All development follows this guide. Project rules require it (`follow-contributing`).

English is required for every public function, class, module, CLI flag, data file, and user-facing document. Other languages or missing documentation are not allowed.

## Architecture

```
CLI (metametro.cli)
  metametro benchbuild NAME  → metametro.bench.cli → metametro.bench.build
      → data/bench/{name}/{assembler}/{properties}/  (cfa/, cdbg/, cgt/, ground_truth/, manifest.yaml, identity.sha256)
  metametro (no subcommand)  → placeholder baseline (metametro.baseline.run_pipeline)
      → JSON result {status, ok, input_path}

src/metametro/
  formats/        CFA, CDBG, CGT: models, loaders, dumpers, validators, schema.yaml
  converters/     CFA → CDBG → CGT, CDBG → CFA, CSR → CGT, CGT → PyG
  contracts/      stage code: simulation and assembly, colouring, classifiers, DS, external commands, CONTRACTS
  bench/          benchmark registry, build, colourings, leakage guard, loaders, scoring
  viz/            CFA colouring and graph-state plots
  edits.py, identity.py, fixtures.py, states.py, digraph_cfa.py, transit_*.py

tests/          mandatory vs optional pytest
examples/toy/   end-to-end run of the current (baseline) tool
scripts/        external baselines that pytest does not run (Samovar, MEGAHIT)
docs/           formats, contracts, implementation, audits, run records
data/           manifests, checksums, and gitignored raw and work trees
cite/           BibTeX for integrated third-party tools
agents/         portable rules and skills (any IDE)
```

The GitHub wiki is a separate repository (`MetaMetro.wiki.git`). Clone it beside this one, not inside it:

```bash
git clone https://github.com/dsmutin/MetaMetro.wiki.git ../MetaMetro.wiki
```

Replace baseline bodies with real implementations. Keep the documented return keys until you change `docs/contracts.md`, the wiki page that mirrors it, and the tests together.

## Testing architecture

| Kind | Marker | Command | When |
|------|--------|---------|------|
| Required | `mandatory` | `pytest -m mandatory` | every commit; GitHub Action `required-tests` |
| Optional | `optional` | `pytest` (all) | push to main, release, or workflow_dispatch; Action `full-tests` |
| Examples | — | `python examples/toy/run.py` | full CI; after features that touch the CLI |
| Extra examples | — | extra `examples/*/run.py` or `pipeline.sh` | full CI when those files exist |
| phage_10 solve | — | job `phage10-solve`: `benchfetch phage_10`, a default build checked against `canonical_identity("phage_10")`, then metaMalevich `solve` at a pinned commit; `scripts/check_phage10_solve.py` requires the phage_10 row to have `build_status` `built`, `solved` true, `kraken2` applied from the MetaMetro pin with no colouring skipped, and metrics, and checks the identity and gold digest of the build solve made, ignoring the exit code and `ok` | push to main, release, workflow_dispatch, and every pull request to main (it is not a required check unless branch protection makes it one) |

Do not mark a contract test `optional`. Optional tests are slow, extra, or nice-to-have.

After **any new feature**, run the **mandatory** suite (and toy if the CLI changed) before you stop.

## Feature checklist (`todo.md`)

Track work in `todo.md` (checkboxes). One line per feature or fix. Check it off only when mandatory tests pass. This is a **feature list**, not a `/do` analysis graph.

## Versioning

Edit **only** `VERSION`. Everything else reads it.

Starting value: `0.0.1`.

| Change | Bump |
|--------|------|
| New feature | **minor** (`0.0.1` → `0.1.0`) |
| Fix or update of an existing feature | **patch** (`0.1.0` → `0.1.1`) |
| Release | **major** (`0.1.1` → `1.0.0`) |

A major bump is made only by the owner's decision, and only for an incompatible change of a contract (`docs/contracts.md`, a file format, or a public function's behaviour that downstream tools rely on). An agent or contributor never makes one on their own; a breaking change is proposed first and waits for that decision.

Every version bump is its own commit, and that commit gets an annotated tag `vX.Y.Z` with the message `MetaMetro X.Y.Z`, matching `VERSION`. Push the tag together with the commit. Never move, rewrite, or delete a pushed tag; a mistake is fixed by the next version. Downstream packages pin a tag, not a branch.

## Install (conda only)

```bash
conda env create -f environment.yml
conda activate metametro
python -m pip install --no-deps --no-build-isolation -e .
```

The `pip` step installs the package into the conda environment (no dependencies are pulled from PyPI) and gives the `metametro` command, as in CI.

`environment.yml` is the required environment and pins exact versions. `environment-pyg.yml` adds the optional PyTorch Geometric backend (`metametro-pyg`); keep it out of the required environment so `required-tests` stays light. That file cannot reuse the `setuptools` pin because conda-forge `pytorch` requires `setuptools <82`.

## GitHub

Never `git push` unless the human explicitly asks. CI runs on GitHub after they push. A push straight to `main` without a pull request (see [Benchmarks](#benchmarks)) also needs that explicit request for the task.

## Benchmarks

New benchmarks are added in this repository, under `src/metametro/bench/`, and reviewed as a pull request. Do not start a second benchmark tree in a downstream tool.

When the owner allows it for a task, a benchmark change may go straight to `main` without a pull request. Then the mandatory tests (`pytest -m mandatory`) must pass before the push, and the second reviewer reviews the pushed commits before the task is closed.

Do not hard-code a machine path (`/mnt`, `/nfs`, `/home`, or a drive letter). A benchmark input is the `benchbuild` output directory, a CLI argument, or an environment variable. If it is missing, stop.

Do not mock a benchmark graph, a taxonomy label, or a metric. A tiny graph written inside a test may stay in that test. A benchmark comes from `metametro benchbuild`.

Do not copy an evaluation target into graph features, colours, or any file a model reads as input. Simulated taxon ids, expected debubbler actions, and gold abundances stay in `ground_truth/` and are used only by the scorer. `benchbuild` rejects a graph that carries those columns.

See [docs/benchmarks.md](docs/benchmarks.md).

Genome-binning scores call `amber.py` from <https://github.com/CAMI-challenge/AMBER> at tag `v2.0.17-beta`. Clone that tag and set `AMBER_SRC`, or put `amber.py` on `PATH`. Do not replace it with a local reimplementation.

## Colourings

New colouring methods belong in this repository, under `src/metametro/bench/colourings.py` and `src/metametro/contracts/`, and are reviewed as a pull request. Downstream tools (BubbleBlower, metaMalevich, ParaGVAE) do not add a second colouring implementation.

Do not mock a colouring. A tiny graph that exists only inside a test may stay in that test. A classifier colouring is applied from a real Kraken2 or Kaiju output, or from a fixture that is that tool's file format, not from invented taxon masks.

A new colouring must declare when it is available, which colour namespace it writes, and whether it auto-applies. `benchbuild` applies every auto colouring that can run. Check that the new namespace does not collide with an existing one and that `filter_colours` / `load_bench_cdbg(..., namespaces=...)` can select it. Simulated taxon ids stay out of colours.

Downstream tools select layers from the ToCUMG:

```python
from metametro.bench import load_bench_cdbg, namespaces_for

graph = load_bench_cdbg(bench_dir, namespaces=namespaces_for(("kraken2", "decaying")))
```

## Citations

Add a `.bib` entry in `cite/` only for tools this package actually integrates. Do not invent papers.
