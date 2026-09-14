# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

OrthoPhyl is an orthology-based phylogenomics pipeline for bacteria (and other datasets) at broad evolutionary scales. It has three layers, from lowest to highest:

1. **`OrthoPhyl.sh`** — the core pipeline. Annotates genomes (Prodigal), runs OrthoFinder, aligns/trims orthogroups, identifies single-copy orthologs, and infers species trees (FastTree / IQ-TREE / RAxML / ASTRAL). This is bash driving external bioinformatics tools.
2. **`ReLeaf.sh`** — "adds leaves" to an existing OrthoPhyl run: places new assemblies into a pre-computed phylogeny using saved HMM profiles instead of rebuilding from scratch. Reuses much of the same function library.
3. **`orthophyl_pipeline_wrapper.py`** — Python orchestrator that routes each input assembly to either ReLeaf (matches an existing database) or OrthoPhyl (novel taxon → download related genomes from NCBI → build new tree → create new database). This is the newest layer and the focus of the unit/integration test suite.

`README.md` documents the Python wrapper; `README_OrthoPhyl_ReLeaf.md` documents the two shell pipelines. Extended docs live in `readmes/`.

## Running

```bash
# Core pipeline from a directory of genome FASTAs
bash OrthoPhyl.sh -g path/to/genomes -s output_dir -t <threads>

# From pre-annotated CDS + protein dirs (comma-separated, transcripts first)
bash OrthoPhyl.sh -a path_to_transcripts,path_to_prots -s output_dir -t <threads>

# Rigor presets override most conflicting params (fast | medium | full)
bash OrthoPhyl.sh -g genomes -s out -R full -t <threads>

# Built-in test datasets (mutually exclusive with -g/-a/-s)
bash OrthoPhyl.sh -T TESTER_fasttest -t <threads>   # small orchid chloroplast set, fastest
bash OrthoPhyl.sh -T TESTER -t <threads>            # full bacterial genomes
bash OrthoPhyl.sh -T TESTER_chloroplast -t <threads>

# End-to-end smoke test (OrthoPhyl then ReLeaf on fasttest data)
bash test.sh <threads> <num_for_OrthoFinder>

# Python wrapper (routing + DB management)
python orthophyl_pipeline_wrapper.py --input assemblies.tsv --database-dir databases/ --output-dir results/ --threads 32
```

Run `bash OrthoPhyl.sh -h` for the full flag list (defined in `script_lib/arg_parse.sh`).

## Testing

The Python layer has pytest suites; the shell pipelines are covered by integration tests that actually run the tools.

```bash
conda activate orthophyl          # provides pytest + all bioinformatics tools

pytest tests/                      # unit tests (subprocess mocked, run in seconds)
pytest tests/unit/test_router.py -v   # a single module
pytest tests/ -m "not slow"        # skip network/long tests
pytest tests/ --cov=. --cov-report=term-missing

pytest -m integration -v           # runs REAL OrthoPhyl.sh/ReLeaf.sh (~5-15 min)
ORTHOPHYL_RUN_INTEGRATION=1 pytest tests/integration/   # CI form
```

Integration tests assert topology correctness via Robinson-Foulds distance (RF=0) against reference trees in `TESTER/REFERENCE_TESTER_TREES/`, and auto-skip if tools (iqtree, orthofinder, mafft, hmmbuild/hmmsearch, trimal) are missing.

Note: the target scripts (`assembly_router.py`, `orthophyl_pipeline_wrapper.py`, `create_hierarchical_database.py`) live at the repo root and in `assembly_router/`, which are not import packages. `tests/conftest.py` loads them by path via `importlib` and exposes them as fixtures (`router_module`, `wrapper_module`, `db_creator_module`). Import them that way in tests, not with `import`.

Container-aware temp handling: set `PYTEST_TMP_DIR` or `TMPDIR` when `/tmp` is read-only (common under Singularity). See `pytest.ini` and `tests/CONTAINER_TESTING.md`.

## Architecture notes

**`OrthoPhyl.sh` is a thin driver.** Nearly all logic lives in sourced libraries under `script_lib/`. The top-level `MAIN_PIPE()` in `OrthoPhyl.sh` just calls named functions in order (`SET_UP_DIR_STRUCTURE`, `PRODIGAL_PREDICT`, `ORTHO_RUN`, `TRIM`, `SCO_MIN_ALIGN`, `TREE_BUILD`, ...). To understand or change a step, find its function:
- `script_lib/functions.sh` (~1800 lines) — the core pipeline steps.
- `script_lib/functions_addem.sh` + `arg_parse_addem.sh` — ReLeaf-specific ("add 'em") functions and arg parsing.
- `script_lib/arg_parse.sh` — `USAGE` and `ARG_PARSE` for OrthoPhyl.
- `script_lib/run_setup.sh` — `tester` (test-dataset setup), `test_args`, `SET_RIGOR` (the fast/medium/full presets), `control_c` trap.
- `script_lib/bash_utils_and_aliases.sh` — aliases (`shopt -s expand_aliases` required before sourcing).

**Configuration is layered**, applied in this order (later overrides earlier):
1. `control_file.paths` — paths to external programs (ASTRAL jar, catfasta2phyml, Alignment_Assessment) and `conda activate`. Only sourced when NOT in a container (guarded by `$SINGULARITY_CONTAINER`/`$DOCKER`); in containers these are on `PATH` already.
2. `control_file.defaults` — all pipeline parameter defaults (trimming, tree methods, ANI shortlist size, etc.).
3. Command-line args (`ARG_PARSE`).
4. `-c control_file` if given — **overrides the command line**.
5. `SET_RIGOR` if `-R` given — overrides conflicting params.

Use `control_file.user` for local overrides. ReLeaf uses `control_file.ReLeaf.defaults`.

**ANI shortlisting**: if the input has more sequences than `$ANI_shortlist` (default 20, `-n`), OrthoFinder runs on a MASH-selected diverse subset, then genes are propagated to the full set via HMM profiles (`ANI_ORTHOFINDER_TO_ALL_SEQS`, saved in `OG_alignmentsToHMM/`). These saved HMMs are what ReLeaf reuses.

**Output layout** lives under the `-s` storage dir: `genomes/`, `annots_nucls/`, `annots_prots/`, and `phylo_current/` (the working dir, `$wd`), with final trees in `phylo_current/FINAL_SPECIES_TREES/`.

**Python wrapper**: `PipelineWrapper` in `orthophyl_pipeline_wrapper.py` runs in phases (initialization → routing → releaf → orthophyl → aggregation) with checkpoint/resume. It shells out to three scripts in `assembly_router/`: `assembly_router.py`, `create_hierarchical_database.py`, and `add_releaf_version.py`. These paths are hardcoded in the constructor (`orthophyl_pipeline_wrapper.py:101-105`).

**Subclade partitioning** (`--max-tree-genomes`, default 150): when a novel taxon downloads more genomes than the ceiling, the wrapper caps final tree size by MASH-partitioning the raw set into subclades named `<Taxon>_1`, `<Taxon>_2`, …, and building a tree only for the subclade(s) that actually contain a query. Key points:
- **Partition happens *before* QC.** The download step is split into `gather_filter_asms.sh --download-only` (raw FASTAs land at `<dl>/assemblies_all.TMP/*.fna`, queries staged as clustering leaves) and a per-subclade `gather_filter_asms.sh --qc-only` (CheckM2 runs here). CheckM2 therefore runs only on genomes that will enter a tree. Neither flag alone changes the script's default combined behaviour.
- **Clustering** lives in `python_scripts/subclade_partition.py`: `mash triangle -k 17 -s 5000 -E` (params MUST match everywhere), scipy `average`-linkage (UPGMA), recursive size-bounded split, deterministic size-desc/name numbering. It writes `partition_manifest.json`, a per-subclade `.msh` sketch, and a `.members.txt`. Do NOT reuse `ANI_genome_picking.py` (wrong divergence floor / merge logic).
- **Total-genome guardrail** (`--max-total-genomes`, default 5000): the partitioner builds a DENSE `NxN` MASH distance matrix, O(n²) in time and memory (~20 GB at n=50k). `_enforce_total_genome_ceiling` runs right before `_partition_genomes` in both create and batch modes and hard-stops above the ceiling rather than OOM-killing the node. The planned `--subsample` (one tree from a diverse subset) and `--megatree` (per-subclade trees merged) strategies are the supported way to handle very large taxa and will route around this check with their own bounded handling.
- **Lazy build**: subclades with no query are registered as `built=false` DB entries (via `create_hierarchical_database.py --register-only --is-subclade`) carrying their sketch + raw members + `source_genome_dir`; their tree is built on demand when a query first routes to them. `--taxon` create mode builds *all* subclades.
- **Routing** (`assembly_router.py`): all subclades share one GTDB string, so the router disambiguates by MASH sequence distance — `mash dist query.msh <subclade>.msh`, nearest member (min). An unbuilt target yields a third pipeline decision, `OrthoPhyl_subclade_build`, which the wrapper executes (QC + OrthoPhyl) before ReLeaf; a built target routes straight to ReLeaf.
- Subclade metadata fields (`is_subclade`, `parent_taxon`, `subclade_id`, `built`, `sketch_file`, `members_file`, `source_genome_dir`) load via `config.get()` defaults, so pre-existing databases stay backward-compatible.

## Repo hygiene warnings

- Historically, files were **versioned by copy** rather than git (`.v2`, `.cmd_out3`, `_v2`, `.bak`, `.save`, `.1.1` suffixes). The live scripts have since been renamed to canonical names and the dead copies removed. If you find such suffixed twins reappearing, commit over the original instead — git is the version history. `dev_plans/` still references the old names historically; don't treat those as live paths.
- Large binaries live in the tree: `*.sif` Singularity images (multi-GB), `GPLv3.pdf`, `ncbi_dataset.zip`. `.sif` is gitignored. Do not commit new large artifacts.
- Test outputs under `TESTER/` (e.g. `FULLTEST_OUT.*`) are gitignored ephemera.

## Containers

Singularity recipes: `Singularity.OP.v3.1.0.recipe` (current) and the `OP2.2.1` variants; `Dockerfile.mamba` for Docker. Both use micromamba with a base env (orthofinder, iqtree, fasttree, hmmer, mash, prodigal, trimal, raxml, ete3, pytest) and a separate `gather_genomes` env (checkm2, entrez-direct, ncbi-datasets-cli). Conda env spec: `orthophyl_env.2.2.1.yml`.
