# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

OrthoPhyl is an orthology-based phylogenomics pipeline for bacteria (and other datasets) at broad evolutionary scales. It has three layers, from lowest to highest:

1. **`OrthoPhyl.sh`** — the core pipeline. Annotates genomes (Prodigal), runs OrthoFinder, aligns/trims orthogroups, identifies single-copy orthologs, and infers species trees (FastTree / IQ-TREE / RAxML / ASTRAL). This is bash driving external bioinformatics tools.
2. **`ReLeaf.sh`** — "adds leaves" to an existing OrthoPhyl run: places new assemblies into a pre-computed phylogeny using saved HMM profiles instead of rebuilding from scratch. Reuses much of the same function library.
3. **`orthophyl_pipeline_wrapper.py`** — Python orchestrator with three modes: batch (`--input`, routes each assembly to either ReLeaf if it matches an existing database, or OrthoPhyl to download related genomes from NCBI and build a new one), taxon (`--taxon`, same NCBI-download-and-build path for a single named clade), and local-genome (`--genome-dir`, builds a database from genomes already on disk under a user-supplied clade name). This is the newest layer and the focus of the unit/integration test suite.

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

# Python wrapper: build a database from genomes already on disk (QC runs by default)
python orthophyl_pipeline_wrapper.py --genome-dir path/to/fastas --clade-name MyIsolates \
  --database-dir databases/ --output-dir results/ --threads 32
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

**Large-taxon handling** (`--max-tree-genomes`, default 2000): when a novel taxon downloads more raw genomes than the ceiling, the wrapper caps final tree size. The **default** is **diverse subsampling** — `python_scripts/subsample_genomes.py` sketches all genomes once (`mash sketch`, linear) then greedily picks the most-diverse `--subsample-size` (default 500) via farthest-point (max-min) selection, one `mash dist` call per pick (O(n·N) time, O(n) memory — no dense matrix, scales to huge taxa). Query/must-keep genomes seed the pick so they're always retained. `_subsample_genomes` stages the selected FASTAs into `<orthophyl_dir>/subsample/<taxon>/selected/`, which then flows through the normal single-tree QC+build path. The `greedy_maxmin(names, dist_fn, target, seed_names)` core is pure (inject a distance oracle to unit-test without mash), mirroring `partition_matrix`.

The **partition/megatree** path is opt-in via `--megatree`. Instead of subsampling an oversized taxon down to one tree, it MASH-partitions the raw set into subclades named `<Taxon>_1`, `<Taxon>_2`, …, builds a full tree for **every** subclade, builds a small BACKBONE tree from a few diverse reps per subclade, and GRAFTS each subclade's full tree onto its reps in the backbone — one merged tree containing every genome. `_run_megatree` orchestrates it; the merged tree + conflict report land in `03_results/trees/orthophyl/<taxon>_megatree.nwk` / `_megatree_conflicts.json`. Key points:
- **Engaged only when raw count > `--max-tree-genomes`** and `--megatree` is set; otherwise it collapses to the normal single-tree build. Not mutually exclusive with subsampling (subsample is the default when `--megatree` is absent).
- **Partition happens *before* QC.** The download step is split into `gather_filter_asms.sh --download-only` (raw FASTAs land at `<dl>/assemblies_all.TMP/*.fna`, queries staged as clustering leaves) and a per-subclade `gather_filter_asms.sh --qc-only` (CheckM2 runs here). CheckM2 therefore runs only on genomes that will enter a tree. Neither flag alone changes the script's default combined behaviour.
- **Clustering** lives in `python_scripts/subclade_partition.py`: `mash triangle -k 17 -s 5000 -E` (params MUST match everywhere), scipy `average`-linkage (UPGMA), recursive size-bounded split (megatree passes `--max-size self.subclade_size`, default 150), deterministic size-desc/name numbering. It writes `partition_manifest.json`, a per-subclade `.msh` sketch, and a `.members.txt`. Do NOT reuse `ANI_genome_picking.py` (wrong divergence floor / merge logic).
- **Grafting** lives in `python_scripts/megatree_graft.py` (pure ete3, unit-tested by injecting newick strings). `graft(backbone, subclade_trees, rep_map, min_support)` does MRCA-replace with a monophyly check (non-monophyletic reps → best-effort anchor-graft, foreign leaves preserved, flagged). `induced_conflicts(...)` flags high-support (≥ `--conflict-min-support`, default 90) bipartition disagreements between a subclade tree and the backbone — flagging only, NOT resolution (reconciliation is future work). Backbone reps = `min(subclade_size, --backbone-reps)` (default 5) per subclade via that subclade's MASH greedy max-min. The taxon DB is created from the **backbone** run so ReLeaf has a coherent HMM set. `--taxon` create mode takes the same path (no query → all subclades built).
- **Total-genome guardrail** (`--max-total-genomes`, default 25000): the *partitioner* builds a DENSE `NxN` MASH distance matrix, O(n²) in both time and memory (measured ~8 hours wall-clock and ~45 GB RAM at n=75k, 12 threads). `_enforce_total_genome_ceiling` hard-stops the megatree path above the ceiling rather than running for hours or OOM-killing the node. The **default subsample path never builds this matrix and is not subject to this check** — it handles arbitrarily large taxa. A query-neighborhood-mode ToDo exists for taxa too large even for this path, sidestepping the O(n²) matrix entirely.
- **Routing** (`assembly_router.py`) has three decisions: **ReLeaf** (matches an existing DB), **OrthoPhyl** (novel taxon), or **OrthoPhyl_subclade_build** (query matched a lazily-registered subclade that has no tree yet — build it, then ReLeaf). Every megatree DB shares its parent's taxonomy string with its siblings (the backbone plus every `<Taxon>_N` subclade), so a query can tie at the same specificity across several DBs; `route_assembly` disambiguates ties by `--placement` and MASH:
  - `--placement subclade` (default): among tied subclade DBs, sketch the query (`mash sketch -k 17 -s 5000`, same params as `subclade_partition.py` — `MASH_K`/`MASH_S` module constants in `assembly_router.py` MUST match) and pick the nearest via `mash dist` (`_route_subclade_by_mash` / pure `pick_nearest_subclade`, ties broken by ascending `clade_name`). Falls back to the backbone DB if no subclade sketch is usable.
  - `--placement backbone`: routes straight to the tied parent's `is_backbone` DB (the sparse megatree overview tree) with no MASH call.
  - A DB with `is_subclade=True` and `built=False` (registered but never built — see below) yields `OrthoPhyl_subclade_build` instead of `ReLeaf`, carrying the subclade's own `clade_taxonomy` as `subclade_taxonomy` (NOT the query's) so a later rebuild doesn't overwrite it.
- **Lazy build-on-demand** (`--megatree-lazy`, opt-in, off by default): at partition time, a subclade with no query assigned to it is only *registered* (`create_hierarchical_database.py --register-only`, requires `--is-subclade`) rather than built — `built=false`, a placeholder `phylogeny.nwk`, no `orthophyl_run` symlink, `genome_list.txt` sourced from the subclade's `members_file`. Its sketch/members/`source_genome_dir` are recorded so a later query that MASH-matches it can trigger an on-demand build (`_register_lazy_subclade` → `_phase_subclade_build`/`_process_subclade_build` → `_build_subclade(force=True)`, which promotes the placeholder to `built=true`). Registration and build use *separate* checkpoint keys (`register_<name>` vs `database_<name>`) — sharing one key would let `--resume` mistake a placeholder for a finished build. Without `--megatree-lazy`, every subclade is still built up front (the pre-existing default).
- Subclade/backbone metadata fields (`is_subclade`, `is_backbone`, `parent_taxon`, `subclade_id`, `built`, `sketch_file`, `members_file`, `source_genome_dir`) load via `config.get()` defaults, so pre-existing databases stay backward-compatible. The megatree backbone DB is written with `is_backbone=true` (via `--is-backbone`) so the router can tell it apart from a dense subclade sharing the same parent-level taxonomy.

**Local-genome ingest** (`--genome-dir DIR --clade-name NAME`): builds a routable database from genomes the user already has, instead of downloading from NCBI. `_run_local_genomes_mode` in `orthophyl_pipeline_wrapper.py`:
- **QC runs by default**, skippable with `--skip-qc` (the request explicitly wanted QC on, not off, by default). Skipping sets `qc_applied=false` in the database config and is only ever honored on this path — `_qc_subclade`/the NCBI download path always QC.
- **Taxonomy resolution** (`_resolve_local_taxonomy`) tries, in order: (1) `--clade-taxonomy` verbatim (an escape hatch for a full GTDB string), (2) resolving `--clade-name` against the local NCBI taxdump to render a full lineage, (3) falling back to a bare `<rank>__<name>` (default rank `g`, `--clade-rank`) with a logged warning and a paste-ready `--clade-taxonomy` hint. A name-only fallback taxonomy is **not routable** — `GTDBTaxonomy.is_within_clade` (`assembly_router.py`) compares every rank down to the query's, so upstream `None` ranks never match a fully-specified query. This is why resolution against the taxdump matters, not just cosmetic padding.
- Every input FASTA is normalized to `<stem>.fna` (gunzip `.gz`, symlink otherwise) into a wrapper-owned staging dir — originals in `--genome-dir` are never touched or renamed, because OrthoPhyl.sh rewrites contig names in its `-g` input in place.
- `taxonomy_source` (`"ncbi"` | `"user_supplied"`) and `qc_applied` (bool) are written to every `database_config.json`, both read with `.get()` defaults (`"ncbi"`/`True`) so pre-existing databases stay valid. `assembly_router.py`'s startup log appends `[user-supplied taxonomy: not NCBI-assigned]` for `user_supplied` DBs — provenance is surfaced, not used to gate routing.
- `--genome-dir` + `--megatree` together is rejected by `main()` — the local path doesn't (yet) support per-subclade taxonomy/backbone DB creation.

## Repo hygiene warnings

- Historically, files were **versioned by copy** rather than git (`.v2`, `.cmd_out3`, `_v2`, `.bak`, `.save`, `.1.1` suffixes). The live scripts have since been renamed to canonical names and the dead copies removed. If you find such suffixed twins reappearing, commit over the original instead — git is the version history. `dev_plans/` still references the old names historically; don't treat those as live paths.
- Large binaries live in the tree: `*.sif` Singularity images (multi-GB), `GPLv3.pdf`, `ncbi_dataset.zip`. `.sif` is gitignored. Do not commit new large artifacts.
- Test outputs under `TESTER/` (e.g. `FULLTEST_OUT.*`) are gitignored ephemera.

## Containers

Singularity recipes: `Singularity.OP.v3.1.0.recipe` (current) and the `OP2.2.1` variants; `Dockerfile.mamba` for Docker. Both use micromamba with a base env (orthofinder, iqtree, fasttree, hmmer, mash, prodigal, trimal, raxml, ete3, pytest) and a separate `gather_genomes` env (checkm2, entrez-direct, ncbi-datasets-cli). Conda env spec: `orthophyl_env.2.2.1.yml`.
