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

Note: the target scripts (`assembly_router.py`, `orthophyl_pipeline_wrapper.py`, `OP_database_tool.py`) live at the repo root and in `assembly_router/`, which are not import packages. `tests/conftest.py` loads them by path via `importlib` and exposes them as fixtures (`router_module`, `wrapper_module`, `db_creator_module`). Import them that way in tests, not with `import`.

Container-aware temp handling: set `PYTEST_TMP_DIR` or `TMPDIR` when `/tmp` is read-only (common under Apptainer/Singularity). See `pytest.ini` and `readmes/README_APPTAINER.md`.

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

**Python wrapper**: `PipelineWrapper` in `orthophyl_pipeline_wrapper.py` runs in phases (initialization → routing → releaf → orthophyl → aggregation) with checkpoint/resume. It shells out to three scripts in `assembly_router/`: `assembly_router.py`, `OP_database_tool.py`, and `add_releaf_version.py`. These paths are hardcoded in the constructor (`orthophyl_pipeline_wrapper.py:101-105`).

**Large-taxon handling** (`--max-tree-genomes`, default 2000): when a novel taxon downloads more raw genomes than the ceiling, the wrapper caps final tree size. Default is diverse subsampling (`python_scripts/subsample_genomes.py`, pure `greedy_maxmin` core, mirrors `partition_matrix`). `--megatree` (opt-in) instead partitions into subclades and grafts a full tree per subclade onto a backbone so every genome survives. Full mechanics, flags, and defaults: `readmes/README_ADVANCED_FEATURES.md` §10 (Large-Taxon Handling) and §11 (External HMM Assignment). Implementation notes not covered there:
- Partitioning (`python_scripts/subclade_partition.py`) uses `mash triangle -k 17 -s 5000 -E` + scipy UPGMA; these params MUST match `assembly_router.py`'s `MASH_K`/`MASH_S` constants and `script_lib/functions.sh`'s sketch params exactly, or sketch comparisons are meaningless. Do NOT reuse `ANI_genome_picking.py` here (wrong divergence floor / merge logic).
- Grafting (`python_scripts/megatree_graft.py`) is pure ete3, unit-tested via injected newick strings; `induced_conflicts` only flags high-support disagreements, never resolves them.
- `--megatree-lazy` registration and build use *separate* checkpoint keys (`register_<name>` vs `database_<name>`) — sharing one would let `--resume` mistake a placeholder for a finished build. Registration also stages backbone reps immediately (only the subclade's own tree is deferred) — otherwise a pure `--taxon`/`--genome-dir` create-mode run (no query anywhere) would starve the backbone of input and crash the grafter; `_run_megatree` fails fast with a `RuntimeError` if that still happens.
- `--megatree-hmm-reuse` requires building the backbone before any subclade (reversing the default order) so `OG_alignmentsToHMM/hmms_final/` exists to pass as `--hmm-assign-dir`; falls back to independent builds with a warning if the backbone never crosses `--ani-shortlist`. Leftover (unmatched-HMM) genes get their own OrthoFinder pass by default, filtered through `OG_sco_filter.py` (threshold 1 — a fresh clustering has no paralog guarantee) and renamed `OG0_LFT_*` to avoid clobbering the external HMM set's OG IDs.
- Routing (`assembly_router.py`): a megatree backbone and its subclades share one taxonomy string, so ties are broken by `--placement` (`subclade` default: MASH-nearest via `pick_nearest_subclade`; `backbone`: no MASH call). A registered-but-unbuilt lazy subclade routes to `OrthoPhyl_subclade_build` instead of `ReLeaf`.

**Tree-building speed controls** (`--tree-method`, `--iqtree-model`, `--iqtree-mset`, `--iqtree-fast`, `--iqtree-options`): apply to OrthoPhyl builds only — every one of them (`--taxon`, `--genome-dir`, megatree subclade/backbone builds, novel-taxon builds in batch mode) funnels through the single `_run_orthophyl` call site, which is where these are wired in; ReLeaf placement (`_run_releaf`) is untouched. `--tree-method` maps to OrthoPhyl.sh's `-p` (`iqtree`/`fasttree`/`"fasttree iqtree"`). Any `--iqtree-*` flag makes `_tree_control_file` write a one-line control file overriding `IQtree_CDS_speciestree_options` (CDS only — the wrapper always passes `-o CDS`) and passes it as `-c`, exploiting the `control_file.defaults`-before-`-c` ordering noted above (`OrthoPhyl.sh:79-116`) — so the file can append to the existing default via a literal `$IQtree_CDS_speciestree_options` reference rather than needing to know its contents. `--iqtree-options` replaces that base outright instead of appending. `--iqtree-model` (skip ModelFinder with a fixed model, the single biggest IQ-TREE speedup) and `--iqtree-mset` (constrain which models ModelFinder tests) are mutually exclusive. `--tree-method fasttree`/`both` logs a warning: ReLeaf can't yet extend a FastTree species tree (`functions_addem.sh`'s `ADD_2_fasttree` is unfinished), and FastTree branch support is on a 0–1 scale so `--megatree`'s `--conflict-min-support` (90-scale) never flags a FastTree subclade/backbone conflict. `--iqtree-fast` is rejected (both in `__init__` and `main()`, against the same `_DEFAULT_IQTREE_CDS_OPTIONS` constant) unless the base string it appends `--fast` to is free of `-B`/`-bb` — verified against real IQ-TREE 2.2.0.3, which hard-errors ("Ultrafast bootstrap does not work with -fast...") on `--fast` + UltraFast Bootstrap, and the default `IQtree_CDS_speciestree_options` always includes `-B 1000`, so `--iqtree-fast` alone would otherwise crash every run; the fix is pairing it with `--iqtree-options` that drops `-B`/`-bb` (e.g. `--iqtree-options "--symtest --lmap 2000 -t PARS --ninit 2"`).

**Local-genome ingest** (`--genome-dir DIR --clade-name NAME`): builds a routable database from genomes the user already has, instead of downloading from NCBI. `_run_local_genomes_mode` in `orthophyl_pipeline_wrapper.py`:
- **QC runs by default**, skippable with `--skip-qc` (the request explicitly wanted QC on, not off, by default). Skipping sets `qc_applied=false` in the database config and is only ever honored on this path — `_qc_subclade`/the NCBI download path always QC.
- **Taxonomy resolution** (`_resolve_local_taxonomy`) tries, in order: (1) `--clade-taxonomy` verbatim (an escape hatch for a full GTDB string), (2) resolving `--clade-name` against the local NCBI taxdump to render a full lineage, (3) falling back to a bare `<rank>__<name>` (default rank `g`, `--clade-rank`) with a logged warning and a paste-ready `--clade-taxonomy` hint. A name-only fallback taxonomy is **not routable** — `GTDBTaxonomy.is_within_clade` (`assembly_router.py`) compares every rank down to the query's, so upstream `None` ranks never match a fully-specified query. This is why resolution against the taxdump matters, not just cosmetic padding.
- Every input FASTA is normalized to `<stem>.fna` (gunzip `.gz`, symlink otherwise) into a wrapper-owned staging dir — originals in `--genome-dir` are never touched or renamed, because OrthoPhyl.sh rewrites contig names in its `-g` input in place.
- `taxonomy_source` (`"ncbi"` | `"user_supplied"`) and `qc_applied` (bool) are written to every `database_config.json`, both read with `.get()` defaults (`"ncbi"`/`True`) so pre-existing databases stay valid. `assembly_router.py`'s startup log appends `[user-supplied taxonomy: not NCBI-assigned]` for `user_supplied` DBs — provenance is surfaced, not used to gate routing.
- `--genome-dir` + `--megatree` **is supported**: when the staged genome count exceeds `--max-tree-genomes`, `_run_local_genomes_mode` forks into `_run_megatree` exactly like `--taxon` create mode's fork (no query assemblies either way, so every subclade is built/registered immediately). `taxonomy_source='user_supplied'`/`qc_applied=True` are threaded through every DB-creation call `_run_megatree` makes (backbone, every subclade's `_build_subclade`, every `_register_lazy_subclade`) via new optional kwargs on each — defaults keep every existing `--taxon`/`--input` caller byte-identical. `--genome-dir --megatree --skip-qc` together is rejected by `main()` (every subclade/backbone build QCs its own genomes in this path; no `skip_qc` passthrough exists yet).

**Geographic figures** (`--geo`, opt-in): after a tree is published, generates a country-annotated tree + world map + patristic-distance pair lines via `python_scripts/annotate_tree_by_country.py` (run as a subprocess by `_generate_geo_figures`) into `<results-dir>/figures/`. Full flags, outputs, and query-highlighting behavior: `readmes/README_ADVANCED_FEATURES.md` §12. Implementation notes not covered there:
- Hooked at all four tree-publish points (megatree, local-genome, taxon-create, batch aggregation). `_generate_geo_figures` itself is non-fatal; only a missing `matplotlib` is checked fatally, at startup (`_validate_dependencies`).
- Map rendering is plain matplotlib against a bundled, public-domain Natural Earth 110m dataset (`python_scripts/data/`) — deliberately not cartopy/geopandas, which pull in proj/shapely/pyproj and fetch coastline data over the network, breaking offline/container runs. Country-name resolution against that dataset (`resolve_country_latlon()`) runs as a fixed-point multi-pass loop, not a single pass, since a `"country"`-kind alias can point at another alias rather than straight at a Natural Earth name.
- `--genome-dir` mode has no NCBI metadata, so `--geo-tsv genome_id<TAB>location` is mandatory there (hard `parser.error` in `main()`).
- The BioSample `lat_lon` attribute does not yet feed these figures (`--geo-tsv`/NCBI `geo_loc_name` only) — see **Geo-provenance checks** below, which does extract and use it, but only for flagging, not figure rendering.

**Geo-provenance checks** (automatic during NCBI download for `--taxon`/`--input` novel-taxon builds; flag-only, never drops or relabels a genome): flags a submitter typing their own institution's location into a sample's geolocation fields. Full column contract and flag semantics: `readmes/README_ADVANCED_FEATURES.md` §13. Implementation notes not covered there:
- `utils/gather_filter_asms.sh`'s `annotate_geo_provenance` step (`python_scripts/geo_provenance.py annotate`) is non-fatal — a failure warns and leaves the unannotated 13-column `all_asm_acc_metadata` in place.
- `country_at()` reuses the same bundled NE 110m dataset the `--geo` figures use (ray-casting + nearest-vertex fallback within 50 km for coastal points).
- `--geo-ror`/`--ror-lookup` is independent of `--geo` and deliberately does not imply it — it affects metadata gathering, not figure rendering.
- Not yet done: `--ror-lookup`/`--geo-ror` resolution only feeds the flag columns, not `resolve_submitter_country`'s offline text scan; a live network smoke test against real NCBI/ROR data.

**Viral annotation** (`--annotator viral`, opt-in; only `--genome-dir` is supported): annotates raw eukaryotic viral assemblies with NCBI **VADR** (model-based classification) instead of Prodigal, falling back to **VIGOR4** (auto-detected reference DB) when no VADR model library classifies the genome. `--taxon`/`--input` novel-taxon builds and `--megatree` reject `--annotator viral` outright — not a missing feature, a mismatched combination: those paths detect viral taxa and get NCBI's existing annotations directly instead (see **Viral NCBI download** below).
- `python_scripts/viral_annotate.py` does the per-genome work (VADR first, VIGOR4 fallback, header normalization to `<g>_<n>` matching `$trans`/`$prots`); `script_lib/arg_parse.sh`/`arg_parse_addem.sh` already exposed `--annotator`/`--vadr-models`/`--vigor-db`/`--vigor-db-dir` on both `OrthoPhyl.sh` and `ReLeaf.sh` before the wrapper layer existed.
- **Wrapper threading** (`orthophyl_pipeline_wrapper.py`): `--annotator`/`--vadr-models`/`--vigor-db`/`--vigor-db-dir` flow through `_run_orthophyl` for every OrthoPhyl build. `_run_releaf` takes its `annotator` argument from the **target database's own** recorded value (never this run's CLI default) — a query ReLeaf'd onto a pre-existing viral DB is always re-annotated viral-style, mirroring how `tree_method`/`tree_data` are read from the matched DB rather than `self.*`. `--vadr-models`/`--vigor-db-dir` on the CLI still describe this run's own tool/DB locations on disk (not per-database), so `_run_releaf` reuses `self.vadr_model_libs`/`self.vigor_db`/`self.vigor_db_dir` when the matched DB is viral and raises a clear `RuntimeError` if they weren't supplied.
- **DB metadata chain**: `OP_database_tool.py --annotator prodigal|viral` writes `annotator` into `database_config.json` (`.get()` default `"prodigal"`, backward-compatible); `assembly_router.py`'s `_load_single_database`/`_route_to_releaf` read it back and surface it on the routing-decision dict, which the wrapper's `_parse_routing_results` carries into `releaf_batch` entries.
- `--genome-dir --annotator viral` auto-forces `self.skip_qc = True` (logged) since CheckM2 is bacteria-only — no need to also pass `--skip-qc`.
- **Not yet done**: the `viral_annot` micromamba env (VADR, exonerate, openjdk, diamond) and VIGOR4 release unpack in `Singularity.OP.v3.3.recipe`, and `tests/integration/test_viral_annotation.py` with committed `TESTER/viral_test/` fixtures — both need the actual VADR/VIGOR4 binaries to build/capture against, unavailable in a dev sandbox.

**Viral NCBI download** (`--taxon`/`--input` novel-taxon builds, automatic — no flag needed): a taxon whose NCBI taxid chain includes 10239 (Viruses) is detected via `is_viral_taxid` in `utils/taxon_assembly_gatherer.py` (`TaxonAssemblyGatherer.is_viral`, checked by rank-independent taxid membership, not a rank-name string — NCBI's root rank label for Viruses has changed over time) and routed to a completely separate download path instead of the bacteria-only CheckM2/assembly pipeline. `render_gtdb_lineage(..., is_viral=True)` forces the rendered lineage's domain slot to `d__Viruses` regardless of what NCBI's own rank names say.
- `utils/gather_viral_annots.py` downloads NCBI's **existing CDS+protein annotations** for the taxon (`datasets download virus genome taxon ... --include cds,protein,annotation`) rather than raw assemblies — there is nothing to annotate, so VADR/VIGOR4 never run on this path. It splits the pooled `cds.fna`/`protein.faa` into one `annots_nucls/<id>.fna` + `annots_prots/<id>.faa` pair per genome, grouping multi-segment records by isolate name when present (segmented-virus packages carry no BioSample), else by accession — meaning a genome's filename stem is not always a literal NCBI accession. `resolve_cds_strand` brute-forces frame×orientation against `protein.faa` per CDS because `cds.fna` does not reliably reverse-complement minus-strand/ambisense segments and some ranges start short of a ribosomal-frameshift junction (both confirmed against real NCBI Virus packages).
- `orthophyl_pipeline_wrapper.py`: `_download_genomes` branches to `_download_viral_genomes` for a viral taxon; `_run_taxon_create_mode_viral` is the dedicated create-mode path (`-a <nucls>,<prots>` instead of `-g` to `OrthoPhyl.sh`, `-n` computed from the CDS-file count, `annotator='viral'`/`input_type='annotations'` recorded on the new DB). Oversized viral taxa reuse the normal subsampling flow; `_stage_selected_viral_prots` symlinks each selected genome's matching `.faa` alongside its `.fna`. `--megatree` for viral taxa is rejected with an explicit error (out of scope).
- **Update mode** (`_run_taxon_update_mode_viral`, dispatched when the existing DB's `annotator == 'viral'`): there's no bacteria-style `assembly_summary.txt` for viruses, so this re-downloads the full current NCBI Virus set via the same `_download_viral_genomes` used at creation time and diffs by `annots_nucls/*.fna` filename stem against the DB's stored `assembly_accessions` — the creation-time ID scheme, not a raw NCBI accession, since segment/isolate grouping means they can differ. New genomes are already NCBI-annotated, so they're added via `_run_releaf`'s `annotation_dirs` param (`-a`, no VADR/VIGOR4).
- `_run_releaf` gained `annotation_dirs: Optional[Tuple[Path, Path]]` (mutually exclusive with `input_genomes`) precisely to let pre-annotated genomes flow through ReLeaf's pre-existing `-a` flag.
- **Router**: `assembly_router.py`'s `_route_to_orthophyl` detects `is_viral` (via `GTDBTaxonomy.get_rank('d') == 'Viruses'`) for a novel query and branches its advisory text — Step 1 points at `utils/gather_viral_annots.py` instead of `gather_filter_asms.sh`, Step 2's `OrthoPhyl.sh` invocation uses `-a` instead of `-g`. `is_viral` is surfaced on the routing-decision dict.
- **Not yet done**: the docs pass beyond this file, and a slow/network-marked end-to-end integration test (`--taxon "<small viral species>"` through to a tree, then a raw query routed to ReLeaf). CheckV-based QC for viral genomes remains unimplemented.

## Repo hygiene warnings

- Historically, files were **versioned by copy** rather than git (`.v2`, `.cmd_out3`, `_v2`, `.bak`, `.save`, `.1.1` suffixes). The live scripts have since been renamed to canonical names and the dead copies removed. If you find such suffixed twins reappearing, commit over the original instead — git is the version history. `dev_plans/` still references the old names historically; don't treat those as live paths.
- Large binaries live in the tree: `*.sif` Singularity images (multi-GB), `GPLv3.pdf`, `ncbi_dataset.zip`. `.sif` is gitignored. Do not commit new large artifacts.
- Test outputs under `TESTER/` (e.g. `FULLTEST_OUT.*`) are gitignored ephemera.

## Containers

Apptainer/Singularity is the only supported container path (the Docker recipe has been removed). `Singularity.OP.v3.3.recipe` is current; `Singularity.OP.v3.1.0.recipe`/`Singularity.OP2.2.1.recipe` are historical and not maintained. It uses micromamba with a base env (orthofinder, iqtree, fasttree, hmmer, mash, prodigal, trimal, raxml, ete3, pytest) and a separate `gather_genomes` env (checkm2, entrez-direct, ncbi-datasets-cli). Conda env spec: `orthophyl_env.2.2.1.yml`. Building requires `--build-arg checkm2_db_path=/path/to/uniref100.KO.1.dmnd` (no network fallback). The default `%runscript` just drops into a shell — use `apptainer exec <image> bash /opt/gits/OrthoPhyl/OrthoPhyl.sh ...`, not `apptainer run <image> <args>`. Full guide: `readmes/README_APPTAINER.md` (covers building, bind-mount semantics, running tests in-container, and moving databases between host and container).

**Megatree database layout**: a subclade/backbone DB (`is_subclade`/`is_backbone` + `parent_taxon`) nests under its parent taxon directory — `<database-dir>/<Taxon>/<Taxon>_db` (backbone), `<database-dir>/<Taxon>/<Taxon>_N_db` (subclades). Non-megatree DBs stay flat at `<database-dir>/<name>_db`. Every path a database persists (symlinks, sketch/members files, `orthophyl_runs.tsv` entries) is stored relative to a nearby anchor so a database built under one bind mount (host or container) works unmodified under another, as long as the bind is an identity bind (host path == container path). Databases predating this, or built with a mismatched non-identity bind, can be repaired with `python assembly_router/OP_database_tool.py --migrate [--strip-prefix PREFIX] [--dry-run]` — infers the old mount prefix, nests megatree DBs, rewrites paths relative to the new location, relativizes `orthophyl_runs.tsv`, regenerates the index, and reports (without deleting) dangling symlinks. Idempotent.

Other readmes live in `readmes/`: `README_APPTAINER.md` (containers), `README_USAGE.md` (worked examples for all three modes), `README_ADVANCED_FEATURES.md`, `README_PIPELINE_PHASES.md` (includes the routing workflow diagram), `README_TECHNICAL_REFERENCE.md`, `README_TROUBLESHOOTING.md`.
