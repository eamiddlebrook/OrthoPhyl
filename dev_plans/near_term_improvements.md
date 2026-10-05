# Near-term OrthoPhyl improvements (not yet implemented)

Four items the user wants to address soon. Captured with current-code grounding
(file:line) so they're actionable later without re-deriving context.

## 1. Fail fast when the gather script is missing/non-viable

**Ask**: Check for a viable gather script early. Currently the pipeline goes
through the taxonomy decision-tree stuff first; want it to fail earlier.

**Current behavior**: `--taxon` mode's call chain is `main()` → `run()` →
`_phase_initialization()` → `_validate_dependencies()`
(`orthophyl_pipeline_wrapper.py:570-593`). That check is:
```python
if self.gather_script and not self.gather_script.exists():
    logger.warning(f"Genome download script not found: {self.gather_script}")
    logger.warning("  Will generate manual download instructions instead")
    self.gather_script = None
```
Only `Path.exists()` — no executability check, no smoke-test/`--help` invocation.
No such viability-check pattern exists anywhere else in the repo to mirror.

Taxonomy resolution work happens AFTER this check but BEFORE the gather script is
actually invoked: `_run_taxon_create_mode()` constructs `TaxonAssemblyGatherer(...)`
at `:2447`, whose `__init__` calls `_ensure_taxonomy_database()`
(`utils/taxon_assembly_gatherer.py:234→365`), which downloads a ~50-60MB NCBI
taxdump tarball over HTTPS if not cached, then resolves the taxon's rank/lineage.
Only after all that does the gather script get re-checked (`:2472`) and finally
invoked via `subprocess.run` at `:1600` (inside `_download_genomes` →
`_download_raw`).

**Fix shape**: add a real viability check (exists + `os.access(X_OK)`, maybe a
quick `--help`/dry invocation) in `_validate_dependencies()` or earlier in `run()`,
so a bad/missing gather script fails before the taxdump download and taxonomy
resolution, not after.

## 2. Default `--gather-script` to `./utils/gather_filter_asms.sh`

**Ask**: Bake in the current gather script as the default.

**Current behavior**: no default exists. Constructor: `self.gather_script =
Path(gather_script) if gather_script else None` (`:146`). Argparse has no
`default=` (`:3099-3102`): `parser.add_argument('--gather-script', help='Path to
gather_filter_asms.sh for genome downloading')`. The user must pass
`--gather-script` explicitly every invocation or the pipeline falls back to
manual-download-instructions mode.

**Fix shape**: default to `self.script_dir / "utils" / "gather_filter_asms.sh"`
(`script_dir` is already computed as `Path(__file__).parent` at `:253`) when
`--gather-script` is omitted, rather than requiring explicit opt-in every time.

## 3. Megatree subclades reuse the backbone's HMM gene models only

**Ask**: Make megatree subclades use the precomputed "representative subsample"
clade's (i.e. the backbone's) HMM gene models, similar to ReLeaf — but hold off on
reusing anything else from that run (alignments, trim, evolutionary/tree model),
since those might legitimately differ by subclade (gene content, divergence).

**Current behavior**: `_run_megatree` (`orthophyl_pipeline_wrapper.py:1068-1222`)
builds each subclade independently via `_build_subclade` (its own full
`OrthoPhyl.sh` run, own OrthoFinder orthogroup inference, own `hmms_final/` if it
exceeds its own `-n`), and separately builds a backbone over diverse per-subclade
reps (`_subsample_genomes` + its own `_run_orthophyl`). **No sharing happens at
all today** — each of these OrthoFinder runs is fully independent.

The reusable primitive already exists and does exactly the "assign genes into
existing orthogroups via HMM" job: `HMM_search`
(`script_lib/functions_addem.sh:183-263`) — the same mechanism ReLeaf uses to
place a new query genome. For each `.hmm` file in `hmm_dir` it runs `hmmsearch -T
25 --tblout ... $I $all_prots`, applies a per-OG no-paralog score threshold
(`:241-251`), then pulls matching sequences via `filterbyname.sh`
(`:253-256`). It keys purely on the `.hmm` filename stem as the group label —
no hardcoded dependency on OrthoFinder's `OG\d+` ID format.

The internal analog, `OG_hmm_search` inside `ANI_ORTHOFINDER_TO_ALL_SEQS`
(`script_lib/functions.sh:591-796`, HMM build/search at `:610-656`), is what
produces `$OG_alignmentsToHMM/hmms_final/*.hmm` in the first place — this is the
backbone's output ReLeaf (and this proposed feature) would consume. The DB
created from the backbone run (`:1219-1222`, `is_backbone=True`) is already the
canonical ReLeaf target per CLAUDE.md's "Backbone reps... coherent HMM set" note
— confirming the backbone's HMMs are the right thing to point at.

**Fix shape**: for each subclade build, instead of (or before) its own
independent OrthoFinder orthogroup inference, call `HMM_search` against the
backbone's `hmms_final/*.hmm` to assign the subclade's genes into the backbone's
orthogroup IDs. Then let the subclade run its OWN alignment, trimming, and
tree-model/branch-length steps on those assigned genes — explicitly NOT reusing
the backbone's alignments or evolutionary model, per the user's caveat.

## 4. Scaffold support for an external HMM set (e.g. BUSCO)

**Ask** (related to #3, bigger lift): let OrthoPhyl use an external HMM set (like
BUSCO) to identify "easily identifiable" homologs directly, then route only the
remaining unclassified genes through the standard OrthoFinder path.

**Current behavior**: no existing scaffolding. Only mention anywhere in the repo
is doc text (`README_OrthoPhyl_ReLeaf.md:358`) noting OrthoPhyl "does not robustly
compute SCOs (like BUSCOs)... does not do any modeling to ensure species
tree-like behavior." No code, flags, or TODOs reference BUSCO or external HMMs.

This shares its core mechanism with #3: `HMM_search` doesn't care where its
`hmm_dir` came from — in principle it could point at a directory of
externally-sourced HMMs (e.g. BUSCO's, renamed to `<BUSCO_ID>.hmm`) and run
unmodified to classify "easy" genes by homology search alone.

**Open design question (not solved, needs real design work later)**:
`HMM_search`'s output feeds `ADD_2_ALIGNMENTS`
(`script_lib/functions_addem.sh:265-326`), which expects a matching "old
alignment" per OG basename (`$OLD_alignments/OG*`) to append new hits onto. For
an externally-sourced HMM set there is no pre-existing "old alignment" for a
BUSCO ID — so this downstream step needs either synthetic/seed alignment
stand-ins per external-HMM ID, or a dedicated code path instead of reusing
`ADD_2_ALIGNMENTS` as-is. This is the key unresolved piece, not a detail to paper
over when this gets picked up.

**Fix shape** (sketch only): run the external HMM set against all genomes' predicted
proteins first (homology search, no paralog threshold reused from #3/ReLeaf's
pattern), set aside genes that get classified as single-copy orthologs this way,
then run the standard OrthoFinder pipeline only on each genome's remaining
(unclassified) gene set, and combine both orthogroup sets downstream for
alignment/tree-building.
