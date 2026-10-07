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

## 3. Megatree subclades reuse the backbone's HMM gene models only [DONE]

**Status**: Implemented. New `OrthoPhyl.sh --hmm-assign-dir <path>` flag
(`script_lib/arg_parse.sh`) branches `MAIN_PIPE` to a new
`HMM_ASSIGN_FROM_EXTERNAL` function (`script_lib/functions.sh`) that calls
`HMM_search` (reused verbatim from `script_lib/functions_addem.sh`, now also
sourced by `OrthoPhyl.sh`) against the given external HMM dir, aligns each
assigned OG from scratch with plain `mafft --quiet` (NOT the extend-only
`ADD_2_ALIGNMENTS`/`mafft --add --keeplength`), and drops the result into
`$wd/AlignmentsProts/<OG>.faa` -- the same slot `TRIM`/`SCO_MIN_ALIGN`/
`TRIMAL_backtrans`/`TREE_BUILD` already consume unmodified. OGs with zero
hits are logged to `$wd/hmm_assign_unmatched_OGs.txt` (the seam for item 4
below) rather than erroring.

`orthophyl_pipeline_wrapper.py`'s `_run_megatree` now has an opt-in
`--megatree-hmm-reuse` flag (default off, so plain `--megatree` is
byte-identical to before): when set, `_build_megatree_subclades_hmm_reuse`
QCs every subclade and pools backbone reps FIRST, builds the backbone tree
(producing `hmms_final/`), then builds every subclade passing
`hmm_assign_dir=<backbone's hmms_final/>` into `_build_subclade` ->
`_run_orthophyl` (`--hmm-assign-dir`). Falls back to the default independent
per-subclade-OrthoFinder path with a warning if the backbone never produces
`hmms_final/` (e.g. too few pooled reps to cross `--ani-shortlist`).
Default path's exact call order/behavior is preserved in a sibling
`_build_megatree_subclades_independent` method and covered by the
pre-existing `TestRunMegatree` suite, which still passes unmodified.

Verified end-to-end with real tools (chloroplast genomes): built a real
8-genome backbone (forced through the HMM-building branch via `-n 3`), then
ran a disjoint 6-genome subclade with `--hmm-assign-dir` pointed at the
backbone's `hmms_final/`. Confirmed no `OrthoFinder`/`Results_*` directory
anywhere in the subclade's output, a valid tree was produced, 84/85 backbone
OGs got assigned (1 correctly logged to `hmm_assign_unmatched_OGs.txt`), and
every assigned OG ID in the subclade's `AlignmentsProts/` is a subset of the
backbone's HMM filenames (no foreign IDs). Full existing integration suite
(`tests/integration/test_megatree_lazy.py`, 9 tests, real subprocesses) still
passes with the new code present but the flag unused, confirming the default
path is unaffected.

Original ask (kept for context): Make megatree subclades use the precomputed "representative subsample"
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

**Align-from-scratch gap (shared with item 4): SOLVED as part of item 3.**
`HMM_search`'s pre-existing downstream consumer, `ADD_2_ALIGNMENTS`
(`script_lib/functions_addem.sh:265-326`), only ever *extends* a pre-existing
alignment (`mafft --add --keeplength`, `:316`) — no use here, since there's
no backbone alignment to extend. Item 3's new `HMM_ASSIGN_FROM_EXTERNAL`
(`script_lib/functions.sh`) instead aligns each HMM-assigned OG from scratch
with plain `mafft --quiet` (mirroring `OG_hmm_search`'s own precedent) and
drops it straight into `$wd/AlignmentsProts/<OG>.faa`. This is exactly the
"align fresh, no old alignment required" primitive item 4 also needs — reuse
it directly rather than re-solving.

## 4. Scaffold support for an external HMM set (e.g. BUSCO) [DONE]

**Status**: Implemented, then made the default. New
`OrthoPhyl.sh --hmm-assign-leftover-orthofinder` flag (`script_lib/arg_parse.sh`),
only meaningful alongside `--hmm-assign-dir`, landed first as opt-in
(default off). Per a follow-up request, leftover-OrthoFinder routing is now
**ON by default** whenever `--hmm-assign-dir` is used — the flag itself is
kept as a no-op for explicitness/symmetry, and a new
`--skip-hmm-assign-leftover` flag restores the old log-and-drop behavior
(mirrors the existing `--skip-qc` pattern in this codebase: on by default,
explicit opt-out). The default-flip is a one-line change
(`script_lib/functions.sh`'s `${hmm_assign_leftover_orthofinder:-true}`
check), everything else below is unchanged.

At the wrapper level, `--megatree-hmm-reuse` gets this for free (every
subclade build already passes `hmm_assign_dir`, so it now also gets
leftover-routing with no wrapper change) — a matching
`--megatree-hmm-reuse-skip-leftover` wrapper flag was added so a wrapper
user can opt back out too, threading `--skip-hmm-assign-leftover` through
`_run_orthophyl`/`_build_subclade` for every subclade build in that path.

When set, `HMM_ASSIGN_FROM_EXTERNAL`'s tail (`script_lib/functions.sh`) now:
1. Computes the gene-level leftover set — all `all_prots.nm.fa` names MINUS
   every name that appears in any matched `$new_OG_prots/<OG>.faa` (NOT the
   OG-level `hmm_assign_unmatched_OGs.txt`, which is coarser: an OG can have
   *some* hits without covering every gene in every genome) — via
   `comm -23` on sorted name lists, mirroring the pre-existing
   `script_lib/taxa_missing_OGs.sh` idiom (the only prior art for this shape
   of set-difference anywhere in the repo).
2. Materializes a new per-genome `hmm_assign_leftover_prots/<genome>.faa`
   directory (`filterbyname.sh` against `$prots.fixed/<genome>.faa`, whose
   headers are byte-identical to `all_prots.nm.fa`'s — same transform, same
   loop, in `FIX_PROTS_NAMES` — so no header translation needed). Guards:
   skip genomes with zero leftover genes; abort (log + skip, not error) if
   fewer than 2 genomes have any leftovers at all (OrthoFinder needs ≥2
   inputs).
3. Runs a real `ORTHO_RUN` (reused verbatim) on that leftover pool.
4. New sibling function `REALIGN_LEFTOVER_ORTHOGROUP_PROTS` turns the result
   into alignments — **critically, filtered through `OG_sco_filter.py`
   (threshold 1) first**, rejecting any orthogroup with paralogs, before
   realigning. (Externally-assigned OGs are already paralog-free by
   construction, via `HMM_search`'s own no-paralog score filter — a fresh
   OrthoFinder clustering has no such guarantee, and `SCO_MIN_ALIGN`'s
   `ANI=true` branch only counts total headers per alignment file, not
   distinct genomes, so an unfiltered multi-copy OG silently corrupts SCO
   membership and crashes `catfasta2phyml` downstream on mismatched
   per-genome sequence lengths — caught during manual verification below,
   not a hypothetical.) Surviving OGs are renamed from OrthoFinder's native
   `OG0000001`-style IDs to `OG0_LFT_0000001` before landing in
   `$wd/AlignmentsProts/` — avoiding a silent filename clobber against the
   external HMM set's own OG IDs (a real risk whenever `hmm_assign_dir`
   came from a prior OrthoPhyl run, e.g. item 3's `--megatree-hmm-reuse`,
   since OrthoFinder always numbers fresh runs from `OG0000001` with no
   per-run salt), while the `OG0_LFT_` prefix still satisfies
   `SCO_MIN_ALIGN`'s tighter `$alignment_dir/OG0*` glob (stricter than the
   bare `OG*` every other downstream consumer uses).
5. `GET_OG_NAMES` runs as before, now seeing both ID namespaces.

Verified end-to-end with real tools (chloroplast genomes, extending item
3's own scenario): built an 8-genome backbone, took only 40 of its 85
`hmms_final/` HMMs as a deliberately-partial external set, ran a disjoint
6-genome subclade against it with the new flag. Confirmed: 386 leftover
genes found; a real `ORTHO_RUN` executed on the 6-genome leftover pool
(`hmm_assign_leftover_prots/OrthoFinder/Results_ortho/` materialized); 70
single-copy `OG0_LFT_*` orthogroups survived SCO-filtering (down from 85 raw
OrthoFinder OGs — several rejected for paralogs, confirming the filter
engages); 40 external OGs; zero filename collisions between the two
namespaces; both `SCO_strict`/`SCO_2`-derived trees built successfully.
Also reran the identical scenario with the flag OMITTED and confirmed
byte-identical-to-item-3 behavior (exactly 40 OGs, zero `OG0_LFT_*` entries,
no leftover directory created at all). Full Python unit suite (357 tests)
unaffected, as expected for a pure bash-layer change.

**Non-goals, deliberately out of scope for this pass**: no
`orthophyl_pipeline_wrapper.py` CLI flag (pure `OrthoPhyl.sh`-level
capability for now, same staged approach item 3 itself used); no
BUSCO-specific scaffolding (fetching/converting a BUSCO HMM set) —
`--hmm-assign-dir` already accepts any directory of `<id>.hmm` files
regardless of origin, this item only changed what happens to *unmatched*
sequences.
