# OrthoPhyl Wrapper & Assembly Router — Unit Test Plan

## Scope

This plan covers the **pipeline-wrapper / assembly-router layer** documented in
`README.v2.md` and `PIPELINE_SUMMARY.md`. It is **complementary to** the existing
`UNIT_TEST_PLAN.md`, which covers the *core* OrthoPhyl scientific pipeline
(`ANI_genome_picking.py`, `OG_sco_filter.py`, `Newick2FastTreeConstraints.py`, and the
Bash functions in `functions.sh` / `functions_addem.sh`). There is no overlap: this
document tests orchestration, routing, and database management — not phylogenetics.

Components in scope:

| # | Component | File | Status |
|---|-----------|------|--------|
| 1 | Pipeline wrapper (orchestrator) | `orthophyl_pipeline_wrapper.v2.py` | exists |
| 2 | Multi-database assembly router | `assembly_router/assembly_router_multi.cmd_out3.py` | exists |
| 3 | Hierarchical database creator | `assembly_router/create_hierarchical_database_v2.py` | exists |
| 4 | NCBI assembly stats | `utils/ncbi_assembly_stats_v2.py` | exists |
| 5 | Genome gather/filter | `utils/gather_filter_asms.sh` | exists |
| 6 | ReLeaf version manager | `assembly_router/add_releaf_version.py` | exists |
| 7 | Taxon assembly gatherer | `utils/taxon_assembly_gatherer.py` | exists |

---

## ⚠️ Pre-existing issues surfaced while writing this plan

These were defects the tests below are designed to catch. **B1–B4 are now FIXED**
(commit pending); each has a regression test guarding it. Descriptions kept for context.

### B1. `_run_releaf` returned before ever running ReLeaf — ✅ FIXED (`orthophyl_pipeline_wrapper.v2.py:475`)
```python
        if self.dry_run:
            logger.info(f"  [DRY RUN] Would run ReLeaf for {database_name}")
        return          # <-- WAS indented at 'if' level → ALWAYS returned
```
The `return` sat at the same indent as `if self.dry_run`, so the method returned
unconditionally. **ReLeaf was never executed**, and everything below (log handling,
version creation) was dead code in real runs. Fixed by indenting the `return` into the
`if self.dry_run:` block.
Regression test: `test_wrapper_batch.py::TestReleafPhase::test_releaf_actually_invoked_when_not_dry_run`.

### B2. `_create_releaf_version` returned before doing anything — ✅ FIXED (`:529`)
```python
        if not db_dir:
            logger.warning(f"  ⚠ Could not find database for {database_name}")
        return          # <-- WAS always returning, even when db_dir was found
```
Fixed by indenting the `return` into the `if not db_dir:` block.

### B3. Second unconditional `return` in `_create_releaf_version` — ✅ FIXED (`:542`)
```python
        if self.dry_run:
            logger.info(f"  [DRY RUN] Would create new database version")
        return          # <-- WAS always returning
```
Same fix. B2/B3 regression test:
`test_wrapper_batch.py::TestReleafPhase::test_version_creation_runs_when_db_found`.

### B4. Wrapper ↔ `TaxonAssemblyGatherer` API mismatch — ✅ FIXED (taxon mode)
The wrapper's taxon-mode paths called an interface the class did **not** implement.
Comparing `orthophyl_pipeline_wrapper.v2.py` against `utils/taxon_assembly_gatherer.py`:

| Wrapper called | Was | Fix |
|---|---|---|
| `TaxonAssemblyGatherer(taxon_name=..., ...)` | ctor kwarg is `taxon=` → `TypeError` | wrapper now passes `taxon=` |
| `gatherer.query_ncbi()` | no such method (only `gather_assemblies()`) | added `query_ncbi()` alias on gatherer |
| `gatherer.download_assemblies(...)` | did not exist | **implemented** on gatherer (FTP fetch + gunzip) |
| `gatherer.get_taxonomy_string()` | did not exist | **implemented** on gatherer (NCBI lineage → GTDB string) |
| `gatherer.rank` | attribute is `.taxon_rank` | wrapper now reads `.taxon_rank` |

**Was:** `_run_taxon_create_mode` / `_run_taxon_update_mode` raised `TypeError`/
`AttributeError` outside `--dry-run` — taxon mode was non-functional. Reconciled by
fixing the wrapper's kwarg/attribute names *and* adding the two genuinely-missing methods
to the gatherer (where the README already documented that behavior).
Regression/contract tests: `test_taxon_assembly_gatherer.py::TestWrapperContract` plus
unit coverage of the two new methods.

Other robustness gaps worth encoding as tests:
- `_download_genomes` counts genomes with `*.fna`/`*.fasta`, but `_run_orthophyl`'s
  total count (`:613`) only globs `*.fna` — `.fasta` query genomes are silently
  undercounted.
- `main()` sets root logger to `DEBUG` on any `-v`, but per-subprocess stdout/stderr
  routing keys off the numeric `verbose` level — worth a test that level→stream mapping
  is correct.
- `add_releaf_version.count_genomes_in_releaf` returns `(0, len(all_genomes))` from a
  `genome_list` (the "new genomes" count is a hardcoded `0` with a `TODO`) — a test
  should document that `added_genomes` is currently always `0` when a `genome_list`
  exists.

---

## Testing frameworks

| Layer | Framework | Rationale |
|-------|-----------|-----------|
| Python units | `pytest` + `pytest-mock` + `pytest-cov` | mock `subprocess.run`, isolate filesystem |
| Bash units | `bats-core` + `bats-mock` | stub `datasets`, `checkm`, `wget`, `statswrapper.sh` |
| Fixtures | `tmp_path` (pytest) / `mktemp -d` (bats) | no network, no real bioinformatics tools |

Install:
```bash
pip install pytest pytest-mock pytest-cov
# bats
git clone https://github.com/bats-core/bats-core.git && (cd bats-core && ./install.sh /usr/local)
```

Proposed layout:
```
tests/
├── unit/
│   ├── test_gtdb_taxonomy.py
│   ├── test_router.py
│   ├── test_db_creator.py
│   ├── test_wrapper_batch.py
│   ├── test_wrapper_taxon.py
│   ├── test_ncbi_stats.py
│   ├── test_gather_filter_asms.bats
│   ├── test_add_releaf_version.py
│   └── test_taxon_assembly_gatherer.py
├── fixtures/
│   ├── databases/                        # fake *_db dirs + database_index.json
│   ├── orthophyl_run/                    # skeleton OP output (HMMs, trees, lists)
│   ├── assemblies.tsv
│   ├── orthophyl_runs.tsv
│   ├── routing_decision_*.json
│   ├── assembly_summary_refseq.txt       # truncated NCBI header + few rows
│   └── taxdump_mini/                      # tiny nodes.dmp + names.dmp
└── conftest.py                            # shared fixtures & fakes
```

A note on imports: the module filenames contain dots (`...cmd_out3.py`,
`...wrapper.v2.py`), so they can't be `import`ed normally. Load them in `conftest.py`
via `importlib.util.spec_from_file_location` and expose the module objects as fixtures.

---

## 1. `GTDBTaxonomy` (shared class — highest ROI)

`GTDBTaxonomy` is defined **twice** (in the router and in the db creator) with near-
identical logic. It's the purest, most testable unit and underpins all routing. Test
both copies (parametrize over the two module objects) and consider deduplicating into a
shared module afterward.

**Suite:** `tests/unit/test_gtdb_taxonomy.py`

| Test | Input | Expected |
|------|-------|----------|
| `test_parse_full_lineage` | `d__Bacteria;p__X;c__Y;o__Z;f__F;g__G;s__S` | all 7 ranks populated |
| `test_parse_empty_species` | `...;g__Escherichia;s__` | `get_rank('s')` is `None`, `get_rank('g')=='Escherichia'` |
| `test_parse_strips_whitespace` | `d__Bacteria; p__X` (space after `;`) | still parses `p` correctly |
| `test_parse_ignores_malformed` | `garbage;p__X` | only `p` captured, no crash |
| `test_get_most_specific_rank` | genus defined, species empty | returns `'g'` |
| `test_get_most_specific_rank_none` | all empty | returns `None` |
| `test_get_rank_name` | `'o'` | `'order'` |
| `test_is_within_clade_match` (router only) | query genus in same family, clade at `f` | `True` |
| `test_is_within_clade_mismatch` (router only) | different phylum | `False` |
| `test_is_within_clade_higher_query` (router only) | query only to order, clade at genus | `False` (can't be within a deeper clade) |
| `test_get_taxonomy_string_at_rank` (router only) | rank `'f'` | string truncated at family, empties rendered `x__` |

Edge cases: quotes around the string (router strips them upstream, class does not),
duplicate rank prefixes, and a leading header token.

---

## 2. Assembly Router — `assembly_router_multi.cmd_out3.py`

**Suite:** `tests/unit/test_router.py`

Fixtures: a temp `database_dir` with two fake `*_db` dirs (a family-level and a genus-
level database sharing a lineage) plus a `database_index.json`. Each `database_config.json`
carries `clade_taxonomy`, `clade_rank`, `clade_rank_name`, `n_genomes`,
`available_tree_methods`, `available_data_types`.

### 2.1 `_load_databases` / `_load_single_database`
1. `test_load_from_index` — loads exactly the databases listed in `database_index.json`.
2. `test_load_by_glob_fallback` — no index file present → scans `*_db` dirs instead.
3. `test_empty_dir_raises` — no databases → `ValueError("No databases found")`.
4. `test_skips_index_entries_with_missing_dir` — index points at a non-existent dir → skipped, not crashed.

### 2.2 `find_matching_databases`
1. `test_single_match` — query within one family DB → 1 match.
2. `test_most_specific_first` — query matches both family & genus DBs → genus (higher
   specificity index) sorts first.
3. `test_no_match` — novel phylum → empty list.

### 2.3 `route_assembly` → ReLeaf (`_route_to_releaf`)
1. `test_routes_to_releaf_on_match` — decision `pipeline == 'ReLeaf'`, correct
   `matched_database`, `matched_rank`.
2. `test_prefers_iqtree_cds` — config offers `[fasttree, iqtree]` / `[PROT, CDS]` →
   picks `iqtree` + `CDS`.
3. `test_falls_back_when_preferred_absent` — config offers only `fasttree` / `PROT` →
   picks those, appends the "non-default" note to the command.
4. `test_copies_assembly_to_releaf_input` — asserts `releaf_input/{id}.fna` created.
5. `test_writes_decision_json_and_summary` — both `routing_decision_{id}.json` and
   `routing_summary_{id}.txt` exist and are valid.

### 2.4 `route_assembly` → OrthoPhyl (`_route_to_orthophyl`)
1. `test_novel_taxon_routes_to_orthophyl` — `pipeline == 'OrthoPhyl'`.
2. `test_download_rank_uses_most_specific` — genus defined → `download_rank == 'genus'`,
   `download_value == <genus>`.
3. `test_species_empty_falls_back_to_genus` — `s__` empty → downloads at genus.
4. `test_genus_empty_falls_back_to_family` — `g__` and `s__` empty → downloads at family.
5. `test_command_includes_gather_script_when_provided` — gather script path present →
   emitted command calls it with `download_value output_dir threads`.
6. `test_command_manual_instructions_when_no_script` — no script → manual-download text.
7. `test_download_taxonomy_string_truncated_at_rank`.

### 2.5 `batch_route`
1. `test_batch_mixed` — TSV of matched + novel → correct split counts.
2. `test_batch_header_detection` — first line starting with `assembly`/`#` is skipped;
   a data first line is **not** skipped (`f.seek(0)` path).
3. `test_batch_tilde_expansion` — `~/g.fna` path is expanded via `expanduser()`.
4. `test_batch_strips_quotes_from_taxonomy`.
5. `test_batch_skips_short_lines` — `<2` fields → warning, row skipped, no crash.
6. `test_batch_summary_written` — `batch_routing_summary.txt` has correct DB grouping and totals.

### 2.6 `main` / CLI
1. `test_requires_taxonomy_with_assembly` — `--assembly` without `--taxonomy` → parser error.
2. `test_mutually_exclusive_assembly_batch`.
3. `test_exit_code_on_init_failure` — empty DB dir → returns `1`.

---

## 3. Database Creator — `create_hierarchical_database_v2.py`

**Suite:** `tests/unit/test_db_creator.py`

Fixture `orthophyl_run/`: a skeleton OP output containing
`OG_alignmentsToHMM/hmms_final/*.hmm`, `FINAL_SPECIES_TREES/*.iqtree.treefile`,
`phylo_current/AlignmentsProts.trm/*.fa`, and a `genome_list` file.

### 3.1 `validate_orthophyl_run`
1. `test_valid_run` — full skeleton → `valid == True`, `has_hmms`, `has_trees` set.
2. `test_hmm_location_priority` — HMMs only in the 2nd candidate location → still found.
3. `test_missing_hmms_warns` — no HMM dir → `has_hmms False`, warning recorded.
4. `test_prefers_iqtree_tree` — dir has both fasttree & iqtree trees → picks iqtree file.
5. `test_genome_count_from_genome_list` — counts non-comment, non-blank lines.
6. `test_genome_count_fallback_to_faa` — no list files → counts `annots_prots/*.faa`.
7. `test_invalid_when_no_genomes` — HMMs present but 0 genomes → `valid == False`.
8. `test_missing_dir_raises` — `FileNotFoundError`.

### 3.2 `database_exists` / `get_existing_databases`
1. `test_detects_existing_db` — `{safe_name}_db/database_config.json` present.
2. `test_safe_name_sanitization` — clade `"Foo Bar/Baz"` → `Foo_Bar_Baz_db`.
3. `test_get_existing_ignores_bad_config` — unreadable/invalid JSON is skipped silently.

### 3.3 `create_database_for_run`
1. `test_creates_config_tree_genomelist_readme` — all four artifacts written.
2. `test_symlink_to_source` — `orthophyl_run` symlink points at resolved source.
3. `test_symlink_fallback_copies_hmms` — simulate `symlink_to` raising `OSError`
   (monkeypatch) → HMMs copied into `orthophyl_run/hmms` instead.
4. `test_placeholder_tree_when_no_tree` — no tree file → placeholder `phylogeny.nwk`.
5. `test_raises_on_existing_without_force` — `FileExistsError`.
6. `test_force_rebuilds` — `--force` removes and recreates.

### 3.4 `parse_input_table`
1. `test_parses_three_columns`.
2. `test_skips_header_and_comments`.
3. `test_skips_rows_with_fewer_than_3_fields` — warns, continues.
4. `test_skips_empty_fields`.

### 3.5 `create_master_index` + `main`
1. `test_master_index_contents` — `n_databases`, per-DB rank fields correct.
2. `test_summary_txt_written`.
3. `test_update_mode_skips_existing` — existing clade in TSV is skipped; index still
   refreshed to include it.
4. `test_update_no_new_still_updates_index`.
5. `test_update_and_force_together_errors` — parser rejects the combination.
6. `test_exit_code_reflects_failures` — a failing entry → return `1`.

---

## 4. Pipeline Wrapper — `orthophyl_pipeline_wrapper.v2.py`

This is the orchestrator. **All external calls (`subprocess.run`) and the two child
scripts must be mocked** — unit tests must never invoke ReLeaf/OrthoPhyl/gather/NCBI.
Patch `subprocess.run` to a fake returning `returncode=0` and record `cmd` argv.

**Suites:** `tests/unit/test_wrapper_batch.py`, `tests/unit/test_wrapper_taxon.py`

### 4.1 Construction & path wiring
1. `test_output_subdirs_computed` — `routing_dir`, `releaf_dir`, `orthophyl_dir`,
   `results_dir`, `logs_dir`, `checkpoint_dir` derived from `output_dir`.
2. `test_script_paths_relative_to_wrapper` — `assembly_router`, `database_creator`,
   `releaf_versioner`, `orthophyl_script`, `releaf_script` resolve under `script_dir`.
3. `test_taxon_mode_flag` — `taxon` set → `taxon_mode == True`.

### 4.2 `_phase_initialization`
1. `test_creates_all_dirs`.
2. `test_resume_skips_when_flag_present`.
3. `test_missing_scripts_raise` — remove `OrthoPhyl.sh` fixture → `FileNotFoundError`
   listing the missing script(s).
4. `test_missing_gather_script_downgrades_to_none_with_warning` (not fatal).
5. `test_no_db_and_no_orthophyl_runs_raises`.
6. `test_creates_initial_db_when_tsv_given` — calls db creator with `--input/--output-dir`.

### 4.3 `_validate_dependencies`
1. `test_all_present_ok`.
2. `test_reports_each_missing_script`.

### 4.4 `_phase_routing` + `_parse_routing_results`
1. `test_builds_router_command` — argv contains `--batch`, `--database-dir`,
   `--output-dir`, `--threads`; adds `--gather-filter-script` only when set.
2. `test_verbose_stream_routing` — v0 → stdout+stderr to logfile; v1 → stderr to file,
   stdout inherited; v2 → both inherited. (Assert on the `subprocess.run` kwargs.)
3. `test_parse_releaf_decisions` — reads fixture `routing_decision_*.json`, builds
   `releaf_batch` with `tree_method`/`tree_data` defaults.
4. `test_parse_orthophyl_decisions_grouped_by_download_value`.
5. `test_routing_failure_raises` — mocked non-zero return → `RuntimeError`.
6. `test_resume_loads_prior_results`.

### 4.5 `_phase_releaf` / `_run_releaf`  ⚠️ (targets bug **B1**)
1. `test_groups_by_database`.
2. `test_copies_assemblies_into_input_genomes`.
3. **`test_releaf_actually_invoked_when_not_dry_run`** — assert `subprocess.run` is
   called with the `ReLeaf.sh` argv. *Fails on current code (B1).*
4. `test_releaf_command_shape` — `--store <db>/orthophyl_run`, `--input_genomes`,
   `-t`, `--tree_method`, `--TREE_DATA`.
5. `test_dry_run_does_not_invoke_releaf` — no `subprocess.run` for ReLeaf.
6. `test_releaf_nonzero_raises` — mocked failure → `RuntimeError` w/ log path.
7. `test_checkpoint_written_per_database`.
8. `test_resume_skips_completed_database`.

### 4.6 `_create_releaf_version`  ⚠️ (targets bugs **B2/B3**)
1. **`test_version_created_when_db_found`** — with `add_releaf_version.py` present and a
   matching `*_db`, assert the versioner subprocess is invoked. *Fails on current code.*
2. `test_warns_when_db_not_found`.
3. `test_skipped_when_versioner_missing`.
4. `test_dry_run_no_version_subprocess`.

### 4.7 `_phase_orthophyl` (novel taxa)
1. `test_download_then_add_queries_then_run_then_db` — ordering of the 4 stages.
2. `test_query_genomes_copied_into_genomes_to_keep`.
3. `test_skip_download_requires_existing_genomes` — `--skip-download` with no
   `genomes_to_keep/` → `FileNotFoundError`.
4. `test_orthophyl_command_shape` — `-g`, `-s`, `-t`, `-p iqtree`, `-o CDS`.
5. `test_resume_skips_completed_stages` (download / orthophyl / database).
6. `test_fasta_vs_fna_count_consistency` — put a `.fasta` query in the set; assert the
   "total genomes" count matches what was actually copied. *(Exposes the `*.fna`-only
   glob at `:613`.)*

### 4.8 `_download_genomes` / `_verify_download_complete`
1. `test_gather_command_argv` — `[script, taxon, out, threads]`.
2. `test_low_ram_appends_reduced_tree`.
3. `test_use_bbmap_appends_flag` (and that it takes precedence over `--low-ram`).
4. `test_no_gather_script_warns_and_returns`.
5. `test_zero_genomes_after_download_raises` — empty `genomes_to_keep/` → `RuntimeError`
   with the actionable message.
6. `test_writes_download_complete_marker`.
7. `test_verify_requires_checkpoint_marker_and_files` — all three conditions; each
   missing piece → `False`.

### 4.9 `_run_orthophyl` + `_verify_queries_in_tree`
1. `test_verifies_each_query_in_treefile` — treefile containing/omitting an id →
   correct ✓/⚠ logging and `all_found`.
2. `test_missing_treefile_warns_not_crashes`.

### 4.10 `_create_database_entry`
1. `test_appends_to_orthophyl_runs_tsv`.
2. `test_calls_db_creator_update_mode`.
3. `test_nonzero_raises`.

### 4.11 `_phase_aggregation` + `_generate_summary_report`
1. `test_collects_releaf_trees` — copies `ReLeaf_results/phylogeny_with_new_genomes.nwk`.
2. `test_collects_orthophyl_trees` — copies `FINAL_SPECIES_TREES/SCO_strict.CDS.iqtree.treefile`.
3. `test_summary_report_lists_both_routes`.
4. `test_status_phase_counts`.

### 4.12 Taxon mode (`test_wrapper_taxon.py`)  ⚠️ (targets bug **B4**)
> `utils/taxon_assembly_gatherer.py` now exists, but the wrapper calls an interface it
> doesn't implement (see **B4**). Two complementary approaches:
> - **Contract tests** — inject a stub `TaxonAssemblyGatherer` via
>   `sys.modules['taxon_assembly_gatherer']` exposing exactly what the wrapper calls
>   (`taxon_name=` kwarg, `query_ncbi()`, `download_assemblies()`, `get_taxonomy_string()`,
>   `.rank`). These isolate wrapper logic from the gatherer.
> - **Integration guard** — one test that imports the *real* class and asserts the
>   wrapper's call signature resolves (constructor kwargs + method names exist). This is
>   the test that currently fails on B4 and should gate any "taxon mode works" claim.
1. `test_existing_db_without_update_flag_errors` — returns `1`, clear message.
2. `test_existing_db_with_update_runs_update_path`.
3. `test_no_existing_db_runs_create_path`.
4. `test_create_mode_queries_downloads_runs_op_creates_db` (all mocked).
5. `test_update_mode_diffs_accessions` — only accessions not already in config are
   downloaded; "up to date" short-circuit when the diff is empty.
6. `test_update_metadata_merges_accessions_and_updates_counts`.
7. `test_dry_run_short_circuits_both_modes`.

### 4.13 `main` / CLI validation
1. `test_taxon_and_input_mutually_exclusive`.
2. `test_requires_one_of_taxon_or_input`.
3. `test_update_existing_requires_taxon`.
4. `test_verbose_count_sets_debug_logging`.
5. `test_run_returns_1_on_uncaught_exception` — `run()` wraps in try/except, writes
   failed status JSON, returns `1`.

### 4.14 Checkpoints & status
1. `test_write_then_check_checkpoint_roundtrip`.
2. `test_save_final_status_writes_json` (and honors dry-run "would save" path).

---

## 5. NCBI Assembly Stats — `ncbi_assembly_stats_v2.py`

**Suite:** `tests/unit/test_ncbi_stats.py`

Mock all network I/O (`urllib.request.urlretrieve`) and use a tiny `taxdump_mini/`
(`nodes.dmp`, `names.dmp`) plus a truncated `assembly_summary_refseq.txt`.

### 5.1 `NCBITaxonomy`
1. `test_load_nodes_and_names` — counts match fixture.
2. `test_names_keeps_only_scientific` — synonyms/common names ignored.
3. `test_get_lineage_walks_to_root` — returns the ranks of interest along the path.
4. `test_get_lineage_cycle_guard` — self-parent / cycle → terminates (uses `visited`).
5. `test_get_rank_name_unknown` — taxid absent → `'Unknown'`.

### 5.2 `NCBIAssemblyStatsV2`
1. `test_download_taxdump_skips_when_present`.
2. `test_extract_only_needed_members` — only `nodes.dmp`/`names.dmp` extracted.
3. `test_parse_summary_finds_header` — single-`#` header located, `##` lines ignored.
4. `test_parse_summary_column_indices` — missing required column → graceful skip.
5. `test_parse_skips_short_rows`.
6. `test_summarize_by_taxonomy_counts` — combined = refseq + genbank per rank.
7. `test_summarize_requires_loaded_taxonomy` — `RuntimeError` if not loaded.
8. `test_generate_rank_table_min_assemblies_filter` — rows below threshold dropped,
   sorted by total desc.
9. `test_summary_report_written`.

---

## 6. `gather_filter_asms.sh` (Bash)

**Suite:** `tests/unit/test_gather_filter_asms.bats`

Stub external tools on `PATH` (`datasets`, `unzip`, `checkm`, `statswrapper.sh`,
`esearch`, `esummary`, `xtract`, `wget`) with fake scripts that emit canned fixtures.
Because the script calls `source ~/.bash_profile` and `conda activate` at the top, tests
should either shim those or refactor the script to guard them behind a testability flag.

### 6.1 Argument & path handling
1. `test_absolute_path_used_verbatim`.
2. `test_tilde_expansion`.
3. `test_relative_path_made_absolute`.
4. `test_reduced_tree_flag_parsed` — sets CheckM `--reduced_tree`.
5. `test_use_bbmap_flag_parsed` — routes to `get_asm_stats` (bbmap) branch.

### 6.2 `filter_NCBI_genomes` (RefSeq/GenBank dedup)
1. `test_prefers_refseq_over_genbank` — accession in both → GCF kept, GCA dropped.
2. `test_genbank_only_kept`.
3. `test_union_and_uniq_lists_correct` — validate the `comm -12/-13/-23` outputs.

### 6.3 Stats reformatting
1. `test_bbmap_stats_reformatted_to_checkM_columns` — placeholders (completeness=100,
   contamination=0, dup=0) inserted; genome-size/N50/GC mapped to correct columns.
2. `test_checkM_aggregation_column_math` — the `awk` duplication-ratio formula.

### 6.4 `filter_asm_by_stats`
1. `test_default_thresholds_from_mean_stddev` — `default` → mean ± 3·SD computed.
2. `test_explicit_thresholds_respected`.
3. `test_assembly_below_min_len_removed`; `above_max_len`; `low_n50`; `gc_out_of_range`;
   `high_contam`; `low_completeness` — one row per predicate.
4. `test_too_few_args_exits` — missing 8th arg → error + exit.

### 6.5 `filter_for_redundancy`
1. `test_multiple_versions_collapsed` — `GCA_X.1` and `GCA_X.2` → keeps highest version.
2. `test_removed_list_respected` — accessions in `final_assemblies_to_remove` excluded
   from `genomes_to_keep/`.

> Integration-level (network-dependent, out of unit scope, mark `@slow`):
> `get_NCBI_genomes` retry/backoff logic and 503 handling — best exercised with a mock
> `datasets` that fails N times then succeeds.

---

## 7. Version manager & taxon gatherer (both now present)

Both scripts exist in the repo. Tests below are keyed to their **actual** interfaces.

### 7.1 `assembly_router/add_releaf_version.py`  (`test_add_releaf_version.py`)

Fixtures: a `{clade}_db/` with a `v1_orthophyl_initial/` version dir (containing
`version_info.json`, `database_config.json`, `genome_list.txt`, and an `orthophyl_run/`
skeleton with `OG_alignmentsToHMM/hmms_final`, `phylo_current/{trimmed_columns,SCO_*,
SpeciesTree/iqtree}`) and a `current → v1_orthophyl_initial` symlink. A fake ReLeaf
output dir with `new_prot_alignments.trm.nm`, `new_CDS_alignments.trm.nm`, and
`new_trees/*.tree*`.

**`get_current_version`**
1. `test_reads_current_symlink` — returns the target of `current`.
2. `test_fallback_to_v1_when_no_current`.
3. `test_raises_when_indeterminate`.

**`get_next_version_number`**
1. `test_first_version_is_1` — no `v*` dirs → 1.
2. `test_increments_from_max` — `v1`,`v2` present → 3.
3. `test_ignores_unparseable_version_dirs`.

**`validate_releaf_output`**
1. `test_valid_with_alignments_and_trees` → `valid True`.
2. `test_missing_alignments_invalid` — only trees present → `valid False`, warning.
3. `test_missing_trees_invalid`.
4. `test_missing_dir_raises`.

**`count_genomes_in_releaf`**
1. `test_counts_from_genome_list` — returns `(0, N)`. *(Documents the hardcoded `0`
   "added" count / open `TODO`.)*
2. `test_counts_from_annots_when_no_list` — returns `(n_faa, None)`.
3. `test_returns_zero_none_when_nothing`.

**`create_composite_orthophyl_dir`** (symlink wiring)
1. `test_hmms_linked_from_base` — `OG_alignmentsToHMM/hmms_final` → base.
2. `test_prot_and_cds_alignments_linked_from_releaf` → `AlignmentsProts.trm`,
   `AlignmentsCDS.trm`.
3. `test_trimmed_columns_and_sco_from_base`.
4. `test_trees_linked_and_addasm_suffix_stripped` — `x.addasm.treefile` → `x.treefile`.
5. `test_iqtree_info_linked_from_base`.
6. `test_missing_base_orthophyl_raises`.

**`create_releaf_version`** (end to end, filesystem-only)
1. `test_creates_version_dir_and_info_json` — `version_info.json` fields:
   `version_number`, `version_type=='releaf'`, `parent_version`, `base_genomes`,
   `added_genomes`, `total_genomes`, `source_path`.
2. `test_base_and_releaf_symlinks` — `base_version`, `releaf_output` links.
3. `test_orthophyl_run_points_to_composite`.
4. `test_current_symlink_repointed_to_new_version`.
5. `test_phylogeny_prefers_iqtree_cds_tree`.
6. `test_genome_list_written_or_copied`.
7. `test_invalid_releaf_output_raises` (ValueError).
8. `test_duplicate_version_name_raises` (FileExistsError).
9. `test_missing_parent_version_info_raises` — no `version_info.json` in current
   version → the `open(...)` at `:243` raises. *(Note: a v1 created by
   `create_hierarchical_database_v2.py` has **no** `version_info.json` — so the very
   first ReLeaf version creation on a freshly-built DB will fail. Encode this as an
   xfail/known-issue test; it's the practical path the wrapper drives.)*

**`list_versions` / `main` / CLI**
1. `test_list_versions_marks_current`.
2. `test_requires_releaf_output_or_list` — parser error otherwise.
3. `test_version_name_override_respected` — `--version-name` used verbatim.
4. `test_wrapper_cli_contract` — the wrapper invokes only `--database-dir` +
   `--releaf-output` (`:533-537`); assert those two flags alone produce a valid run
   (auto-named version), so wrapper and script stay in sync.

### 7.2 `utils/taxon_assembly_gatherer.py`  (`test_taxon_assembly_gatherer.py`)

Mock all network I/O (`urllib.request.urlretrieve`); use `taxdump_mini/` and a truncated
`assembly_summary_refseq.txt` (must include the columns the parser requires:
`assembly_accession, taxid, organism_name, asm_name, assembly_level, seq_rel_date,
ftp_path`).

**`NCBITaxonomy`** (this file's own copy — distinct from the ncbi_stats copy)
1. `test_load_nodes_names_and_name_to_taxid` — reverse index populated (lowercased).
2. `test_resolve_taxon_by_name` — case-insensitive name → taxid.
3. `test_resolve_taxon_by_id` — numeric string present in nodes → itself.
4. `test_resolve_unknown_returns_none`.
5. `test_get_rank_get_name`.
6. `test_get_lineage_cycle_guard`.

**`TaxonAssemblyGatherer`**
1. `test_ctor_resolves_taxid_and_sets_name_rank` — sets `.taxid`, `.taxon_name`,
   `.taxon_rank`.
2. `test_ctor_unresolvable_taxon_raises` (ValueError).
3. `test_ensure_taxonomy_skips_when_present`.
4. `test_parse_summary_filters_by_descendant` — `_is_descendant_of_taxon` keeps only
   rows under the query taxid.
5. `test_since_date_filter` — rows before `since_date` dropped.
6. `test_parse_missing_column_raises`.
7. `test_is_descendant_true_and_false` (incl. cycle guard).
8. `test_remove_redundancy_prefers_refseq` — GCF kept over GCA of same base acc.
9. `test_remove_redundancy_latest_version` — GCA_x.1 vs GCA_x.2 → keeps `.2`.
10. `test_gather_assemblies_combines_and_dedups` (mock both summaries).
11. `test_generate_output_tsv_header_and_rows`.
12. `test_generate_metadata_json_fields` — `source_taxon_name`, `source_taxid`,
    `source_rank`, `assembly_accessions`, `n_assemblies`, `quality_filters`.
13. `test_main_returns_1_when_no_assemblies`.

> **Cross-component contract (B4):** add `test_matches_wrapper_call_signature` asserting
> the wrapper's usage resolves against the real class (ctor kwarg name, method names,
> attribute names). This is the guard that fails today and should be fixed before taxon
> mode is advertised as working.

---

## 8. Shared fixtures (`conftest.py`)

- `load_module(path)` helper — `importlib` loader for dotted filenames.
- `fake_subprocess` — records argv, returns configurable `returncode`; asserts no real
  bioinformatics binary is ever invoked.
- `db_dir` — factory building `{clade}_db/database_config.json` (+ optional index) with
  parametrized ranks / tree methods / data types.
- `orthophyl_run_skeleton` — minimal valid OP output tree (HMMs, trees, alignments,
  genome_list) for the db creator and validators.
- `routing_decisions` — writes ReLeaf and OrthoPhyl decision JSONs into a temp
  `00_routing/`.
- `taxdump_mini` — 5–10 node `nodes.dmp` + `names.dmp`.
- `stub_bin` (bats) — puts fake `datasets`/`checkm`/`wget`/`statswrapper.sh` first on `PATH`.

---

## 9. Coverage targets & CI

| Area | Target |
|------|--------|
| `GTDBTaxonomy`, router decision logic | ≥ 95% |
| db creator, wrapper phase methods | ≥ 85% |
| ncbi stats parsing | ≥ 85% |
| bash functions (function-level) | ≥ 80% |

`.github/workflows/wrapper-tests.yml` stages:
1. **Lint** — `ruff`/`pylint` on the Python scripts, `shellcheck` on `gather_filter_asms.sh`.
2. **Python units** — `pytest tests/unit --cov`.
3. **Bash units** — `bats tests/unit/*.bats`.
4. (nightly) end-to-end dry-run of the wrapper against a fixture database — asserts
   directory structure, routing split, and `pipeline_status.json` without invoking any
   real tool.

Matrix: Python 3.7 / 3.9 / 3.11 (wrapper advertises 3.7+), Ubuntu 22.04.

---

## 10. Implementation priority

1. **Phase 1 (foundation).** `GTDBTaxonomy` (§1), router decision logic (§2.2–2.4),
   db-creator validation (§3.1). Pure logic, no mocking, catches the most bugs per line.
2. **Phase 2 (orchestration).** Wrapper phases with mocked subprocess (§4) — **start
   with §4.5/§4.6 to lock down bugs B1–B3**, then routing/orthophyl phases.
3. **Phase 3 (I/O layers).** NCBI stats (§5), `gather_filter_asms.sh` bats (§6),
   version manager + taxon gatherer (§7) — both scripts exist, so these are
   characterization tests, not TDD stubs. **Add the B4 contract test early** (§7.2 last
   item) since it gates taxon mode.
4. **Phase 4 (reconcile).** Fix the wrapper↔gatherer API mismatch (B4) and the
   first-ReLeaf-version `version_info.json` gap (§7.1 test 9); make the contract/xfail
   tests pass.
5. **Ongoing.** Regression test for every future bug; CI gating (§9).
