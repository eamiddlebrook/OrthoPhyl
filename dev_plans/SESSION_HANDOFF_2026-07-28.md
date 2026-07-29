# Session Handoff: OrthoPhyl Project State (2026-07-28)

## Purpose & How to Use This Document

This document provides a comprehensive snapshot of the OrthoPhyl project state as of July 28, 2026. It's designed to allow a fresh AI assistant (Opus/Sonnet) to quickly understand:

1. **What was just completed** (integration test fixes)
2. **Current project structure** (test suite, main components)
3. **In-flight work** (wrapper v2 fixes, taxon mode, router v3)
4. **Key domain facts** (test data, char-stripping rules, CLI flags)
5. **Where to look next** (authoritative docs, open TODOs)

**For a new AI instance:** Read sections 1-4 to get oriented, then jump to section 6 for in-flight work details. Use section 7 as a quick reference for file locations.

---

## 1. Project Overview

### What is OrthoPhyl?

OrthoPhyl is a phylogenomics pipeline that takes genome assemblies and produces phylogenetic trees. The workflow:

1. **Annotate** genomes (Prokka/Prodigal)
2. **Orthologize** proteins (OrthoFinder)
3. **Align** orthogroups (MAFFT)
4. **Trim** alignments (trimAl)
5. **Concatenate** single-copy orthologs (SCOs)
6. **Build tree** (IQ-TREE or FastTree)

### Three Main Products

1. **`OrthoPhyl.sh`** — Main pipeline (de novo phylogeny from genomes)
2. **`ReLeaf.sh`** — Add new genomes to an existing OrthoPhyl run (HMM-based placement)
3. **`orthophyl_pipeline_wrapper.v2.py`** — Database manager + router (decides ReLeaf vs. OrthoPhyl based on taxonomic distance)

### Repository Info

- **Remotes:**
  - `gitea: gitea@bio-gitea.lanl.gov:earlm/OrthoPhyl.git`
  - `origin: git@github.com:eamiddlebrook/OrthoPhyl.git`
- **Current commit:** `48a25e19f322e925b3f552bab198ed1f5132ed0e`
- **Working directory:** `/home/earlm/gits/OrthoPhyl`

### Key Documentation Files

- `README.md` / `README.v2.md` — User-facing docs
- `PIPELINE_SUMMARY.md` — Pipeline architecture overview
- `INTEGRATION_TEST_PLAN.md` — Integration test design + rationale
- `WRAPPER_UNIT_TEST_PLAN.md` — Unit test design for wrapper
- `TAXON_MODE_GUIDE.md` — User guide for taxon-based routing
- `Taxon_mode_implementation_notes.md` — Developer notes on taxon mode
- `work_plan_7_23_26.md` — Detailed plan for wrapper v2 fixes (3 critical bugs)
- `RELEAF_DATABASE_DESIGN.md` — Hierarchical database design for ReLeaf versioning

---

## 2. Environment & How to Run

### Conda Environment

- **Active env:** `orthophyl5`
- **Spec file:** `orthophyl_env.2.2.1.yml`
- **Activation:** `conda activate orthophyl5`

### Required Tools on PATH

The pipeline requires these external tools:
- `iqtree` — Tree inference
- `orthofinder` — Ortholog detection
- `mafft` — Multiple sequence alignment
- `hmmer` (hmmsearch, hmmbuild) — HMM-based gene finding (ReLeaf)
- `trimal` — Alignment trimming
- `mash` — Genome distance estimation (for routing)
- `fasttree` — Fast tree inference (alternative to IQ-TREE)
- `prokka` / `prodigal` — Genome annotation

### Test Suite

- **Config:** `pytest.ini`
- **Test deps:** `tests/requirements-test.txt` (pytest, pytest-cov, pytest-mock, ete3)
- **Unit tests:** `pytest tests/unit/` (fast, no external tools)
- **Integration tests:** `pytest tests/integration/ -m integration` (slow, ~26 min, requires full env)
  - Opt-in via `-m integration` or `ORTHOPHYL_RUN_INTEGRATION=1`
  - Runs real OrthoPhyl + ReLeaf on fasttest dataset

---

## 3. This Session's Work: Integration Test Fixes ✅ DONE

### Task Summary

**Goal:** Fix 4 failing integration tests in `tests/integration/test_orthophyl_fasttest.py` and `tests/integration/test_releaf_fasttest.py`

**Status:** ✅ **COMPLETE** — All tests now pass (15 passed, 1 skipped, 26m runtime)

**Verification:** See `integration_test_7_28.out` for full test output

### Root Cause

Tests had **hardcoded taxa counts** that didn't match the actual fasttest data:

| Test Expected | Actual Reality | Issue |
|---------------|----------------|-------|
| 8 OrthoPhyl taxa | 9 inputs | Wrong count |
| 2 ReLeaf additions | 5 additions | Wrong count |
| 6 original taxa | 9 original taxa | Wrong count |
| Raw filenames | Cleaned names | Char-stripping not applied |

**Character stripping issue:** OrthoPhyl's `CLEAN_N_COPY_GENOMES` function (in `script_lib/functions.sh`, lines ~144-147) strips these characters from taxon names:
- `(` and `)` — e.g., `(PE)` → removed
- `:` — e.g., `NC_026568.1:Cattleya` → `NC_026568.1_Cattleya`
- `_=_` — replaced with `_`

Tests were comparing raw filenames to cleaned tree labels, causing mismatches.

### Files Changed (5 total)

#### 1. `tests/integration/helpers/tree_compare.py`
**Added:** `clean_taxon_name()` function
- Replicates OrthoPhyl's char-stripping logic
- Ensures test expectations match pipeline's actual name normalization
- Used by both tests and fixtures

#### 2. `tests/integration/conftest.py`
**Added:** 3 new fixtures + helper function
- `_collect_basenames()` — Helper to scan input directories and apply `clean_taxon_name()`
- `expected_orthophyl_taxa` — Derives from `TESTER/genomes_fasttest/` + `TESTER/annots_prots_fasttest/` (9 taxa)
- `expected_releaf_added_taxa` — Derives from `TESTER/genomes_fasttest_addasm/` + `TESTER/annots_prots_fasttest_addasm/` (5 taxa)
- `expected_releaf_total_taxa` — Union of above (14 taxa)

**Why fixtures?** Tests now self-update if test data changes (no more hardcoded counts).

#### 3. `tests/integration/test_orthophyl_fasttest.py`
**Fixed 3 tests:**
- `test_all_input_taxa_in_tree` — Changed hardcoded 8→9 taxa expectation
- `test_source_inputs_match_tree` — **NEW TEST** using `expected_orthophyl_taxa` fixture (validates against source directories, catches copy/clean bugs)
- `test_genome_list_matches_input` — Now applies `clean_taxon_name()` to both sides of comparison (handles parentheses)

#### 4. `tests/integration/test_releaf_fasttest.py`
**Fixed 4 tests:**
- `test_added_taxa_in_releaf_tree` — Now uses `expected_releaf_added_taxa` fixture (5 taxa, not 2)
- `test_original_taxa_retained` — Now uses `expected_orthophyl_taxa` fixture (9 taxa, not 6)
- `test_releaf_tree_has_8_taxa` → **RENAMED** to `test_releaf_tree_has_expected_taxa` — Uses `expected_releaf_total_taxa` fixture (14 taxa, not 8)
- `test_releaf_preserves_orthophyl_outputs` — Fixed to expect 9 taxa (not 6)

**Added import:** `clean_taxon_name` to imports

### Design Decision: Option C (Best of Both Worlds)

We kept **both** the original output-file-reading test (`test_all_input_taxa_in_tree`) **and** added a new source-directory-reading test (`test_source_inputs_match_tree`). Why?

- **Output-file test** catches bugs in OrthoPhyl's internal bookkeeping (genome_list, pre_annotated_list)
- **Source-dir test** catches bugs in the copy/clean phase (inputs dropped before bookkeeping)

They validate different failure modes, so both are valuable.

### Key Domain Facts (fasttest dataset)

- **Dataset:** Orchid chloroplast genomes (small, fast for testing)
- **OrthoPhyl inputs:** 9 genomes
  - Source: `TESTER/genomes_fasttest/` (6 `.fasta` files) + `TESTER/annots_prots_fasttest/` (3 `.faa` files)
- **ReLeaf additions:** 5 genomes
  - Source: `TESTER/genomes_fasttest_addasm/` (3 `.fasta`) + `TESTER/annots_prots_fasttest_addasm/` (2 `.faa`)
- **Final ReLeaf tree:** 14 taxa (9 original + 5 added)
- **OrthoPhyl run uses `-n 5`** → Triggers MASH shortlist path (not all-vs-all)

### Character Stripping Rules (Critical!)

When comparing taxon names between tests and pipeline outputs, always apply:

```python
def clean_taxon_name(name: str) -> str:
    """Replicate OrthoPhyl's CLEAN_N_COPY_GENOMES char-stripping."""
    name = name.replace('(', '').replace(')', '')  # Remove parens
    name = name.replace(':', '_')                   # Colon → underscore
    name = name.replace('_=_', '_')                 # Normalize underscores
    return name
```

**Location in pipeline:** `script_lib/functions.sh`, function `CLEAN_N_COPY_GENOMES` (lines ~144-147)

---

## 4. Test Suite Architecture

### Directory Structure

```
tests/
├── conftest.py                          # Root fixtures (project_root, etc.)
├── README.md                            # Test suite overview
├── requirements-test.txt                # pytest, ete3, etc.
├── integration/
│   ├── __init__.py
│   ├── conftest.py                      # Integration fixtures (orthophyl_run, expected_*)
│   ├── test_orthophyl_fasttest.py       # 9 tests (8 active, 1 skipped)
│   ├── test_releaf_fasttest.py          # 7 tests (all active)
│   └── helpers/
│       ├── __init__.py
│       └── tree_compare.py              # Tree utilities (load_tree, clean_taxon_name, etc.)
└── unit/
    ├── test_db_creator.py               # Database creation logic
    ├── test_gtdb_taxonomy.py            # GTDB taxonomy parsing
    ├── test_router.py                   # Assembly router (MASH-based)
    ├── test_taxon_assembly_gatherer.py  # Taxon-mode NCBI fetching
    ├── test_wrapper_batch.py            # Wrapper batch-mode logic
    └── test_wrapper_taxon.py            # Wrapper taxon-mode logic
```

### Key Fixtures

#### Session-Scoped (Expensive, Run Once)

- **`orthophyl_run`** (`tests/integration/conftest.py`)
  - Runs full OrthoPhyl.sh on fasttest dataset (9 inputs, `-n 5`)
  - Takes ~15-20 minutes
  - Returns `Path` to output directory
  - Used by all OrthoPhyl + ReLeaf tests

#### Derived Fixtures (Fast, Computed from Disk)

- **`expected_orthophyl_taxa`** — Set of cleaned taxon names from OrthoPhyl input dirs (9 taxa)
- **`expected_releaf_added_taxa`** — Set of cleaned taxon names from ReLeaf addasm dirs (5 taxa)
- **`expected_releaf_total_taxa`** — Union of above (14 taxa)

These fixtures make tests **self-updating**: if test data changes, counts adjust automatically.

### How to Run Tests

```bash
# Unit tests (fast, no external tools required)
pytest tests/unit/ -v

# Integration tests (slow, requires orthophyl5 env + tools on PATH)
conda activate orthophyl5
pytest tests/integration/ -v -m integration

# Run specific test
pytest tests/integration/test_orthophyl_fasttest.py::TestOrthoPhylFasttest::test_all_input_taxa_in_tree -v

# With coverage
pytest tests/unit/ --cov=orthophyl_pipeline_wrapper --cov-report=html
```

### Runtime Expectations

- **Unit tests:** < 1 minute total
- **Integration tests:** ~26 minutes (OrthoPhyl run ~15-20m, ReLeaf run ~5-10m, validation ~1m)

---

## 5. Key File/Path Reference

### Main Pipeline Scripts

| File | Purpose |
|------|---------|
| `OrthoPhyl.sh` | Main pipeline entry point (de novo phylogeny) |
| `ReLeaf.sh` | Add genomes to existing run (HMM-based) |
| `orthophyl_pipeline_wrapper.v2.py` | Database manager + router (taxon-mode + batch-mode) |
| `script_lib/functions.sh` | Core bash functions (CLEAN_N_COPY_GENOMES, etc.) |
| `script_lib/functions_addem.sh` | ReLeaf-specific bash functions |
| `script_lib/arg_parse.sh` | OrthoPhyl.sh argument parser |
| `script_lib/arg_parse_addem.sh` | ReLeaf.sh argument parser |

### Assembly Router / Database

| File | Purpose |
|------|---------|
| `assembly_router/assembly_router.py` | Original router (MASH-based distance) |
| `assembly_router/assembly_router_hierarchical.py` | Hierarchical DB router (v2) |
| `assembly_router/create_hierarchical_database_v3.py` | Latest DB creator |
| `assembly_router/add_releaf_version.py` | ReLeaf versioner (timestamps new runs) |
| `RELEAF_DATABASE_DESIGN.md` | Design doc for hierarchical DB |

### Utilities

| File | Purpose |
|------|---------|
| `utils/taxon_assembly_gatherer.py` | Fetch assemblies from NCBI by taxon name |
| `utils/ncbi_assembly_stats.py` | Parse NCBI assembly stats |
| `utils/gather_filter_asms.sh` | Bash wrapper for assembly gathering |

### Test Data

| Path | Contents |
|------|----------|
| `TESTER/genomes_fasttest/` | 6 orchid chloroplast genomes (`.fasta`) |
| `TESTER/annots_prots_fasttest/` | 3 pre-annotated protein sets (`.faa`) |
| `TESTER/genomes_fasttest_addasm/` | 3 genomes for ReLeaf addition |
| `TESTER/annots_prots_fasttest_addasm/` | 2 pre-annotated for ReLeaf addition |
| `TESTER/REFERENCE_TESTER_TREES/` | Reference trees for topology comparison |

### Critical Code Locations

- **Char-stripping:** `script_lib/functions.sh`, function `CLEAN_N_COPY_GENOMES` (lines ~144-147)
- **ReLeaf CLI parser:** `script_lib/arg_parse_addem.sh` (lines 31-229)
- **ReLeaf output location:** `ReLeaf.sh` line 113 (`addasm_dir=$store/ReLeaf_dir`)

---

## 6. In-Flight / Related Work

### 6.1 Wrapper v2 Fixes (CRITICAL BUGS) — In Progress

**Source:** `work_plan_7_23_26.md` (detailed 383-line plan)

**Status:** Pre-implementation analysis complete, implementation not started

#### Three Critical Issues

1. **False success messages** — Pipeline reports "All assemblies placed!" even when ReLeaf/OrthoPhyl fail
2. **Results directory confusion** — Need to default to database-dir and name results by taxon
3. **ReLeaf invocation failure** — Wrong CLI flags cause immediate exit, empty results

#### Issue 3: ReLeaf Invocation Failure (PRIMARY BUG)

**Problem:** Wrapper calls ReLeaf.sh with incorrect flags:

```python
# WRONG (current code)
cmd = [
    str(self.releaf_script),
    '--store', str(database_dir / 'orthophyl_run'),      # ❌ Should be -s or --storage_dir
    '--input_genomes', str(input_genomes),               # ❌ Should be -g or --genome_dir
    '--tree_method', tree_method,                        # ❌ Should be -p or --phylo_tool
    '--TREE_DATA', tree_data                             # ❌ Should be -o or --omics
]
```

**Correct flags** (from `script_lib/arg_parse_addem.sh`):
- `-s|--storage_dir` (not `--store`)
- `-g|--genome_dir` (not `--input_genomes`)
- `-p|--phylo_tool` (not `--tree_method`)
- `-o|--omics` (not `--TREE_DATA`, expects `CDS|PROT|BOTH`)

**Result:** ReLeaf.sh hits `*) USAGE; exit 1` and exits immediately without creating any output.

**Additional issue:** ReLeaf hardcodes output to `$store/ReLeaf_dir` (line 113 of ReLeaf.sh), ignoring the wrapper's `cwd=output_dir`. Aggregation looks in the wrong place.

#### Implementation Plan (from work_plan)

**Phase 1:** Fix ReLeaf invocation (Issue 3) — CRITICAL
- Fix command flags in `_run_releaf()` (lines 450-506)
- Fix output path expectations (look in `database_dir / 'orthophyl_run' / 'ReLeaf_dir'`)
- Add post-run verification (check return code + expected files exist)
- Handle stale ReLeaf_dir

**Phase 2:** Implement real error handling (Issue 1)
- Add `pipeline_status` tracking dict to `__init__()`
- Track successes/failures in `_run_releaf()` and `_run_orthophyl()`
- Rewrite `_generate_summary_report()` to show actual counts
- Return non-zero exit code if failures occurred

**Phase 3:** Default to database dir (Issue 2)
- Make `--output-dir` optional (default to `<database-dir>/.pipeline_runs/<timestamp>`)
- Results named by taxon, placed under database-dir

#### Testing Checklist

**Pre-implementation (all ✅ done):**
- [x] Confirmed ReLeaf.sh arg parser only accepts `-s`, `-g`, `-t`, `-p`, `-o`
- [x] Confirmed ReLeaf writes to `$store/ReLeaf_dir`
- [x] Confirmed versioner expects `new_prot_alignments.trm.nm`, `new_CDS_alignments.trm.nm`, `new_trees/`
- [x] Confirmed summary always prints success message

**Post-implementation (all ❌ not started):**
- [ ] Dry-run mode works without errors
- [ ] Single ReLeaf-matched taxon runs successfully
- [ ] Single OrthoPhyl novel taxon runs successfully
- [ ] Intentional failure caught and reported
- [ ] Mixed batch (ReLeaf + OrthoPhyl) works
- [ ] Resume from checkpoint works
- [ ] Stale ReLeaf_dir detected and handled

### 6.2 Taxon Mode — Implemented, Testing in Progress

**Source:** `TAXON_MODE_GUIDE.md`, `Taxon_mode_implementation_notes.md`

**What it is:** Wrapper mode that fetches assemblies from NCBI by taxon name, then routes them through the pipeline.

**Key files:**
- `utils/taxon_assembly_gatherer.py` — NCBI Datasets API wrapper
- `tests/unit/test_taxon_assembly_gatherer.py` — Unit tests
- `tests/unit/test_wrapper_taxon.py` — Wrapper taxon-mode tests

**Status:** Implementation complete, unit tests passing, integration testing ongoing.

### 6.3 Assembly Router / Hierarchical Database — v3 in Development

**Source:** `RELEAF_DATABASE_DESIGN.md`, `assembly_router/` directory

**What it is:** Hierarchical database system for organizing OrthoPhyl runs by taxonomy, enabling efficient ReLeaf routing.

**Key files:**
- `assembly_router/create_hierarchical_database_v3.py` — Latest DB creator
- `assembly_router/add_releaf_version.py` — ReLeaf versioner
- `assembly_router/assembly_router_hierarchical.py` — Router using hierarchical DB

**Status:** v3 in active development (v2 exists but being superseded).

**Design:** Each database entry is a completed OrthoPhyl run. ReLeaf adds genomes to the closest-matching database (by MASH distance). Versioner creates timestamped snapshots of ReLeaf outputs.

---

## 7. Known TODOs / Open Threads

### Integration Tests

- [ ] **Generate reference tree** for `test_topology_matches_reference` (currently skipped)
  - Need: `TESTER/REFERENCE_TESTER_TREES/fasttest_n5.SCO_strict.CDS.tree`
  - How: Run OrthoPhyl.sh on fasttest with `-n 5`, save tree as reference
  - Why: Enables topology regression testing (RF distance = 0)

### Wrapper v2

- [ ] **Implement Phase 1** (fix ReLeaf CLI flags) — CRITICAL
- [ ] **Implement Phase 2** (real error handling)
- [ ] **Implement Phase 3** (default to database-dir)
- [ ] **Complete post-implementation testing checklist** (7 items)

### Test Coverage

- [ ] **Add integration test for taxon mode** (currently only unit tests)
- [ ] **Add integration test for wrapper batch mode** (end-to-end with real router)
- [ ] **Consider adding unit tests for** `script_lib/functions.sh` bash functions (if feasible)

### Housekeeping

- [ ] **Clean up duplicate router scripts** (`assembly_router_multi.cmd_out*.py` — 3 versions)
- [ ] **Archive old test outputs** (`TESTER/FULLTEST_OUT.*` — 20+ directories)
- [ ] **Consolidate ncbi_stats directories** (`ncbi_stats/`, `ncbi_stats2/`, ..., `ncbi_stats5/`)

---

## 8. Gotchas / Non-Obvious Facts

### Integration Tests Take ~26 Minutes

The `orthophyl_run` fixture runs the full OrthoPhyl.sh pipeline on 9 genomes. This is intentional (real-world validation), but means:
- Don't run integration tests in a tight dev loop
- Use unit tests for rapid iteration
- Integration tests are opt-in (`-m integration`)

### Character Stripping Happens in the Pipeline

Any test comparing taxon names (from filenames vs. tree labels) **must** apply `clean_taxon_name()`. The pipeline strips `()`, `:`, `_=_` during the copy/clean phase.

**Where:** `script_lib/functions.sh`, function `CLEAN_N_COPY_GENOMES`

### ReLeaf Overwrites Output (No Built-In Versioning)

ReLeaf.sh writes to `$store/ReLeaf_dir` and **overwrites** it on each run. The versioner (`add_releaf_version.py`) creates timestamped snapshots in the database, but the raw ReLeaf output is ephemeral.

**Implication:** If you run ReLeaf twice on the same database without versioning, the first run's output is lost.

### fasttest Uses `-n 5` (MASH Shortlist Path)

The integration tests run OrthoPhyl.sh with `-n 5`, which triggers the MASH shortlist path (not all-vs-all OrthoFinder). This is faster but produces different intermediate files than the default path.

**Why it matters:** If you're debugging test failures, check whether the issue is specific to the `-n 5` path.

### Wrapper CLI Flag Mismatch is Silent

The wrapper's incorrect ReLeaf flags (`--store`, `--input_genomes`, etc.) cause ReLeaf.sh to print usage and exit with code 1. But the wrapper doesn't check the return code, so it appears to succeed. This is **Issue 3** in the wrapper v2 work plan.

### Test Fixtures are Session-Scoped

The `orthophyl_run` fixture is `scope="session"`, meaning it runs **once** for all tests in a session. If you modify the fasttest input data, you must restart pytest (or use `--setup-show` to see fixture reuse).

---

## 9. Quick Start for a New AI Instance

If you're a fresh Opus/Sonnet instance picking up this project:

1. **Read this doc** (you're doing it!)
2. **Skim the key docs** (README.v2.md, PIPELINE_SUMMARY.md, work_plan_7_23_26.md)
3. **Understand the test suite** (section 4 above)
4. **Check in-flight work** (section 6 above)
5. **Run the tests** to verify your environment:
   ```bash
   conda activate orthophyl5
   pytest tests/unit/ -v                          # Should pass quickly
   pytest tests/integration/ -v -m integration    # Takes ~26 min
   ```
6. **Pick up where we left off:**
   - If working on wrapper v2: Start with Phase 1 (ReLeaf CLI flags) in `work_plan_7_23_26.md`
   - If working on tests: Check the TODOs in section 7
   - If working on taxon mode: See `TAXON_MODE_GUIDE.md`

---

## 10. Session Metadata

- **Date:** 2026-07-28
- **Primary work:** Integration test fixes (4 failing tests → all passing)
- **Time spent:** ~2 hours (diagnosis, implementation, verification)
- **Test results:** 15 passed, 1 skipped, 0 failed (26m runtime)
- **Verification file:** `integration_test_7_28.out`
- **Git commit at session end:** `48a25e19f322e925b3f552bab198ed1f5132ed0e`

---

**End of handoff document. Good luck! 🚀**
