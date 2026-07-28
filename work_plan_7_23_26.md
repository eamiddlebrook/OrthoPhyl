# Work Plan: Pipeline Wrapper Fixes (7/23/2026)

## Overview
This plan addresses three critical issues discovered during testing of `orthophyl_pipeline_wrapper.v2.py`:

1. **False success messages** - Pipeline reports success even when ReLeaf/OrthoPhyl fail
2. **Results directory confusion** - Need to default to database-dir and name results by taxon
3. **ReLeaf invocation failure** - Wrong CLI flags cause immediate exit, empty results

## Root Cause Analysis

### Issue 3: ReLeaf Invocation Failure (PRIMARY BUG)
**Location:** `orthophyl_pipeline_wrapper.v2.py`, lines 450-506 (`_run_releaf()`)

**Problem:** The wrapper calls ReLeaf.sh with incorrect flags:
```python
cmd = [
    str(self.releaf_script),
    '--store', str(database_dir / 'orthophyl_run'),
    '--input_genomes', str(input_genomes),
    '-t', str(self.threads),
    '--tree_method', tree_method,
    '--TREE_DATA', tree_data
]
```

**But** `script_lib/arg_parse_addem.sh` only accepts:
- `-s|--storage_dir` (not `--store`)
- `-g|--genome_dir` (not `--input_genomes`)
- `-p|--phylo_tool` (not `--tree_method`)
- `-o|--omics` (not `--TREE_DATA`, and expects `CDS|PROT|BOTH`)

**Result:** ReLeaf.sh hits the `*) USAGE; exit 1` case and exits immediately without creating any output.

**Additional issue:** ReLeaf hardcodes output to `addasm_dir=$store/ReLeaf_dir` (line 113 of ReLeaf.sh), ignoring the wrapper's `cwd=output_dir`. So even if it ran:
- Aggregation looks for `01_releaf_only/<db>/ReLeaf_results/phylogeny_with_new_genomes.nwk` → never exists
- Versioner is passed the empty `output_dir` instead of `$store/ReLeaf_dir` → validation fails

### Issue 1: False Success Messages
**Location:** `orthophyl_pipeline_wrapper.v2.py`, lines 902-948 (`_generate_summary_report()`)

**Problem:** Line 945 unconditionally writes:
```python
f.write("All assemblies have been placed in phylogenetic trees!\n")
```

No verification that:
- ReLeaf/OrthoPhyl actually succeeded (return code could be non-zero but caught)
- Expected output files exist
- Per-database runs completed successfully

**Result:** User sees "All assemblies have been placed" even when everything failed.

### Issue 2: Results Directory Design
**Current behavior:** 
- `--output-dir` is required
- Results scattered: routing in `00_routing/`, ReLeaf in `01_releaf_only/`, OrthoPhyl in `02_orthophyl_novel/`, aggregated in `03_results/`
- ReLeaf actually writes to `$store/ReLeaf_dir` (inside database), not the wrapper's output dir

**Desired behavior:**
- Wrapper is a **database manager**: `--database-dir` is the primary anchor
- `--output-dir` becomes optional (defaults to database-dir)
- Results named by taxon and placed under database-dir:
  - **ReLeaf route:** `<database-dir>/<matched_db>/releaf_<taxon>/`
  - **OrthoPhyl route:** `<database-dir>/<taxon>/` (becomes new DB entry)
- Working dirs (routing, logs, checkpoints) under `<database-dir>/.pipeline/` or `--output-dir` if provided

## Implementation Plan

### Phase 1: Fix ReLeaf Invocation (Issue 3) - CRITICAL
**File:** `orthophyl_pipeline_wrapper.v2.py`

**Changes to `_run_releaf()` (lines 450-506):**

1. **Fix command flags:**
```python
cmd = [
    str(self.releaf_script),
    '-s', str(database_dir / 'orthophyl_run'),  # --storage_dir
    '-g', str(input_genomes),                    # --genome_dir
    '-t', str(self.threads),                     # --threads
    '-p', tree_method,                           # --phylo_tool
    '-o', tree_data                              # --omics (CDS|PROT|BOTH)
]
```

2. **Fix output path expectations:**
   - ReLeaf writes to `$store/ReLeaf_dir`
   - Update aggregation to look there: `database_dir / 'orthophyl_run' / 'ReLeaf_dir'`
   - Pass correct path to versioner: `releaf_output_dir = database_dir / 'orthophyl_run' / 'ReLeaf_dir'`

3. **Add post-run verification:**
```python
# After subprocess.run()
if result.returncode != 0:
    raise RuntimeError(f"ReLeaf failed for {database_name}. Check log: {log_file}")

# Verify expected outputs exist
releaf_output = database_dir / 'orthophyl_run' / 'ReLeaf_dir'
expected_files = [
    releaf_output / 'new_prot_alignments.trm.nm',
    releaf_output / 'new_CDS_alignments.trm.nm',
    releaf_output / 'new_trees'
]
missing = [f for f in expected_files if not f.exists()]
if missing:
    raise RuntimeError(
        f"ReLeaf completed but missing expected outputs:\n" +
        "\n".join(f"  - {f}" for f in missing) +
        f"\nCheck log: {log_file}"
    )
```

4. **Handle stale ReLeaf_dir:**
```python
# Before running ReLeaf
releaf_dir = database_dir / 'orthophyl_run' / 'ReLeaf_dir'
if releaf_dir.exists() and not self.resume:
    logger.warning(f"  ⚠ Removing stale ReLeaf output: {releaf_dir}")
    shutil.rmtree(releaf_dir)
```

**Changes to `_create_releaf_version()` (lines 514-558):**
```python
def _create_releaf_version(self, database_name: str, database_dir: Path):
    """Create a new database version from ReLeaf output."""
    logger.info(f"\n  Creating new database version from ReLeaf output...")
    
    # ReLeaf writes to $store/ReLeaf_dir
    releaf_output_dir = database_dir / 'orthophyl_run' / 'ReLeaf_dir'
    
    if not releaf_output_dir.exists():
        logger.warning(f"  ⚠ ReLeaf output not found: {releaf_output_dir}")
        return
    
    # ... rest of function uses releaf_output_dir
```

**Changes to `_phase_aggregation()` (lines 860-900):**
```python
# Collect ReLeaf trees - look in database's ReLeaf_dir
releaf_trees = []
for db_dir in self.database_dir.glob("*_db"):
    releaf_output = db_dir / 'orthophyl_run' / 'ReLeaf_dir'
    tree_file = releaf_output / "new_trees" / "SCO_strict.CDS.iqtree.treefile.addasm"
    if tree_file.exists():
        dst = trees_dir / "releaf" / f"{db_dir.name}_phylogeny.nwk"
        shutil.copy(tree_file, dst)
        releaf_trees.append(db_dir.name)
        logger.info(f"  ✓ Collected ReLeaf tree: {db_dir.name}")
```

### Phase 2: Implement Real Error Handling (Issue 1)
**File:** `orthophyl_pipeline_wrapper.v2.py`

**Changes to `PipelineWrapper.__init__()` (lines 46-110):**
```python
# Add tracking for successes/failures
self.pipeline_status = {
    'start_time': datetime.now().isoformat(),
    'phases': {},
    'summary': {},
    'successes': [],
    'failures': []
}
```

**Changes to `_run_releaf()` (lines 450-506):**
```python
# After verification
self.pipeline_status['successes'].append({
    'type': 'releaf',
    'database': database_name,
    'assemblies': len(assemblies)  # pass this in
})

# In exception handler
except Exception as e:
    self.pipeline_status['failures'].append({
        'type': 'releaf',
        'database': database_name,
        'error': str(e),
        'log': str(log_file)
    })
    raise  # or log and continue depending on desired behavior
```

**Changes to `_run_orthophyl()` (lines 733-792):**
```python
# After verification
self.pipeline_status['successes'].append({
    'type': 'orthophyl',
    'taxon': taxon_name,
    'assemblies': len(assemblies)
})

# In exception handler
except Exception as e:
    self.pipeline_status['failures'].append({
        'type': 'orthophyl',
        'taxon': taxon_name,
        'error': str(e),
        'log': str(log_file)
    })
    raise
```

**Rewrite `_generate_summary_report()` (lines 902-948):**
```python
def _generate_summary_report(self, releaf_trees: List[str], orthophyl_trees: List[str]):
    """Generate human-readable summary report with actual success/failure counts."""
    report_file = self.results_dir / "pipeline_summary.txt"
    
    successes = self.pipeline_status.get('successes', [])
    failures = self.pipeline_status.get('failures', [])
    
    with open(report_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("ORTHOPHYL PIPELINE - SUMMARY REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"Input File: {self.input_file}\n")
        f.write(f"Database Directory: {self.database_dir}\n")
        f.write(f"Output Directory: {self.output_dir}\n\n")
        
        f.write("RESULTS:\n")
        f.write("-" * 70 + "\n\n")
        
        # Success counts
        f.write(f"✓ Successful placements: {len(successes)}\n")
        for s in successes:
            if s['type'] == 'releaf':
                f.write(f"  - ReLeaf: {s['database']} ({s.get('assemblies', '?')} assemblies)\n")
            else:
                f.write(f"  - OrthoPhyl: {s['taxon']} ({s.get('assemblies', '?')} assemblies)\n")
        f.write("\n")
        
        # Failure counts
        if failures:
            f.write(f"✗ Failed placements: {len(failures)}\n")
            for fail in failures:
                if fail['type'] == 'releaf':
                    f.write(f"  - ReLeaf: {fail['database']}\n")
                else:
                    f.write(f"  - OrthoPhyl: {fail['taxon']}\n")
                f.write(f"    Error: {fail['error']}\n")
                f.write(f"    Log: {fail['log']}\n")
            f.write("\n")
        
        f.write("OUTPUT LOCATIONS:\n")
        f.write(f"  Trees: {self.results_dir}/trees/\n")
        f.write(f"  Logs: {self.logs_dir}/\n")
        f.write(f"  Routing decisions: {self.routing_dir}/\n")
        f.write("\n")
        
        # Only show success banner if no failures
        if not failures:
            f.write("=" * 70 + "\n")
            f.write("All assemblies have been placed in phylogenetic trees!\n")
            f.write("=" * 70 + "\n")
        else:
            f.write("=" * 70 + "\n")
            f.write(f"PIPELINE COMPLETED WITH {len(failures)} FAILURE(S)\n")
            f.write("=" * 70 + "\n")
    
    logger.info(f"  Summary report: {report_file}")
    
    # Return failure count for exit code
    return len(failures)
```

**Update `_run_batch_mode()` (lines 164-185):**
```python
def _run_batch_mode(self) -> int:
    """Run original batch mode workflow."""
    # ... existing code ...
    
    # Phase 4: Results aggregation
    self._phase_aggregation()
    
    # Check for failures
    failures = len(self.pipeline_status.get('failures', []))
    
    logger.info("=" * 70)
    if failures == 0:
        logger.info("PIPELINE COMPLETE!")
    else:
        logger.error(f"PIPELINE COMPLETED WITH {failures} FAILURE(S)")
    logger.info("=" * 70)
    
    self._save_final_status()
    return 1 if failures > 0 else 0
```

### Phase 3: Default to Database Dir (Issue 2)
**File:** `orthophyl_pipeline_wrapper.v2.py`

**Changes to `main()` argparse (lines 1337-1492):**
```python
parser.add_argument(
    '--output-dir',
    help='Output directory for results (default: uses --database-dir)'
)
# Remove required=True
```

**Changes to `PipelineWrapper.__init__()` (lines 46-110):**
```python
def __init__(
    self,
    input_file: Optional[Path] = None,
    database_dir: Path = None,
    output_dir: Optional[Path] = None,  # Now optional
    # ... rest of params
):
    self.input_file = Path(input_file) if input_file else None
    self.database_dir = Path(database_dir) if database_dir else None
    
    # Default output_dir to database_dir if not provided
    if output_dir:
        self.output_dir = Path(output_dir)
    else:
        self.output_dir = self.database_dir / '.pipeline_runs' / datetime.now().strftime('%Y%m%d_%H%M%S')
        logger.info(f"No --output-dir provided, using: {self.output_dir}")
    
    # ... rest of init
```

**Changes to result path computation:**
- For ReLeaf: results already in `<database-dir>/<matched_db>/orthophyl_run/ReLeaf_dir/`
- For OrthoPhyl: create under `<database-dir>/<taxon>/` and register as new DB
- Working dirs under `self.output_dir` (which defaults to `.pipeline_runs/`)

## Testing Checklist

### Pre-implementation verification
- [x] Confirmed ReLeaf.sh arg parser only accepts `-s`, `-g`, `-t`, `-p`, `-o`
- [x] Confirmed ReLeaf writes to `$store/ReLeaf_dir`
- [x] Confirmed versioner expects `new_prot_alignments.trm.nm`, `new_CDS_alignments.trm.nm`, `new_trees/`
- [x] Confirmed summary always prints success message

### Post-implementation testing
- [ ] Dry-run mode works without errors
- [ ] Single ReLeaf-matched taxon:
  - [ ] ReLeaf runs successfully (check log for correct flags)
  - [ ] Results appear in `<database-dir>/<matched_db>/orthophyl_run/ReLeaf_dir/`
  - [ ] Versioner creates new DB version successfully
  - [ ] Summary reports 1 success, 0 failures
- [ ] Single OrthoPhyl novel taxon:
  - [ ] OrthoPhyl runs successfully
  - [ ] Results appear in `<database-dir>/<taxon>/`
  - [ ] New database entry created
  - [ ] Summary reports 1 success, 0 failures
- [ ] Intentional failure (bad assembly path):
  - [ ] Pipeline catches failure
  - [ ] Summary reports 0 successes, 1 failure with log path
  - [ ] Exit code is non-zero
- [ ] Mixed batch (some ReLeaf, some OrthoPhyl):
  - [ ] All routes execute correctly
  - [ ] Summary shows correct counts
- [ ] Resume from checkpoint works
- [ ] Stale ReLeaf_dir is detected and handled

## Files Modified
1. `orthophyl_pipeline_wrapper.v2.py` - main fixes
2. `work_plan_7_23_26.md` - this document

## References
- ReLeaf arg parser: `script_lib/arg_parse_addem.sh` lines 31-229
- ReLeaf output location: `ReLeaf.sh` line 113 (`addasm_dir=$store/ReLeaf_dir`)
- Versioner validation: `assembly_router/add_releaf_version.py` lines 66-104
- Wrapper ReLeaf invocation: `orthophyl_pipeline_wrapper.v2.py` lines 450-506
- Wrapper summary: `orthophyl_pipeline_wrapper.v2.py` lines 902-948

## Implementation Order
1. Phase 1 (Issue 3) - Unblocks ReLeaf execution
2. Phase 2 (Issue 1) - Provides honest feedback
3. Phase 3 (Issue 2) - Cleans up directory structure

## Notes
- ReLeaf.sh overwrites `$store/ReLeaf_dir` on each run (no built-in versioning)
- The versioner (`add_releaf_version.py`) creates timestamped versions in the database
- Consider adding a `--force` flag to allow overwriting stale ReLeaf_dir without prompting
