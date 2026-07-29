# Genome Download Harmonization — Summary

**Date:** 2026-07-29  
**Scope:** Taxon-mode create path in `orthophyl_pipeline_wrapper.v2.py`

## Problem

The wrapper had **two divergent genome-gathering code paths**:

1. **Batch mode / OrthoPhyl route** → `_download_genomes()` → `utils/gather_filter_asms.sh`
   - Uses NCBI `datasets` CLI (multithreaded, robust gz compression, retry/backoff)
   - **Runs CheckM/BBMap QC filtering** (completeness ≥95%, contamination ≤1%, N50 stats)
   - Filtered genomes land in `genomes_to_keep/`

2. **Taxon mode (create/update)** → `TaxonAssemblyGatherer.download_assemblies()`
   - Uses Python `urllib` FTP downloads (single-threaded, prone to corruption)
   - **No QC filtering** — incomplete assemblies leak into trees

This meant taxon mode was both **slower** and **lower quality** than the main pipeline.

## Solution (Approach A — Create Mode Only)

Per user decision, we harmonized **create mode only** (update mode deferred):

### Changes to `orthophyl_pipeline_wrapper.v2.py`

#### 1. `_run_taxon_create_mode()` (lines ~1238–1275)
**Before:**
```python
gatherer.download_assemblies(assemblies, download_dir)
self._run_orthophyl(input_dir=download_dir, ...)
```

**After:**
```python
# Hard error if gather_script not provided
if not self.gather_script or not self.gather_script.exists():
    logger.error("ERROR: Genome download script is required for taxon create mode")
    logger.error("  Please provide --gather-script utils/gather_filter_asms.sh")
    return 1

# Download via gather_filter_asms.sh (with QC filtering)
if not self.skip_download and not self._verify_download_complete(download_dir, self.taxon):
    self._download_genomes(self.taxon, download_dir)
    self._write_checkpoint(f"download_{self.taxon}")

# Verify filtered genomes exist
genomes_to_keep = download_dir / "genomes_to_keep"
if not genomes_to_keep.exists() or not list(genomes_to_keep.glob("*.fna")):
    logger.error("ERROR: No genomes passed QC filtering")
    return 1

# Use filtered genomes for OrthoPhyl
self._run_orthophyl(input_dir=genomes_to_keep, ...)  # Not raw download_dir
```

**Benefits:**
- NCBI `datasets` multithreaded download (faster, more robust)
- Automatic retry/backoff on network errors
- **QC filtering** — only assemblies with completeness ≥95%, contamination ≤1% proceed to tree
- Respects `--low-ram` (`--reduced_tree`) and `--use-bbmap` flags
- Checkpoint/resume support (skip re-downloading on re-runs)

#### 2. `_run_taxon_update_mode()` (lines ~1327–1330)
**Added TODO comment:**
```python
# TODO: New assemblies downloaded here via TaxonAssemblyGatherer are NOT yet
# QC-filtered (completeness/contamination/N50). Harmonize with gather_filter_asms.sh
# filtering in a future update. See create-mode for the filtered path.
gatherer.download_assemblies(new_assemblies, download_dir)
```

**No behavioral change** — update mode still uses `TaxonAssemblyGatherer.download_assemblies()` (unfiltered). This is a **known limitation** to be addressed in a future update.

### Files Modified
- `orthophyl_pipeline_wrapper.v2.py` — create-mode download path + TODO note

### Files NOT Modified
- `utils/gather_filter_asms.sh` — no changes (no accession-list mode needed for Approach A)
- `utils/taxon_assembly_gatherer.py` — no deprecation note (still used by update mode)

## Testing Plan (Deferred)

Update `tests/unit/test_wrapper_taxon.py`:
- **Create mode:** assert it calls `_download_genomes` (not `gatherer.download_assemblies`), assert OrthoPhyl receives `genomes_to_keep/`, assert hard-error when `gather_script` missing
- **Update mode:** keep existing expectation that it calls `gatherer.download_assemblies` (unchanged behavior)
- Keep contract test requiring `query_ncbi` / `get_taxonomy_string` on the fake gatherer

## Usage

**Before (broken — no QC filter):**
```bash
python orthophyl_pipeline_wrapper.v2.py \
  --taxon "Methylorubrum" \
  --database-dir databases/ \
  --output-dir results/
# ❌ Downloads via urllib FTP, no filtering, slow
```

**After (fixed — QC filtered):**
```bash
python orthophyl_pipeline_wrapper.v2.py \
  --taxon "Methylorubrum" \
  --database-dir databases/ \
  --output-dir results/ \
  --gather-script utils/gather_filter_asms.sh \
  --threads 32
# ✅ Downloads via NCBI datasets, QC-filtered (completeness/contamination/N50), fast
```

## Future Work

1. **Update mode harmonization** — route `_run_taxon_update_mode()` through `gather_filter_asms.sh` as well (requires either downloading full taxon + subsetting, or adding `--accession-list` mode to the script).
2. **Test suite** — update `tests/unit/test_wrapper_taxon.py` per the testing plan above.
3. **Optional:** Add `--min-completeness` / `--max-contamination` flags to wrapper to expose the QC thresholds (currently hardcoded in `gather_filter_asms.sh` as 95% / 1%).

## Verification

To verify the fix works:
```bash
# Create a new taxon database (should now use gather_filter_asms.sh)
python orthophyl_pipeline_wrapper.v2.py \
  --taxon "Methylorubrum" \
  --database-dir test_dbs/ \
  --output-dir test_run/ \
  --gather-script utils/gather_filter_asms.sh \
  --threads 8

# Check logs — should see:
#   "Downloading and QC-filtering N assemblies..."
#   "Using: utils/gather_filter_asms.sh"
#   "X genomes passed QC filtering"

# Verify genomes are in genomes_to_keep/
ls test_run/downloaded_assemblies/genomes_to_keep/*.fna
```

---

**Summary:** Taxon create mode now uses the same robust, QC-filtered download path as batch mode. Update mode remains unfiltered (TODO for future work).
