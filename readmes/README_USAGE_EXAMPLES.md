# Usage Examples - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---


### Example 1: Basic Run with Existing Databases

```bash
python orthophyl_pipeline_wrapper.py \
    --input my_assemblies.tsv \
    --database-dir /data/taxonomy_databases/ \
    --output-dir /results/run_2025-01-15/ \
    --threads 32
```

**Scenario**: You have pre-built databases and want to place new assemblies

**What happens**:
- Routes assemblies to existing databases
- Runs ReLeaf for matches
- Skips OrthoPhyl (no novel taxa)

---

### Example 2: Complete Run with Genome Downloading

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    -v
```

**Scenario**: You have novel taxa and want automatic genome downloading

**What happens**:
- Routes assemblies
- Runs ReLeaf for matches
- Downloads related genomes for novel taxa
- Runs OrthoPhyl on expanded genome sets
- Creates new databases

---

### Example 3: Initial Setup (Create Databases from Scratch)

```bash
# Step 1: Create orthophyl_runs.tsv
cat > orthophyl_runs.tsv << EOF
Rhizobiaceae	/data/orthophyl_runs/rhizobiaceae	d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae
Enterobacteriaceae	/data/orthophyl_runs/enterobacteriaceae	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae
EOF

# Step 2: Run wrapper
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --orthophyl-runs orthophyl_runs.tsv \
    --threads 32
```

**Scenario**: First time setup, no databases exist yet

**What happens**:
- Creates initial databases from orthophyl_runs.tsv
- Then proceeds with normal routing

---

### Example 4: Resume After Interruption

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    --resume
```

**Scenario**: Pipeline was interrupted (power failure, timeout, etc.)

**What happens**:
- Checks checkpoint flags
- Skips completed phases
- Resumes from last incomplete step

---

### Example 5: Dry Run (Preview Mode)

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    --dry-run \
    -vv
```

**Scenario**: Want to see what would happen without actually running

**What happens**:
- Shows all commands that would be executed
- Creates directory structure
- No actual computation
- Useful for debugging and planning

---

### Example 6: Low RAM Mode

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    --low-ram
```

**Scenario**: Running on a machine with limited RAM

**What happens**:
- Passes CheckM2 `--lowmem` to the gather script
- Halves DIAMOND RAM (CheckM2 already uses far less RAM than legacy CheckM1)
- Slightly slower quality assessment

---

### Example 7: Fast Mode (Skip CheckM2)

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    --use-bbmap
```

**Scenario**: Need fast results, less concerned about genome quality

**What happens**:
- Uses bbmap statswrapper instead of CheckM2
- Much faster (no completeness/contamination modeling)
- Only basic assembly statistics
- No completeness/contamination filtering

---

### Example 8: Skip Download (Use Pre-Downloaded Genomes)

```bash
# Pre-download genomes manually
mkdir -p results/02_orthophyl_novel/downloads/NovelGenus/genomes_to_keep/
cp /data/genomes/*.fna results/02_orthophyl_novel/downloads/NovelGenus/genomes_to_keep/

# Run wrapper
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --threads 32 \
    --skip-download
```

**Scenario**: You've already downloaded genomes or want to use custom genome sets

**What happens**:
- Skips genome downloading step
- Uses genomes in genomes_to_keep/ directories
- Proceeds with OrthoPhyl analysis

---

### Example 9: Require Specific Genomes to Pass QC

```bash
# Comma-separated list...
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    --must-keep GCF_000001.1,GCF_000002.1

# ...or a file with one accession per line
python orthophyl_pipeline_wrapper.py \
    ... \
    --must-keep required_accessions.txt
```

**Scenario**: Downstream steps depend on particular reference accessions being present.

**What happens**:
- Query/input genomes are run **through the same QC filter** as the downloads.
- If a **query genome fails QC**, the run **aborts** by default with a clear report
  naming the genome and the failed metric (e.g. `completeness=82.0 < MIN_completeness=95`).
  Add `--keep-failing-query` to downgrade this to a warning and force the query in.
- If any `--must-keep` accession is dropped by QC, the run **always aborts** (not
  overridable) — those genomes are required downstream.
- Per-genome failure reasons are written to `qc_removal_reasons.txt` in the download dir.

---

### Example 10: Cap Tree Size for Oversized Taxa (`--max-tree-genomes`)

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --max-tree-genomes 150 \
    --threads 32
```

**Scenario**: A query routes to a novel genus (e.g. `Andreesenella`) that has far
more assemblies on NCBI than OrthoPhyl can put into one tree.

**What happens**:
- The wrapper downloads the **raw** candidate set (`--download-only`, *before* the
  expensive CheckM2 QC).
- If the raw count exceeds `--max-tree-genomes`, MASH partitions the set into
  size-bounded subclades named `Andreesenella_1`, `Andreesenella_2`, …
- A tree is built **only for the subclade containing the query** — QC (CheckM2)
  runs *only* on that subclade's genomes, then OrthoPhyl. The other subclades are
  registered `built=false` and built lazily the first time a future query routes to
  them.
- Future queries are matched to a subclade by **MASH sequence distance** (nearest
  member), since all subclades share one GTDB taxonomy string.

`--taxon` create mode instead builds **all** subclades (there is no single query to
target). See `readmes/README_ADVANCED_FEATURES.md` § *Subclade Partitioning* for the
full flow and caveats.

---

[← Back to Main README](../README.md)
