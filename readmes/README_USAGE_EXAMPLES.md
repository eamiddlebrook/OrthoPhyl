# Usage Examples - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---


### Example 1: Basic Run with Existing Databases

```bash
python orthophyl_pipeline_wrapper.v2.py \
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
python orthophyl_pipeline_wrapper.v2.py \
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
python orthophyl_pipeline_wrapper.v2.py \
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
python orthophyl_pipeline_wrapper.v2.py \
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
python orthophyl_pipeline_wrapper.v2.py \
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
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    --low-ram
```

**Scenario**: Running on a machine with limited RAM

**What happens**:
- Uses CheckM --reduced_tree option
- Reduces RAM from ~40 GB to ~16 GB
- Slightly less accurate quality assessment

---

### Example 7: Fast Mode (Skip CheckM)

```bash
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    --use-bbmap
```

**Scenario**: Need fast results, less concerned about genome quality

**What happens**:
- Uses bbmap statswrapper instead of CheckM
- Much faster (no marker gene analysis)
- Only basic assembly statistics
- No completeness/contamination filtering

---

### Example 8: Skip Download (Use Pre-Downloaded Genomes)

```bash
# Pre-download genomes manually
mkdir -p results/02_orthophyl_novel/downloads/NovelGenus/genomes_to_keep/
cp /data/genomes/*.fna results/02_orthophyl_novel/downloads/NovelGenus/genomes_to_keep/

# Run wrapper
python orthophyl_pipeline_wrapper.v2.py \
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

---

[← Back to Main README](../README.md)
