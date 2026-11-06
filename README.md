# OrthoPhyl Pipeline Wrapper v2

## Quick Navigation

📚 **[Jump to Detailed Documentation](#detailed-documentation)** | 🚀 **[Quick Start](#quick-start)** | 🔍 **[Troubleshooting](readmes/README_TROUBLESHOOTING.md)**

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Quick Start](#quick-start)
4. [Common Usage Examples](#common-usage-examples)
5. [Essential Configuration](#essential-configuration)
6. [Common Issues & Solutions](#common-issues--solutions)
7. [Performance Considerations](#performance-considerations)
8. [Detailed Documentation](#detailed-documentation)
9. [Citation & Support](#citation--support)

---

## Overview

The **OrthoPhyl Pipeline Wrapper v2** (`orthophyl_pipeline_wrapper.py`) is an automated orchestration system that intelligently routes genome assemblies to the appropriate phylogenetic placement pipeline. It seamlessly integrates two complementary approaches:

- **ReLeaf Route**: For assemblies matching existing databases (fast, adds to pre-computed phylogenies)
- **OrthoPhyl Route**: For novel taxa requiring new phylogenetic analyses (comprehensive, creates new databases)

### Key Features

✅ **Intelligent Routing**: Automatically determines the best pipeline for each assembly  
✅ **Database Management**: Creates and versions hierarchical taxonomy databases  
✅ **Checkpoint/Resume**: Robust recovery from interruptions  
✅ **Quality Control**: Integrated genome filtering with CheckM2 or bbmap  
✅ **Batch Processing**: Handles multiple assemblies efficiently  
✅ **Flexible Configuration**: Supports various tree methods and data types  

### What It Does

1. **Routes** assemblies by querying multiple taxonomy databases
2. **Executes** ReLeaf for assemblies matching existing databases
3. **Downloads** related genomes for novel taxa from NCBI
4. **Runs** OrthoPhyl to create comprehensive phylogenies for novel taxa
5. **Creates** new database entries for future use
6. **Aggregates** all results into a unified output structure

---

## Architecture

### Workflow Diagram

```
INPUT: assemblies.tsv (assembly_path, taxonomy, [id])
  │
  ├─► PHASE 1: INITIALIZATION
  │   ├─ Create directory structure
  │   ├─ Validate dependencies
  │   └─ Load/create databases
  │
  ├─► PHASE 2: ASSEMBLY ROUTING
  │   └─ assembly_router.py  (--placement subclade|backbone, default subclade)
  │       ├─ Query all databases
  │       ├─ Find best taxonomic match
  │       │   └─ Megatree ties (backbone + its subclades share one taxonomy)
  │       │       are broken by --placement, then MASH distance
  │       └─ Generate routing decisions
  │           ├─► ReLeaf batch (matched, built)
  │           ├─► OrthoPhyl batch (novel taxon)
  │           └─► Subclade-build batch (matched an unbuilt --megatree-lazy
  │                                      subclade; build it, then ReLeaf)
  │
  ├─► PHASE 3A: RELEAF ROUTE (Matched Databases)
  │   └─ For each database:
  │       ├─ ReLeaf.sh
  │       │   ├─ Add genomes to existing phylogeny
  │       │   └─ Generate updated trees
  │       └─ add_releaf_version.py
  │           └─ Create new database version
  │
  ├─► PHASE 3B: ORTHOPHYL ROUTE (Novel Taxa)
  │   └─ For each taxon:
  │       ├─ gather_filter_asms.sh --download-only
  │       │   ├─ Download raw genomes from NCBI (pre-QC)
  │       │   └─ Stage query genomes into the raw set
  │       ├─ subsample_genomes.py   (default, if raw count > --max-tree-genomes)
  │       │   └─ Diverse MASH subsample down to --subsample-size (one tree)
  │       ├─ subclade_partition.py + megatree_graft.py   (opt-in: --megatree)
  │       │   ├─ MASH-partition raw set into <Taxon>_1, <Taxon>_2, …
  │       │   ├─ Build a full tree per subclade + a backbone tree
  │       │   │   (or, with --megatree-lazy, only build subclades holding a
  │       │   │    query; register the rest as built=false placeholders)
  │       │   └─ Graft subclade trees onto backbone → one merged megatree
  │       ├─ gather_filter_asms.sh --qc-only   (per subclade)
  │       │   ├─ Run CheckM2/bbmap QC (subclade members + queries)
  │       │   └─ Filter by quality; abort if a required genome fails
  │       ├─ OrthoPhyl.sh
  │       │   ├─ Annotate genomes
  │       │   ├─ Run OrthoFinder
  │       │   ├─ Build alignments
  │       │   └─ Infer phylogeny
  │       └─ OP_database_tool.py
  │           └─ Create new database entry
  │
  ├─► PHASE 3C: SUBCLADE-BUILD ROUTE (Lazy Subclades, --megatree-lazy only)
  │   └─ For each unbuilt subclade a query matched:
  │       ├─ QC + OrthoPhyl.sh on the subclade's own raw members (no query)
  │       ├─ OP_database_tool.py --force
  │       │   └─ Promote the built=false placeholder to built=true
  │       └─ ReLeaf.sh the waiting query assemblies onto the new tree
  │
  └─► PHASE 4: RESULTS AGGREGATION
      ├─ Collect all trees
      ├─ Generate summary report
      └─ Create unified output structure

OUTPUT: results/03_results/
  ├─ trees/
  │   ├─ releaf/
  │   └─ orthophyl/
  └─ pipeline_summary.txt
```

**📖 For detailed phase descriptions, see [Pipeline Phases Documentation](readmes/README_PIPELINE_PHASES.md)**

---

## Quick Start

### Prerequisites

- Python 3.9+
- The OrthoPhyl toolchain, either via the **Apptainer/Singularity image** (recommended) or a **conda/mamba environment** (see [Installation](#installation) below)
- Access to taxonomy databases (or `orthophyl_runs.tsv` for initial setup)

### Installation

The wrapper orchestrates `OrthoPhyl.sh` and `ReLeaf.sh`, so it needs the full OrthoPhyl toolchain (OrthoFinder, IQ-TREE, MAFFT, HMMER, Prodigal, trimal, etc.). You have two options. The full walkthrough — including per-tool versions, ASTRAL/catfasta2phyml/Alignment_Assessment setup, and `control_file.paths` editing — lives in the [OrthoPhyl & ReLeaf Guide](README_OrthoPhyl_ReLeaf.md#GettingStarted); container specifics (building the image yourself, bind-mount semantics, moving databases between host and container) live in [readmes/README_APPTAINER.md](readmes/README_APPTAINER.md).

**Option A — Apptainer/Singularity (recommended, avoids dependency management)**

Grab the prebuilt container (~2.7 GB; `apptainer` and `singularity` are drop-in replacements for each other in every command below):

```bash
singularity_images=~/singularity_images/
mkdir -p ${singularity_images}
cd ${singularity_images}
apptainer pull library://earlyevol/default/orthophyl
mv orthophyl_latest.sif OrthoPhyl.v3.3.sif
```

Verify the image runs (writes 12 trees under `FINAL_SPECIES_TREES/`; `-s` output dir is **required**). Note the default `%runscript` just drops you into a shell — `exec` the script directly, don't `run` the image with pipeline args:

```bash
apptainer exec ${singularity_images}/OrthoPhyl.v3.3.sif \
    bash /opt/gits/OrthoPhyl/OrthoPhyl.sh -T TESTER_fasttest -s ./tester_fasttest_output -t 4
```

**Option B — conda/mamba environment**

Clone the repo (this pulls large test files, so it takes a minute) and create the `orthophyl` environment from the versioned spec file (replace `XXX` with the current version present in the repo, e.g. `orthophyl_env.2.2.1.yml`):

```bash
git clone https://github.com/eamiddlebrook/OrthoPhyl.git
cd OrthoPhyl

# mamba recommended; conda works with identical commands
mamba env create -n orthophyl -f orthophyl_env.XXX.yml
mamba activate orthophyl
```

A few external tools (ASTRAL, catfasta2phyml, Alignment_Assessment) are installed separately and pointed to via `control_file.paths` — see the [Manual Install section](README_OrthoPhyl_ReLeaf.md#ManualInstall) for those steps and troubleshooting (R, GNU parallel, `libnsl`, conda init).

### Basic Usage

```bash
# Simple run with existing databases
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --threads 32

# With genome downloading for novel taxa
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32

# Initial setup (create databases from scratch)
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --orthophyl-runs orthophyl_runs.tsv \
    --threads 32

# Build a database from genomes already on disk (QC runs by default)
python orthophyl_pipeline_wrapper.py \
    --genome-dir /data/my_isolates/ \
    --clade-name MyIsolates \
    --database-dir databases/ \
    --output-dir results/ \
    --threads 32
```

### Input File Format

**assemblies.tsv** (tab-separated):
```
/path/to/genome1.fnad__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__genome1_id
/path/to/genome2.fnad__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coligenome2_id
```

Columns:
1. **assembly_path**: Full path to genome FASTA file
2. **taxonomy**: GTDB-format taxonomy string
3. **assembly_id** (optional): Identifier (defaults to filename stem)

---

## Common Usage Examples

### Example 1: Standard Run with Existing Databases

```bash
python orthophyl_pipeline_wrapper.py \
    --input my_assemblies.tsv \
    --database-dir /data/phylo_databases/ \
    --output-dir /results/my_project/ \
    --threads 32
```

**Use case**: You have pre-built databases and want to place new assemblies

### Example 2: Enable Genome Downloading for Novel Taxa

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32
```

**Use case**: Assemblies may include novel taxa requiring NCBI genome downloads

### Example 3: Resume from Checkpoint

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --resume \
    --threads 32
```

**Use case**: Pipeline was interrupted; resume from last checkpoint

### Example 4: Low-Memory Mode

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --low-ram --use-bbmap \
    --threads 16
```

**Use case**: Running on systems with limited RAM (<64 GB)

### Example 5: Build a Database from Genomes Already on Disk

```bash
python orthophyl_pipeline_wrapper.py \
    --genome-dir /data/my_isolates/ \
    --clade-name Pseudomonas \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32
```

**Use case**: You have your own assemblies (unpublished isolates, a curated set) and want
a searchable, routable OrthoPhyl database out of them, without going through NCBI. See
**[readmes/README_USAGE.md](readmes/README_USAGE.md)** for QC-skipping, naming clades
that aren't formally assigned by NCBI, and taxonomy-routability details.

**📝 For worked examples across all three modes (batch/taxon/local-genome), see
[readmes/README_USAGE.md](readmes/README_USAGE.md)**

---

## Essential Configuration

Flags are grouped below by subject, matching `--help`'s layout.

### Mode selection (mutually exclusive, one required)

| Argument | Description | Default |
|----------|-------------|---------|
| `--input` | Input TSV: assembly_path, taxonomy, [id]. Batch mode. | - |
| `--taxon` | Taxon name for NCBI auto-gather mode (e.g. "Methylorubrum"). | - |
| `--genome-dir` | Directory of genomes already on disk (local mode; requires `--clade-name`). | - |

### Core

| Argument | Description | Default |
|----------|-------------|---------|
| `--database-dir` | Directory containing `*_db` databases (required) | - |
| `--output-dir` | Output directory for all results | `<database-dir>/.pipeline_runs/<taxon>_<timestamp>` |

### Run control

| Argument | Description | Default |
|----------|-------------|---------|
| `--resume` | Resume from last checkpoint | False |
| `--dry-run` | Show what would be executed without running anything | False |
| `--skip-download` | Skip genome downloading (use existing genomes) | False |
| `--update-existing` | Update existing database with new assemblies (taxon mode only) | False |

### Performance

| Argument | Description | Default |
|----------|-------------|---------|
| `--threads` | Number of threads | 8 |
| `--low-ram` | Reduced-memory CheckM2 mode (`--lowmem`) | False |
| `--use-bbmap` | Use bbmap instead of CheckM2 for genome stats (faster, less RAM, no completeness/contamination filtering) | False |
| `--ani-shortlist` | OrthoFinder MASH-shortlist size (`-n`); forced so small databases still get HMMs built | 20 |

### Gather / QC

| Argument | Description | Default |
|----------|-------------|---------|
| `--gather-script` | Path to genome download script | `utils/gather_filter_asms.sh` (next to the wrapper) |
| `--orthophyl-runs` | TSV for initial database creation | None |
| `--must-keep` | Accessions that must survive QC, or the run aborts with a per-metric report | None |
| `--keep-failing-query` | Let a failing query genome through with a warning instead of aborting | False |
| `--skip-qc` | Skip CheckM2 QC on `--genome-dir` genomes | False |

### Taxon mode (`--taxon`)

| Argument | Description | Default |
|----------|-------------|---------|
| `--taxon-rank` | Taxonomic rank for the `--taxon` query (species/genus/family/order/class/phylum) | auto-detect |

### Local genome-ingest mode (`--genome-dir`)

| Argument | Description | Default |
|----------|-------------|---------|
| `--clade-name` | Names the clade/database (required with `--genome-dir`) | - |
| `--clade-taxonomy` | Full GTDB taxonomy string, used verbatim, for clades that don't resolve against NCBI | None |
| `--clade-rank` | GTDB rank letter (`d`..`s`) an unresolvable `--clade-name` is attached at | `g` (genus) |

### Large-taxon handling

| Argument | Description | Default |
|----------|-------------|---------|
| `--max-tree-genomes` | Single-tree ceiling before diverse subsampling (or `--megatree`) kicks in | 2000 |
| `--subsample-size` | Target genome count for the diverse MASH subsample | 500 |
| `--max-total-genomes` | Guardrail for the opt-in megatree partitioner (O(n²) distance array) | 25000 |

### Megatree (opt-in, `--megatree`)

| Argument | Description | Default |
|----------|-------------|---------|
| `--megatree` | Partition an oversized taxon into subclades + a backbone instead of subsampling to one tree. Supported with `--genome-dir` (not combinable with `--skip-qc`). | False |
| `--backbone-reps` | Diverse representatives each subclade contributes to the backbone | 5 |
| `--subclade-size` | Per-subclade genome ceiling | 150 |
| `--conflict-min-support` | Support threshold for flagging a backbone/subclade bipartition conflict | 90 |
| `--megatree-lazy` | Defer building subclades with no query at partition time; build on demand later. Implies `--megatree`. | False |
| `--megatree-hmm-reuse` | Give every subclade the backbone's orthogroup HMMs instead of its own OrthoFinder run. Implies `--megatree`. | False |
| `--megatree-hmm-reuse-skip-leftover` | With `--megatree-hmm-reuse`, drop genes unmatched by the backbone's HMMs instead of clustering them. Implies `--megatree`. | False |
| `--placement` | Tie-break for a query matching both a subclade and the backbone (`subclade`/`backbone`) | `subclade` |

### Geographic figures (opt-in, `--geo`)

Auto-generates a country-annotated tree, a world map of sample origins, and
lines joining the top-K closest cross-country genome pairs by patristic
distance after each published tree, landing under `<results-dir>/figures/`.
Routed queries (batch mode) or `--geo-highlight`-listed genomes are
highlighted on the figures with their own nearest-neighbor map and TSVs.
Figure generation is non-fatal; a missing `matplotlib` is checked fatally at
startup. `--genome-dir` mode has no NCBI location metadata, so `--geo-tsv
genome_id<TAB>location` is required there.

| Argument | Description | Default |
|----------|-------------|---------|
| `--geo` | Generate the geo tree + map figures after each published tree. | False |
| `--geo-tsv` | `genome_id<TAB>location` TSV, merged on top of (or substituting for) NCBI `geo_loc_name` metadata. Implies `--geo`. | - |
| `--geo-top-k` | Number of closest cross-country genome pairs to draw lines for. | 100 |
| `--geo-format` | Output image format (`png`/`svg`/`pdf`). | `png` |
| `--geo-highlight` | Extra genomes to highlight (comma list or file). Implies `--geo`. | - |
| `--geo-query-top-k` | Nearest-neighbor lines per highlighted genome on the query map. | 25 |

**🔧 For the full writeup (query highlighting, all outputs, per-flag details), see
[Advanced Features § Geographic Figures](readmes/README_ADVANCED_FEATURES.md#12-geographic-figures---geo)**

### Low importance

| Argument | Description | Default |
|----------|-------------|---------|
| `-v`, `--verbose` | `-v` shows subprocess stdout, `-vv` shows stdout and stderr | 0 |

**🔧 For complete configuration reference, see [Technical Reference](readmes/README_TECHNICAL_REFERENCE.md#configuration-options)**

---

## Common Issues & Solutions

### 1. "No databases found" Error

**Problem**: Pipeline cannot find any database directories

**Solution**:
```bash
# Check database directory structure
ls -la databases/
# Should contain *_db directories with database_config.json files

# If empty, initialize with orthophyl_runs.tsv
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --orthophyl-runs orthophyl_runs.tsv
```

### 2. Checkpoint Resume Not Working

**Problem**: `--resume` flag doesn't skip completed phases

**Solution**:
```bash
# If you passed --output-dir explicitly, check it has checkpoint files
ls results/checkpoints/

# If missing, checkpoints may have been deleted
# Re-run without --resume to start fresh
```

If `--output-dir` was omitted on the original run, `--resume` alone (still
without `--output-dir`) automatically reuses the newest matching prior run
dir under `<database-dir>/.pipeline_runs/<taxon>_<timestamp>/` -- no need to
dig up the exact timestamped path. If no matching run is found, it logs a
warning and starts a fresh run instead of resuming nothing. Pass
`--output-dir` explicitly to pin down a specific run.

### 3. Out of Memory During CheckM2

**Problem**: CheckM2 fails with memory errors

**Solution**:
```bash
# CheckM2 (DIAMOND + ML models, no reference tree/pplacer) already uses far less
# RAM than legacy CheckM1. Use --low-ram to halve DIAMOND RAM, or --use-bbmap to skip.
python orthophyl_pipeline_wrapper.py \
    --low-ram --use-bbmap \
    --threads 16  # Reduce threads
```

### 4. NCBI Download Failures

**Problem**: Genome downloads timeout or fail

**Solution**:
- Check internet connection
- NCBI servers may be temporarily down
- Use `--resume` to retry failed downloads
- Consider downloading genomes manually

### 5. ReLeaf Fails with "HMM profiles not found"

**Problem**: Database missing HMM profiles

**Solution**:
```bash
# Check database structure
ls databases/MyDatabase_db/orthophyl_run/OG_alignmentsToHMM/hmms_final/

# If missing, database may be incomplete
# Rebuild database from OrthoPhyl output
```

**🔍 For comprehensive troubleshooting, see [Troubleshooting Guide](readmes/README_TROUBLESHOOTING.md)**

---

## Performance Considerations

### Resource Requirements

**Minimum**:
- 16 GB RAM
- 8 CPU cores
- 100 GB disk space

**Recommended**:
- 64 GB RAM
- 32 CPU cores
- 500 GB disk space (for large datasets)

### Optimization Tips

1. **Thread Allocation**:
   ```bash
   # Use 80% of available cores
   --threads 32  # on 40-core system
   ```

2. **Low-Memory Systems**:
   ```bash
   --low-ram --use-bbmap
   ```

3. **Large Datasets** (>100 assemblies):
   - Split input into batches
   - Process sequentially
   - Merge results manually

4. **Network Optimization**:
   - Run downloads during off-peak hours
   - Use local genome cache if available

### Runtime Estimates

| Dataset Size | ReLeaf Route | OrthoPhyl Route (with download) |
|--------------|--------------|----------------------------------|
| 1-10 assemblies | 1-3 hours | 6-12 hours |
| 10-50 assemblies | 3-8 hours | 12-24 hours |
| 50-100 assemblies | 8-16 hours | 24-48 hours |

*Estimates based on 32-core system with 64 GB RAM*

---

## Detailed Documentation

### 🧬 Core Pipelines (OrthoPhyl & ReLeaf)

- **[OrthoPhyl & ReLeaf Guide](README_OrthoPhyl_ReLeaf.md)** - Full documentation for the underlying pipelines the wrapper orchestrates
  - Installation (Singularity & manual/conda)
  - Running OrthoPhyl and ReLeaf directly
  - Command-line reference and examples
  - Notes, known errors, and citation

### 📖 Understanding the Pipeline

- **[Pipeline Phases](readmes/README_PIPELINE_PHASES.md)** - Detailed walkthrough of all pipeline phases
  - Phase 1: Initialization
  - Phase 2: Assembly Routing
  - Phase 3a: ReLeaf Route (Matched Databases)
  - Phase 3b: OrthoPhyl Route (Novel Taxa)
  - Phase 3c: Subclade-Build Route (Lazy Subclades)
  - Phase 4: Results Aggregation
  - Routing Workflow Diagram: flowchart of `--taxon` create vs. `--input` query
    routing across taxon scale (subsample / `--megatree` / `--megatree-lazy` / `--placement`)

- **[Technical Reference](readmes/README_TECHNICAL_REFERENCE.md)** - Comprehensive technical documentation
  - Architecture & Component Interaction
  - Script Dependencies
  - Input/Output Specifications
  - Configuration Options
  - Checkpoint System
  - File Format Specifications

### 🚀 Running & Configuring

- **[Usage Guide](readmes/README_USAGE.md)** - Worked examples for all three modes
  - Batch mode: standard runs, initial setup, resume, low-memory, QC requirements
  - Taxon mode: create/update databases directly from a taxon name
  - Local genome-ingest mode: build a database from genomes already on disk, QC
    by default (`--skip-qc` to opt out), naming a clade NCBI doesn't recognize
    (`--clade-name`, `--clade-taxonomy`), taxonomy provenance and routability

- **[Advanced Features](readmes/README_ADVANCED_FEATURES.md)** - Power user features
  - Custom routing logic
  - Database versioning
  - Parallel batch processing
  - Custom quality filters
  - Integration with HPC schedulers
  - Performance tuning
  - Debugging modes

### 🐳 Containers & Testing

- **[Apptainer/Singularity Containers](readmes/README_APPTAINER.md)** - Container usage guide
  - Pulling from Sylabs Cloud, or building the image yourself
  - Bind-mount semantics (identity binds) and moving databases between host and container
  - Running in containers, including the test suite
  - HPC integration

### 🔍 Reference & Recovery

- **[Troubleshooting Guide](readmes/README_TROUBLESHOOTING.md)** - Comprehensive problem-solving
  - Common errors and solutions
  - Database issues
  - Memory problems
  - Network failures
  - Checkpoint recovery
  - Log file analysis
  - Debug mode usage

---

## Citation & Support

### Citation

If you use OrthoPhyl in your research, please cite:

> Earl A Middlebrook, Robab Katani, Jeanne M Fair, OrthoPhyl—streamlining large-scale, orthology-based phylogenomic studies of bacteria at broad evolutionary scales, *G3 Genes|Genomes|Genetics*, Volume 14, Issue 8, August 2024, jkae119, https://doi.org/10.1093/g3journal/jkae119

### Support

- **Issues**: Report bugs via GitHub Issues
- **Questions**: Contact the development team
- **Documentation**: See detailed guides in `readmes/` directory

### License

OrthoPhyl is distributed under the GNU General Public License v3.0. See [GPLv3.pdf](GPLv3.pdf) for details.

---

**📚 [Back to Top](#orthophyl-pipeline-wrapper-v2)** | **🚀 [Quick Start](#quick-start)** | **📖 [Detailed Documentation](#detailed-documentation)**
