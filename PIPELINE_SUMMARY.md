# OrthoPhyl Pipeline Wrapper - Complete System Summary

## Overview

We've built a complete automated pipeline system for phylogenetic placement of bacterial genomes that intelligently routes assemblies to either ReLeaf (existing databases) or OrthoPhyl (novel taxa) workflows.

---

## Core Components

### 1. Assembly Router (`assembly_router_multi.cmd_out3.py`)

**Purpose:** Automatically determines the best pipeline for each genome assembly

**Key Features:**
- Multi-database querying
- Hierarchical taxonomy matching (most specific match wins)
- Tree method/data type detection (iqtree, fasttree; CDS, PROT)
- Path expansion (`~` support)
- Quote stripping from taxonomy strings
- Integration with genome download scripts

**Decision Logic:**
```
For each assembly + taxonomy:
  Query all databases
  ↓
  Match found?
  ├─→ YES: Route to ReLeaf
  │   - Use most specific matching database
  │   - Detect compatible tree methods
  │   - Generate ReLeaf command
  │
  └─→ NO: Route to OrthoPhyl
      - Determine download taxon (genus/family)
      - Generate 3-step workflow:
        1. Download reference genomes (gather_filter_asms.sh)
        2. Add query genomes + run OrthoPhyl
        3. Create new database
```

### 2. Database Creator (`create_hierarchical_database_v2.py` / `v3.py`)

**Purpose:** Create queryable databases from OrthoPhyl runs

**Key Features:**
- Validates OrthoPhyl outputs (HMMs, alignments, trees)
- Detects available tree methods and data types
- Creates indexed database structure
- Supports incremental updates (`--update` mode)
- Force rebuild option (`--force`)

**Database Structure:**
```
database_dir/
├── Clade_db/
│   ├── database_config.json       # Metadata
│   ├── orthophyl_run@ → /path/    # Symlink to OrthoPhyl output
│   ├── phylogeny.nwk              # Representative tree
│   └── genome_list.txt            # Genome inventory
└── database_index.json            # Master index
```

### 3. ReLeaf Version Manager (`add_releaf_version.py`)

**Purpose:** Create versioned database snapshots after ReLeaf runs

**Key Features:**
- Version tracking (v1, v2, v3, ...)
- Composite directory synthesis (combines base + ReLeaf outputs)
- Symlink-based (avoids file duplication)
- Automatic version numbering
- 'current' pointer management

**Version Structure:**
```
Clade_db/
├── v1_orthophyl_initial/          # Original OrthoPhyl run
│   ├── version_info.json
│   ├── database_config.json
│   ├── orthophyl_run@ → /original/path
│   └── genome_list.txt (100 genomes)
│
├── v2_releaf_2025-01-15/          # After first ReLeaf
│   ├── version_info.json          # Tracks: +10 genomes, parent=v1
│   ├── base_version@ → ../v1
│   ├── releaf_output@ → /ReLeaf_dir
│   ├── orthophyl_run_composite/   # Synthesized structure
│   │   ├── OG_alignmentsToHMM@ → v1/...
│   │   ├── phylo_current/
│   │   │   ├── AlignmentsProts@ → ReLeaf/new_prot_alignments
│   │   │   ├── AlignmentsCDS@ → ReLeaf/new_CDS_alignments
│   │   │   ├── trimmed_columns@ → v1/...
│   │   │   └── SCO_sets@ → v1/...
│   │   └── FINAL_SPECIES_TREES/
│   │       └── trees@ → ReLeaf/new_trees/
│   ├── orthophyl_run@ → orthophyl_run_composite
│   └── genome_list.txt (110 genomes)
│
├── current@ → v2_releaf_2025-01-15  # Points to latest
└── orthophyl_run@ → current/orthophyl_run  # Backward compatibility
```

**Why This Works:**
- ReLeaf needs specific files from OrthoPhyl runs
- Some files unchanged (HMMs, trimming schemes)
- Some files updated (alignments, trees)
- Composite directory provides correct structure for next ReLeaf run
- Symlinks avoid copying large files

### 4. Pipeline Wrapper (`orthophyl_pipeline_wrapper.py`)

**Purpose:** Orchestrate the entire workflow from raw assemblies to final phylogenies

**Workflow:**
```
1. Initialization
   - Validate dependencies
   - Load/create databases
   - Create output structure

2. Assembly Routing
   - Run assembly_router_multi.cmd_out3.py
   - Parse routing decisions
   - Group by route type

3a. ReLeaf Route (Matched Assemblies)
    - Group by database
    - Prepare input directories
    - Run ReLeaf for each database
    - Create new database version ← NEW!
    - Update 'current' pointer

3b. OrthoPhyl Route (Novel Taxa)
    - Group by download taxon
    - Download reference genomes (gather_filter_asms.sh)
    - Add query genomes to download set ← KEY STEP!
    - Run OrthoPhyl (queries included from start)
    - Create v1 database

4. Results Aggregation
   - Collect all phylogenies
   - Generate summary reports
   - Verify queries in trees
```

**Key Features:**
- Checkpoint system (resume support)
- Dry-run mode (preview without execution)
- Automatic version creation after ReLeaf
- Comprehensive logging
- Status tracking (JSON)

---

## Critical Design Decisions

### 1. Query Inclusion in OrthoPhyl

**Problem:** Original design had queries added AFTER OrthoPhyl, requiring follow-up ReLeaf

**Solution:** Add query genomes BEFORE running OrthoPhyl

```bash
# OLD (inefficient):
Download refs → Run OrthoPhyl → Create DB → Run ReLeaf to place queries

# NEW (efficient):
Download refs → Add queries → Run OrthoPhyl (all genomes) → Create DB
```

**Benefits:**
- One less step
- Queries fully integrated (not just placed)
- Query genes in orthogroups
- Immediate complete results

### 2. Version Control for Databases

**Problem:** ReLeaf modifies alignments/trees, creating new state. How to track?

**Solution:** Version-based database system

**Benefits:**
- Clear provenance
- Can rollback to previous versions
- Run ReLeaf on different versions simultaneously
- Audit trail of database evolution

### 3. Composite Directory Synthesis

**Problem:** ReLeaf needs complete OrthoPhyl structure, but we want to avoid file duplication

**Solution:** Create composite directory with symlinks

**Benefits:**
- No file duplication
- Clear which files from which version
- Fast to create
- Easy to verify structure

---

## Usage Examples

### Initial Setup

```bash
# 1. Create databases from existing OrthoPhyl runs
cat > orthophyl_runs.tsv << EOF
Gaiellales	/data/orthophyl/gaiellales	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales
Escherichia	/data/orthophyl/ecoli	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;g__Escherichia
EOF

python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir databases/
```

### Run Complete Pipeline

```bash
# 2. Prepare assembly list with taxonomies
cat > assemblies.tsv << EOF
/data/genomes/GCF_001.fna	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae	GCF_001
/data/genomes/GCF_002.fna	d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae;g__Rhizobium	GCF_002
EOF

# 3. Run pipeline
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32

# Output:
# - GCF_001 → ReLeaf (matched Gaiellales) → results/01_releaf_only/Gaiellales/
# - GCF_002 → OrthoPhyl (novel Rhizobium) → results/02_orthophyl_novel/Rhizobium/
# - New database created: databases/Rhizobium_db/v1_orthophyl_initial/
# - Gaiellales database updated: v2_releaf_2025-01-15/
```

### Dry Run (Preview)

```bash
# See what would happen without running anything
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --dry-run
```

### Resume After Failure

```bash
# Pipeline uses checkpoints - can resume from where it stopped
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --resume
```

### Add New Database Entry

```bash
# After manual OrthoPhyl run, add to databases
echo "Novitaxon	/data/orthophyl/novitaxon	d__Bacteria;p__..." >> orthophyl_runs.tsv

python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir databases/ \
    --update  # Skip existing, add new only
```

### Manually Create ReLeaf Version

```bash
# If ReLeaf run manually, add as new version
python assembly_router/add_releaf_version.py \
    --database-dir databases/Gaiellales_db/ \
    --releaf-output /path/to/store/ReLeaf_dir/

# List all versions
python assembly_router/add_releaf_version.py \
    --database-dir databases/Gaiellales_db/ \
    --list-versions
```

---

## File Organization

```
OrthoPhyl/
├── orthophyl_pipeline_wrapper.py        # Main orchestrator
│
├── assembly_router/
│   ├── assembly_router_multi.cmd_out3.py    # Routing logic
│   ├── create_hierarchical_database_v2.py   # Database creator
│   ├── create_hierarchical_database_v3.py   # (future: built-in versioning)
│   ├── add_releaf_version.py                # Version manager
│   ├── RELEAF_DATABASE_DESIGN.md            # Design document
│   └── WORKFLOW_DESIGN.md                   # Workflow specification
│
├── utils/
│   └── gather_filter_asms.sh               # Genome download + QC
│
├── OrthoPhyl.sh                            # Core OrthoPhyl
├── ReLeaf.sh                               # Core ReLeaf
├── script_lib/
│   ├── functions.sh
│   └── functions_addem.sh                  # ReLeaf functions
│
└── databases/                              # Database directory
    ├── database_index.json
    ├── Clade1_db/
    │   ├── v1_orthophyl_initial/
    │   ├── v2_releaf_2025-01-15/
    │   └── current@ → v2_releaf_2025-01-15
    └── Clade2_db/
        └── v1_orthophyl_initial/
```

---

## Key Improvements from Original Design

### Before:
1. Manual database creation
2. Manual routing decisions
3. Query genomes added after OrthoPhyl (required extra ReLeaf)
4. No version control
5. No automated workflow

### After:
1. ✅ Automated database creation and updates
2. ✅ Intelligent automatic routing
3. ✅ Query genomes included in initial OrthoPhyl run
4. ✅ Full version control for databases
5. ✅ Complete automated workflow with resume support
6. ✅ Dry-run mode for testing
7. ✅ Tree method/data type compatibility checking
8. ✅ Comprehensive logging and status tracking

---

## Future Enhancements

### Potential v4 Features:

1. **Database merging**: Combine multiple databases at the same rank
2. **Version comparison**: Tools to diff versions and see what changed
3. **Automated testing**: Test database integrity and ReLeaf compatibility
4. **Web interface**: View database hierarchy and versions
5. **Metadata enrichment**: Add isolation info, publication links, etc.
6. **Quality metrics**: Track alignment quality, tree support values
7. **Parallel ReLeaf**: Run multiple ReLeaf jobs simultaneously
8. **Cloud integration**: Support for cloud storage (S3, etc.)

---

## Dependencies

### Python Packages:
- Python 3.7+ (uses `subprocess.run(text=True)`)
- Standard library only (no external packages required)

### External Tools:
- **OrthoPhyl**: Core phylogenomic pipeline
- **ReLeaf**: Add sequences to existing trees
- **Prodigal**: Gene prediction
- **MAFFT**: Multiple sequence alignment
- **IQTree** or **FastTree**: Phylogenetic inference
- **HMMER**: HMM searches
- **CheckM**: Genome quality assessment (for gather_filter_asms.sh)
- **NCBI datasets CLI**: Genome downloading

### Conda Environment:
```bash
conda create -n orthophyl \
    -c bioconda -c conda-forge \
    orthofinder prodigal mafft iqtree fasttree hmmer \
    checkm-genome bbmap entrez-direct ncbi-datasets-cli
```

---

## Testing

### Test Dataset:
```bash
# Small test with known genomes
cat > test_assemblies.tsv << EOF
~/test/GCF_000005845.fna	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia	GCF_000005845
EOF

# Dry run first
python orthophyl_pipeline_wrapper.py \
    --input test_assemblies.tsv \
    --database-dir test_databases/ \
    --output-dir test_results/ \
    --threads 4 \
    --dry-run

# Then real run
python orthophyl_pipeline_wrapper.py \
    --input test_assemblies.tsv \
    --database-dir test_databases/ \
    --output-dir test_results/ \
    --threads 4
```

---

## Citation

If you use this pipeline, please cite:

**OrthoPhyl:**
```
Earl A Middlebrook, Robab Katani, Jeanne M Fair
OrthoPhyl - Streamlining large scale, orthology-based phylogenomic studies of bacteria at broad evolutionary scales
G3 Genes|Genomes|Genetics, 2024; jkae119
https://doi.org/10.1093/g3journal/jkae119
```

---

## Support

For issues or questions:
- GitHub: https://github.com/eamiddlebrook/OrthoPhyl
- Open an issue with:
  - Command used
  - Error messages
  - Log files (`output_dir/logs/`)
  - Pipeline status JSON

---

**Last Updated:** January 2025  
**Pipeline Version:** 3.0  
**Authors:** HGTool/OrthoPhyl Project
