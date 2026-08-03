# Technical Reference - OrthoPhyl Wrapper v2

[← Back to Main README](../README.v2.md)

---

This document provides comprehensive technical details about the OrthoPhyl Pipeline Wrapper architecture, scripts, configuration, and file formats.

## Table of Contents

1. [Architecture](#architecture)
2. [Script Dependencies](#script-dependencies)
3. [Input/Output Specifications](#inputoutput-specifications)
4. [Configuration Options](#configuration-options)
5. [Checkpoint System](#checkpoint-system)
6. [File Format Specifications](#appendix-file-format-specifications)

---

## Architecture

### OrthoPhyl Pipeline Overview

![OrthoPhyl Workflow](/img/OP2.2.1_workflow.png)

*Figure: OrthoPhyl workflow showing the main pipeline stages. Grey boxes indicate processes. Orange, tan, and purple boxes represent user input, intermediate files, and species tree outputs, respectively. Purple arrows show iterative approaches. The workflow is divided into four main tasks: a) annotate assemblies, clean-up files, and remove identical CDSs. If more than "N" assemblies are being analyzed, b1) identify a subset of diversity-spanning assemblies, b2) pass them through OrthoFinder to generate orthogroups, and b3) expand the OrthoFinder-identified orthogroups to the full dataset of assemblies through iterative HMM searches. c) Align full orthogroup protein sets, generate and trim matching codon alignments, then filter orthogroups by taxon representation. Finally, d) estimate species tree topologies with concatenated codon alignment supermatrices along with a gene tree to species tree consensus method.*

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
  │   └─ assembly_router_multi.cmd_out3.py
  │       ├─ Query all databases
  │       ├─ Find best taxonomic match
  │       └─ Generate routing decisions
  │           ├─► ReLeaf batch (matched)
  │           └─► OrthoPhyl batch (novel)
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
  │       ├─ gather_filter_asms.sh
  │       │   ├─ Download genomes from NCBI
  │       │   ├─ Run CheckM/bbmap QC
  │       │   └─ Filter by quality metrics
  │       ├─ Add query genomes
  │       ├─ OrthoPhyl.sh
  │       │   ├─ Annotate genomes
  │       │   ├─ Run OrthoFinder
  │       │   ├─ Build alignments
  │       │   └─ Infer phylogeny
  │       └─ create_hierarchical_database_v2.py
  │           └─ Create new database entry
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

### Component Interaction

```
orthophyl_pipeline_wrapper.v2.py (Main Orchestrator)
    │
    ├─► assembly_router_multi.cmd_out3.py
    │   └─ Queries: database_dir/*_db/database_config.json
    │
    ├─► ReLeaf.sh
    │   └─ Uses: database_dir/*_db/orthophyl_run/
    │
    ├─► gather_filter_asms.sh
    │   └─ Downloads from: NCBI Datasets API
    │
    ├─► OrthoPhyl.sh
    │   └─ Runs: OrthoFinder, MAFFT, trimAl, IQ-TREE
    │
    ├─► create_hierarchical_database_v2.py
    │   └─ Creates: database_dir/*_db/
    │
    └─► add_releaf_version.py
        └─ Versions: database_dir/*_db/v*_releaf_*/
```

---

---

## Script Dependencies

### 1. assembly_router_multi.cmd_out3.py

**Purpose**: Multi-database assembly router with automatic best-match selection

**Key Classes**:

- **GTDBTaxonomy**: Parses and manipulates GTDB taxonomy strings
  ```python
  tax = GTDBTaxonomy("d__Bacteria;p__Pseudomonadota;...")
  tax.get_rank('p')  # Returns 'Pseudomonadota'
  tax.get_most_specific_rank()  # Returns 's', 'g', 'f', etc.
  ```

- **MultiDatabaseRouter**: Routes assemblies by querying multiple databases
  ```python
  router = MultiDatabaseRouter(
      database_dir=Path("databases/"),
      output_dir=Path("routing_output/"),
      gather_filter_script=Path("utils/gather_filter_asms.sh"),
      threads=8
  )
  ```

**Key Methods**:

- `_load_databases()`: Scans database directory, loads all configs
- `find_matching_databases(query_taxonomy)`: Returns all matching databases, sorted by specificity
- `route_assembly(assembly_path, taxonomy, assembly_id)`: Main routing logic
- `_route_to_releaf()`: Generates ReLeaf decision with command
- `_route_to_orthophyl()`: Generates OrthoPhyl decision with download commands
- `batch_route(input_table)`: Processes multiple assemblies from TSV

**Routing Algorithm**:
```python
def route_assembly(assembly, taxonomy):
    matches = find_matching_databases(taxonomy)
    
    if matches:
        # Use most specific match
        best_match = matches[0]  # Sorted by specificity
        return route_to_releaf(assembly, best_match)
    else:
        # No match found
        return route_to_orthophyl(assembly, taxonomy)
```

**Output Files**:
- `routing_decision_{id}.json` - Machine-readable
- `routing_summary_{id}.txt` - Human-readable
- `batch_routing_summary.txt` - Batch overview

---

### 2. create_hierarchical_database_v2.py

**Purpose**: Create/update hierarchical taxonomy databases from OrthoPhyl runs

**Key Functions**:

- **validate_orthophyl_run(orthophyl_dir)**
  - Checks for required files (HMMs, alignments, trees)
  - Counts genomes
  - Returns validation dictionary

- **create_database_for_run(orthophyl_dir, clade_taxonomy, clade_name, output_dir)**
  - Creates database directory structure
  - Copies/links essential files
  - Generates metadata

- **get_existing_databases(output_dir)**
  - Scans for existing databases
  - Used in update mode to skip duplicates

- **create_master_index(databases, output_dir)**
  - Creates `database_index.json`
  - Generates `database_summary.txt`

**Database Structure Created**:
```
{clade_name}_db/
├── database_config.json       # Metadata
├── phylogeny.nwk              # Species tree
├── genome_list.txt            # Genome identifiers
├── orthophyl_run/             # Symlink to OrthoPhyl output
│   ├── OG_alignmentsToHMM/
│   │   └── hmms_final/        # HMM profiles for ReLeaf
│   ├── phylo_current/
│   │   ├── AlignmentsProts.trm/
│   │   └── AlignmentsCDS.trm/
│   └── FINAL_SPECIES_TREES/
└── README.txt
```

**Modes**:
- **Create**: Error if database exists
- **Update** (`--update`): Skip existing, add new only
- **Force** (`--force`): Rebuild all databases

**Input Format** (orthophyl_runs.tsv):
```tsv
clade_name	orthophyl_dir	clade_taxonomy
Rhizobiaceae	/data/rhizo_run	d__Bacteria;p__Pseudomonadota;...;f__Rhizobiaceae
Escherichia	/data/ecoli_run	d__Bacteria;...;g__Escherichia
```

---

### 3. add_releaf_version.py

**Purpose**: Create new database versions from ReLeaf output

**Key Functions**:

- **validate_releaf_output(releaf_dir)**
  - Checks for updated alignments
  - Verifies new trees exist
  - Returns validation status

- **create_composite_orthophyl_dir(base_version_dir, releaf_output_dir, composite_dir)**
  - Combines unchanged files from base with updates from ReLeaf
  - Creates symlink structure:
    - HMMs → base version (unchanged)
    - Alignments → ReLeaf output (updated)
    - Trees → ReLeaf output (updated)

- **create_releaf_version(db_dir, releaf_output, version_name)**
  - Main versioning function
  - Creates new version directory
  - Updates genome list
  - Links to parent version

**Version Structure**:
```
{clade_name}_db/
├── v1_orthophyl_initial/      # Original OrthoPhyl run
│   ├── database_config.json
│   ├── version_info.json
│   ├── phylogeny.nwk
│   └── orthophyl_run/
├── v2_releaf_2025-01-15/      # First ReLeaf update
│   ├── database_config.json
│   ├── version_info.json
│   ├── phylogeny.nwk
│   ├── base_version → ../v1_orthophyl_initial
│   ├── releaf_output → /path/to/ReLeaf_dir
│   ├── orthophyl_run_composite/  # Merged structure
│   └── orthophyl_run → orthophyl_run_composite
├── v3_releaf_2025-02-01/      # Second ReLeaf update
│   └── ...
└── current → v3_releaf_2025-02-01  # Symlink to latest
```

**version_info.json**:
```json
{
  "version": "v2_releaf_2025-01-15",
  "version_number": 2,
  "version_type": "releaf",
  "created": "2025-01-15T14:30:00",
  "parent_version": "v1_orthophyl_initial",
  "base_genomes": 45,
  "added_genomes": 3,
  "total_genomes": 48,
  "source_type": "ReLeaf",
  "source_path": "/path/to/ReLeaf_dir"
}
```

---

### 4. OrthoPhyl.sh

**Purpose**: Main phylogenomic pipeline for comprehensive analysis

**Key Steps**:

1. **Environment Setup**
   - Loads control files and functions
   - Parses command-line arguments
   - Sets up directory structure

2. **Genome Annotation** (Prodigal)
   ```bash
   prodigal -i genome.fna \
       -a proteins.faa \
       -d cds.fna \
       -f gff -o annotations.gff
   ```

3. **Ortholog Identification** (OrthoFinder)
   ```bash
   orthofinder -f annots_prots/ \
       -t {threads} \
       -a {threads}
   ```

4. **SCO Filtering**
   - Selects single-copy orthologs
   - Filters by presence threshold

5. **Alignment** (MAFFT)
   ```bash
   mafft --auto --thread {threads} \
       OG.fa > OG.aln
   ```

6. **Trimming** (trimAl)
   ```bash
   trimal -in OG.aln \
       -out OG.trm \
       -automated1
   ```

7. **Tree Inference** (IQ-TREE)
   ```bash
   iqtree -s concatenated.fa \
       -p partition_file \
       -m MFP \
       -bb 1000 \
       -nt {threads}
   ```

**Output**: Complete phylogenomic analysis in `store/` directory

---

### 5. ReLeaf.sh

**Purpose**: Add new genomes to existing phylogenies

**Key Steps**:

1. **Load Existing Data**
   - HMM profiles from database
   - Old alignments
   - Old trees

2. **Annotate New Genomes** (Prodigal)
   - Same as OrthoPhyl

3. **HMM Search** (HMMER)
   ```bash
   hmmsearch --tblout results.tbl \
       OG.hmm new_proteins.faa
   ```

4. **Extract Sequences**
   - Pulls matching sequences for each OG

5. **Add to Alignments** (MAFFT)
   ```bash
   mafft --add new_seqs.fa \
       --thread {threads} \
       old_alignment.fa > updated.fa
   ```

6. **Re-trim** (trimAl)
   - Uses same columns as original

7. **Update Trees** (IQ-TREE or FastTree)
   ```bash
   iqtree -s updated.aln \
       -m {model} \
       -nt {threads}
   ```

**Output**: Updated phylogenies in `ReLeaf_dir/`

---

### 6. gather_filter_asms.sh

**Purpose**: Download and filter genomes from NCBI

**Key Functions**:

- **get_NCBI_genomes()**
  ```bash
  datasets download genome taxon {taxon} --dehydrated
  unzip ncbi_dataset.zip
  datasets rehydrate --gzip --directory ./
  ```

- **filter_NCBI_genomes()**
  - Removes duplicate RefSeq/GenBank pairs
  - Keeps highest quality

- **get_stats_with_checkM()**
  ```bash
  checkm lineage_wf \
      --reduced_tree \
      -t {threads} \
      assemblies/ checkM_out/
  ```

- **get_asm_stats()** (bbmap alternative)
  ```bash
  statswrapper.sh in=genome.fna
  ```

- **filter_asm_by_stats()**
  - Applies quality thresholds
  - Removes outliers

- **filter_for_redundancy()**
  - ANI-based clustering
  - Keeps representative genomes

**Quality Filters**:
- Completeness ≥ 95% (CheckM only)
- Contamination ≤ 1.0% (CheckM only)
- Duplication ≤ 2% (CheckM only)
- N50, length, GC within 3 SD of mean

**Output**: `genomes_to_keep/` with high-quality, non-redundant genomes

---

---

## Input/Output Specifications

### Input Files

#### 1. assemblies.tsv (Required)

**Format**: Tab-separated values

**Columns**:
1. `assembly_path` - Full path to genome FASTA
2. `taxonomy` - GTDB-format taxonomy string
3. `assembly_id` - Identifier (optional, defaults to filename)

**Example**:
```tsv
/data/genomes/genome1.fna	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__	VAXT01_genome1
/data/genomes/genome2.fna	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coli	Ecoli_strain123
```

**Notes**:
- Paths can be absolute or relative
- Tilde (~) expansion supported
- Taxonomy must follow GTDB format: `rank__name;rank__name;...`
- Empty rank values allowed: `s__` (species undefined)

#### 2. orthophyl_runs.tsv (Optional, for initial setup)

**Format**: Tab-separated values

**Columns**:
1. `clade_name` - Database identifier
2. `orthophyl_dir` - Path to OrthoPhyl output
3. `clade_taxonomy` - GTDB taxonomy for clade

**Example**:
```tsv
Rhizobiaceae	/data/orthophyl_runs/rhizobiaceae	d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae
Enterobacteriaceae	/data/orthophyl_runs/enterobacteriaceae	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae
```

### Output Structure

```
output_dir/
├── 00_routing/
│   ├── routing_decision_{id}.json
│   ├── routing_summary_{id}.txt
│   ├── batch_routing_summary.txt
│   └── routing_{timestamp}.log
│
├── 01_releaf_only/
│   ├── {database_name}/
│   │   ├── input_genomes/
│   │   │   └── {assembly_id}.fna
│   │   ├── ReLeaf_dir/
│   │   │   ├── annots_prots/
│   │   │   ├── hmm_out/
│   │   │   ├── new_prot_alignments/
│   │   │   ├── new_CDS_alignments/
│   │   │   └── new_trees/
│   │   └── ReLeaf_results/
│   │       └── phylogeny_with_new_genomes.nwk
│   └── ...
│
├── 02_orthophyl_novel/
│   ├── downloads/
│   │   └── {taxon_name}/
│   │       ├── assemblies_datasets_uniq/
│   │       ├── checkM_out/
│   │       ├── assembly_stats.tsv
│   │       ├── genomes_to_keep/
│   │       └── .download_complete
│   └── orthophyl_runs/
│       └── {taxon_name}/
│           ├── genomes/
│           ├── annots_prots/
│           ├── annots_nucls/
│           ├── OG_alignmentsToHMM/
│           ├── phylo_current/
│           └── FINAL_SPECIES_TREES/
│               └── SCO_strict.CDS.iqtree.treefile
│
├── 03_results/
│   ├── trees/
│   │   ├── releaf/
│   │   │   └── {database}_phylogeny.nwk
│   │   └── orthophyl/
│   │       └── {taxon}_phylogeny.nwk
│   └── pipeline_summary.txt
│
├── logs/
│   ├── routing.log
│   ├── releaf_{database}.log
│   ├── download_{taxon}.log
│   ├── orthophyl_{taxon}.log
│   └── database_{taxon}.log
│
├── checkpoints/
│   ├── initialization.flag
│   ├── routing.flag
│   ├── releaf_{database}.flag
│   ├── download_{taxon}.flag
│   ├── orthophyl_{taxon}.flag
│   └── database_{taxon}.flag
│
└── pipeline_status.json
```

### Database Structure

```
database_dir/
├── database_index.json
├── database_summary.txt
├── orthophyl_runs.tsv
│
└── {clade_name}_db/
    ├── database_config.json
    ├── phylogeny.nwk
    ├── genome_list.txt
    ├── README.txt
    │
    ├── v1_orthophyl_initial/
    │   ├── database_config.json
    │   ├── version_info.json
    │   ├── phylogeny.nwk
    │   ├── genome_list.txt
    │   └── orthophyl_run/
    │
    ├── v2_releaf_{date}/
    │   ├── database_config.json
    │   ├── version_info.json
    │   ├── phylogeny.nwk
    │   ├── genome_list.txt
    │   ├── base_version → ../v1_orthophyl_initial
    │   ├── releaf_output → /path/to/ReLeaf_dir
    │   ├── orthophyl_run_composite/
    │   └── orthophyl_run → orthophyl_run_composite
    │
    └── current → v2_releaf_{date}
```

---

---

## Configuration Options

### Command-Line Arguments

```bash
python orthophyl_pipeline_wrapper.v2.py [OPTIONS]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--input FILE` | Input TSV file with assemblies and taxonomies |
| `--database-dir DIR` | Directory containing taxonomy databases |
| `--output-dir DIR` | Output directory for all results |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads N` | 8 | Number of CPU threads to use |
| `--gather-script PATH` | None | Path to gather_filter_asms.sh for genome downloading |
| `--orthophyl-runs FILE` | None | TSV for initial database creation |
| `--resume` | False | Resume from last checkpoint |
| `--skip-download` | False | Skip genome downloading (use existing) |
| `--dry-run` | False | Show commands without executing |
| `-v, --verbose` | 0 | Verbose output (use -v or -vv) |
| `--low-ram` | False | Use CheckM --reduced_tree (low RAM mode) |
| `--use-bbmap` | False | Use bbmap instead of CheckM (faster, less stringent) |

### Verbosity Levels

- **Level 0** (default): Minimal output, logs to files
- **Level 1** (`-v`): Show stdout from subprocesses
- **Level 2** (`-vv`): Show stdout and stderr from subprocesses

### Memory Options

**Standard Mode** (default):
- CheckM with full reference tree
- ~40 GB RAM required

**Low RAM Mode** (`--low-ram`):
- CheckM with reduced tree
- ~16 GB RAM required
- Slightly less accurate

**Fast Mode** (`--use-bbmap`):
- bbmap statswrapper instead of CheckM
- ~4 GB RAM required
- No completeness/contamination filtering
- Much faster

---

---

## Checkpoint System

### How It Works

The wrapper uses a checkpoint system to enable robust resumption after interruptions.

**Checkpoint Files**: Simple flag files in `output_dir/checkpoints/`

```
checkpoints/
├── initialization.flag
├── routing.flag
├── releaf_Rhizobiaceae.flag
├── releaf_Enterobacteriaceae.flag
├── download_NovelGenus.flag
├── orthophyl_NovelGenus.flag
└── database_NovelGenus.flag
```

**Checkpoint Content**: ISO timestamp
```
2025-01-15T14:30:45.123456
```

### Checkpoint Hierarchy

```
initialization
    ↓
routing
    ↓
    ├─► releaf_{database_1}
    ├─► releaf_{database_2}
    │
    └─► For each novel taxon:
        ├─► download_{taxon}
        ├─► orthophyl_{taxon}
        └─► database_{taxon}
```

### Resume Behavior

When `--resume` is specified:

1. **Check Each Phase**:
   ```python
   if checkpoint_exists('phase_name') and resume:
       logger.info("✓ Phase already complete (resuming)")
       return
   ```

2. **Skip Completed Work**:
   - Initialization: Skip if flag exists
   - Routing: Load previous results
   - ReLeaf: Skip completed databases
   - Download: Verify genomes exist
   - OrthoPhyl: Skip if output exists
   - Database: Skip if database created

3. **Validation**:
   - For downloads: Checks for `.download_complete` marker
   - For downloads: Verifies `genomes_to_keep/` has files
   - For OrthoPhyl: Checks for tree file

### Manual Checkpoint Management

**View Checkpoints**:
```bash
ls -lh output_dir/checkpoints/
```

**Remove Specific Checkpoint** (to re-run phase):
```bash
rm output_dir/checkpoints/orthophyl_NovelGenus.flag
```

**Clear All Checkpoints** (start fresh):
```bash
rm -rf output_dir/checkpoints/
```

**Partial Reset** (re-run from routing):
```bash
rm output_dir/checkpoints/routing.flag
rm output_dir/checkpoints/releaf_*.flag
rm output_dir/checkpoints/download_*.flag
rm output_dir/checkpoints/orthophyl_*.flag
rm output_dir/checkpoints/database_*.flag
```

---

---

## Appendix: File Format Specifications

### A. database_config.json

```json
{
  "created": "2025-01-15T10:30:00.123456",
  "version": "1.0",
  "database_type": "hierarchical_taxonomy",
  "clade_name": "Rhizobiaceae",
  "clade_taxonomy": "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae",
  "clade_rank": "f",
  "clade_rank_name": "family",
  "orthophyl_source": "/data/orthophyl_runs/rhizobiaceae",
  "n_genomes": 150,
  "has_hmms": true,
  "has_trees": true,
  "has_alignments": true,
  "available_tree_methods": ["iqtree", "fasttree"],
  "available_data_types": ["CDS", "protein"],
  "validation_details": [
    "Found 1234 HMM/alignment files in hmms_final",
    "Found 1234 alignment files in AlignmentsProts.trm",
    "Found tree: SCO_strict.CDS.iqtree.treefile",
    "Counted 150 genomes from genome_list"
  ]
}
```

### B. database_index.json

```json
{
  "created": "2025-01-15T10:30:00.123456",
  "version": "1.0",
  "n_databases": 5,
  "databases": [
    {
      "clade_name": "Rhizobiaceae",
      "clade_taxonomy": "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae",
      "clade_rank": "f",
      "clade_rank_name": "family",
      "database_dir": "/data/databases/Rhizobiaceae_db",
      "n_genomes": 150
    },
    {
      "clade_name": "Enterobacteriaceae",
      "clade_taxonomy": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae",
      "clade_rank": "f",
      "clade_rank_name": "family",
      "database_dir": "/data/databases/Enterobacteriaceae_db",
      "n_genomes": 200
    }
  ]
}
```

### C. version_info.json

```json
{
  "version": "v2_releaf_2025-01-15",
  "version_number": 2,
  "version_type": "releaf",
  "created": "2025-01-15T14:30:45.123456",
  "parent_version": "v1_orthophyl_initial",
  "base_genomes": 150,
  "added_genomes": 5,
  "total_genomes": 155,
  "source_type": "ReLeaf",
  "source_path": "/data/results/01_releaf_only/Rhizobiaceae/ReLeaf_dir"
}
```

### D. pipeline_status.json

```json
{
  "start_time": "2025-01-15T10:00:00.123456",
  "end_time": "2025-01-15T18:30:45.654321",
  "dry_run": false,
  "phases": {
    "initialization": {
      "status": "complete"
    },
    "routing": {
      "status": "complete",
      "releaf_count": 10,
      "orthophyl_count": 3
    },
    "releaf": {
      "status": "complete",
      "databases_processed": 2
    },
    "orthophyl": {
      "status": "complete",
      "taxa_processed": 3
    },
    "aggregation": {
      "status": "complete",
      "releaf_trees": 2,
      "orthophyl_trees": 3
    }
  },
  "summary": {
    "total_assemblies": 13,
    "releaf_assemblies": 10,
    "orthophyl_assemblies": 3,
    "databases_used": 2,
    "databases_created": 3
  }
}
```

---

[← Back to Main README](../README.v2.md)
