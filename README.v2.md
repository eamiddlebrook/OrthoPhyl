# OrthoPhyl Pipeline Wrapper v2 - Complete Walkthrough

## Table of Contents
1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Quick Start](#quick-start)
4. [Pipeline Phases](#pipeline-phases)
5. [Script Dependencies](#script-dependencies)
6. [Input/Output Specifications](#inputoutput-specifications)
7. [Configuration Options](#configuration-options)
8. [Usage Examples](#usage-examples)
9. [Checkpoint System](#checkpoint-system)
10. [Troubleshooting](#troubleshooting)
11. [Advanced Features](#advanced-features)

---

## Overview

The **OrthoPhyl Pipeline Wrapper v2** (`orthophyl_pipeline_wrapper.v2.py`) is an automated orchestration system that intelligently routes genome assemblies to the appropriate phylogenetic placement pipeline. It seamlessly integrates two complementary approaches:

- **ReLeaf Route**: For assemblies matching existing databases (fast, adds to pre-computed phylogenies)
- **OrthoPhyl Route**: For novel taxa requiring new phylogenetic analyses (comprehensive, creates new databases)

### Key Features

✅ **Intelligent Routing**: Automatically determines the best pipeline for each assembly  
✅ **Database Management**: Creates and versions hierarchical taxonomy databases  
✅ **Checkpoint/Resume**: Robust recovery from interruptions  
✅ **Quality Control**: Integrated genome filtering with CheckM or bbmap  
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

## Quick Start

### Prerequisites

- Python 3.7+
- OrthoPhyl conda environment
- Access to taxonomy databases (or orthophyl_runs.tsv for initial setup)

### Basic Usage

```bash
# Simple run with existing databases
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --threads 32

# With genome downloading for novel taxa
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32

# Initial setup (create databases from scratch)
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --orthophyl-runs orthophyl_runs.tsv \
    --threads 32
```

### Input File Format

**assemblies.tsv** (tab-separated):
```
/path/to/genome1.fna	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__	genome1_id
/path/to/genome2.fna	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coli	genome2_id
```

Columns:
1. **assembly_path**: Full path to genome FASTA file
2. **taxonomy**: GTDB-format taxonomy string
3. **assembly_id** (optional): Identifier (defaults to filename stem)

---

## Pipeline Phases

### Phase 1: Initialization

**Purpose**: Set up directory structure and validate dependencies

**Actions**:
1. Creates output directory structure:
   ```
   output_dir/
   ├── 00_routing/          # Routing decisions
   ├── 01_releaf_only/      # ReLeaf outputs
   ├── 02_orthophyl_novel/  # OrthoPhyl outputs
   ├── 03_results/          # Final aggregated results
   ├── logs/                # All log files
   └── checkpoints/         # Resume flags
   ```

2. Validates required scripts:
   - `assembly_router/assembly_router_multi.cmd_out3.py`
   - `assembly_router/create_hierarchical_database_v2.py`
   - `assembly_router/add_releaf_version.py`
   - `OrthoPhyl.sh`
   - `ReLeaf.sh`
   - `utils/gather_filter_asms.sh` (optional)

3. Initializes or validates databases:
   - Checks for `database_dir/database_index.json`
   - If missing and `--orthophyl-runs` provided, creates initial databases
   - Loads database metadata

**Checkpoint**: `initialization.flag`

---

### Phase 2: Assembly Routing

**Purpose**: Determine the appropriate pipeline for each assembly

**Script**: `assembly_router_multi.cmd_out3.py`

**Process**:

1. **Load All Databases**
   - Scans `database_dir/` for `*_db` directories
   - Reads `database_config.json` from each database
   - Parses taxonomy information

2. **For Each Assembly**:
   - Parse query taxonomy (GTDB format)
   - Query all databases for taxonomic matches
   - Find most specific match (species > genus > family > ...)
   - Generate routing decision

3. **Routing Logic**:
   ```python
   if assembly_taxonomy matches database_taxonomy:
       → ReLeaf Route
       - Use existing database
       - Fast phylogenetic placement
   else:
       → OrthoPhyl Route
       - Download related genomes
       - Run full phylogenetic analysis
       - Create new database
   ```

4. **Output Files** (per assembly):
   - `routing_decision_{assembly_id}.json` - Machine-readable decision
   - `routing_summary_{assembly_id}.txt` - Human-readable summary
   - `batch_routing_summary.txt` - Overall batch summary

**Example Routing Decision (ReLeaf)**:
```json
{
  "pipeline": "ReLeaf",
  "reason": "Taxonomy matches Rhizobiaceae at family level",
  "assembly": "/path/to/genome.fna",
  "assembly_id": "genome_id",
  "query_taxonomy": "d__Bacteria;p__Pseudomonadota;...;f__Rhizobiaceae",
  "matched_database": "Rhizobiaceae",
  "matched_rank": "family",
  "database_dir": "/path/to/databases/Rhizobiaceae_db",
  "database_genomes": 150,
  "tree_method": "iqtree",
  "tree_data": "CDS"
}
```

**Example Routing Decision (OrthoPhyl)**:
```json
{
  "pipeline": "OrthoPhyl",
  "reason": "Novel taxonomy not represented in any database",
  "assembly": "/path/to/genome.fna",
  "assembly_id": "genome_id",
  "query_taxonomy": "d__Bacteria;...;g__NovelGenus;s__",
  "download_taxonomy": "d__Bacteria;...;g__NovelGenus",
  "download_rank": "genus",
  "download_value": "NovelGenus",
  "suggestion": "Create new database for genus 'NovelGenus'"
}
```

**Checkpoint**: `routing.flag`

---

### Phase 3a: ReLeaf Route (Matched Databases)

**Purpose**: Add assemblies to existing phylogenies using ReLeaf

**Script**: `ReLeaf.sh`

**Process**:

1. **Group by Database**
   - Assemblies are grouped by matched database
   - Each database processed independently

2. **For Each Database**:

   a. **Prepare Input**
      ```bash
      # Copy assemblies to input directory
      01_releaf_only/{database_name}/input_genomes/
      ```

   b. **Run ReLeaf**
      ```bash
      ./ReLeaf.sh \
          --store {database_dir}/orthophyl_run \
          --input_genomes input_genomes/ \
          -t {threads} \
          --tree_method iqtree \
          --TREE_DATA CDS
      ```

   c. **ReLeaf Steps** (internal):
      - Annotate new genomes (Prodigal)
      - Search against HMM profiles from database
      - Extract orthologous sequences
      - Add to existing alignments
      - Re-trim alignments
      - Update phylogenetic trees (IQ-TREE or FastTree)

   d. **Create Database Version** (if `add_releaf_version.py` exists):
      ```bash
      python add_releaf_version.py \
          --database-dir {database_dir} \
          --releaf-output {releaf_output_dir}
      ```
      - Creates versioned database (e.g., `v2_releaf_2025-01-15`)
      - Links updated alignments and trees
      - Updates genome list
      - Maintains backward compatibility

3. **Output Structure**:
   ```
   01_releaf_only/{database_name}/
   ├── input_genomes/           # Query genomes
   ├── ReLeaf_dir/              # ReLeaf working directory
   │   ├── annots_prots/        # Protein annotations
   │   ├── hmm_out/             # HMM search results
   │   ├── new_prot_alignments/ # Updated alignments
   │   ├── new_CDS_alignments/  # Updated CDS alignments
   │   └── new_trees/           # Updated phylogenies
   └── ReLeaf_results/
       └── phylogeny_with_new_genomes.nwk
   ```

**Checkpoint**: `releaf_{database_name}.flag`

---

### Phase 3b: OrthoPhyl Route (Novel Taxa)

**Purpose**: Create comprehensive phylogenies for novel taxa

**Process**: Multi-stage workflow for each novel taxon

#### Stage 1: Download Genomes

**Script**: `gather_filter_asms.sh`

**Purpose**: Download and filter high-quality genomes from NCBI

**Steps**:

1. **Download from NCBI**
   ```bash
   utils/gather_filter_asms.sh \
       {taxon_name} \
       {output_dir} \
       {threads} \
       [--reduced_tree]  # Low RAM mode
       [--use-bbmap]     # Skip CheckM
   ```

2. **NCBI Datasets API**
   - Downloads all genomes for specified taxon
   - Handles RefSeq and GenBank assemblies
   - Retry logic for server failures

3. **Quality Control** (CheckM or bbmap):
   
   **Option A: CheckM** (default, more stringent):
   ```bash
   checkm lineage_wf \
       --reduced_tree \  # Optional: low RAM
       -t {threads} \
       assemblies/ checkM_out/
   ```
   - Assesses completeness and contamination
   - Uses marker genes
   - Default filters:
     - Completeness ≥ 95%
     - Contamination ≤ 1.0%
     - Duplication ≤ 2%

   **Option B: bbmap statswrapper** (faster, less RAM):
   ```bash
   statswrapper.sh in=genome.fna
   ```
   - Basic assembly statistics only
   - No completeness/contamination filtering
   - Faster for large datasets

4. **Assembly Statistics Filtering**
   - N50, genome length, GC content
   - Removes outliers (mean ± 3 SD)

5. **Redundancy Filtering**
   - Removes RefSeq/GenBank duplicates
   - Keeps highest quality assembly

6. **Output**:
   ```
   02_orthophyl_novel/downloads/{taxon_name}/
   ├── assemblies_datasets_uniq/  # Downloaded assemblies
   ├── checkM_out/                # Quality metrics
   ├── assembly_stats.tsv         # Statistics
   └── genomes_to_keep/           # Filtered, high-quality genomes
   ```

**Checkpoint**: `download_{taxon_name}.flag`

#### Stage 2: Add Query Genomes

**Purpose**: Combine downloaded genomes with query assemblies

```python
# Copy query genomes to filtered set
for assembly in query_assemblies:
    shutil.copy(assembly, genomes_to_keep/)
```

**Result**: Complete genome set for phylogenetic analysis

#### Stage 3: Run OrthoPhyl

**Script**: `OrthoPhyl.sh`

**Purpose**: Comprehensive phylogenomic analysis

**Command**:
```bash
./OrthoPhyl.sh \
    -g {genomes_to_keep}/ \
    -s {output_dir} \
    -t {threads} \
    -p iqtree \
    -o CDS
```

**OrthoPhyl Pipeline Steps**:

1. **Genome Annotation** (Prodigal)
   - Predicts protein-coding genes
   - Extracts protein and CDS sequences

2. **Ortholog Identification** (OrthoFinder)
   - All-vs-all DIAMOND search
   - MCL clustering
   - Identifies orthogroups (OGs)

3. **Single-Copy Ortholog (SCO) Selection**
   - Filters for genes present in ≥ X% of genomes
   - Ensures phylogenetic signal

4. **Multiple Sequence Alignment** (MAFFT)
   - Aligns each SCO independently
   - Protein and/or CDS alignments

5. **Alignment Trimming** (trimAl)
   - Removes poorly aligned regions
   - Improves phylogenetic inference

6. **Phylogenetic Inference** (IQ-TREE)
   - Model selection (ModelFinder)
   - Maximum likelihood tree
   - Bootstrap support (optional)
   - Partitioned analysis (one partition per gene)

7. **Output Structure**:
   ```
   02_orthophyl_novel/orthophyl_runs/{taxon_name}/
   ├── genomes/                    # Input genomes
   ├── annots_prots/               # Protein annotations
   ├── annots_nucls/               # CDS annotations
   ├── annots_prots.fixed/
   │   └── OrthoFinder/
   │       └── Results_ortho/      # OrthoFinder results
   ├── OG_alignmentsToHMM/
   │   └── hmms_final/             # HMM profiles
   ├── phylo_current/
   │   ├── AlignmentsProts.trm/    # Trimmed protein alignments
   │   ├── AlignmentsCDS.trm/      # Trimmed CDS alignments
   │   └── SpeciesTree/            # Tree inference
   └── FINAL_SPECIES_TREES/
       └── SCO_strict.CDS.iqtree.treefile  # Final tree
   ```

**Checkpoint**: `orthophyl_{taxon_name}.flag`

#### Stage 4: Create Database Entry

**Script**: `create_hierarchical_database_v2.py`

**Purpose**: Make OrthoPhyl output queryable for future runs

**Process**:

1. **Update orthophyl_runs.tsv**
   ```tsv
   {taxon_name}	{orthophyl_output}	{taxonomy}
   ```

2. **Run Database Creator**
   ```bash
   python create_hierarchical_database_v2.py \
       --input orthophyl_runs.tsv \
       --output-dir {database_dir} \
       --update
   ```

3. **Database Structure Created**:
   ```
   database_dir/{taxon_name}_db/
   ├── database_config.json      # Metadata
   ├── phylogeny.nwk             # Species tree
   ├── genome_list.txt           # Genome IDs
   ├── orthophyl_run/            # Symlink to OrthoPhyl output
   └── README.txt                # Usage instructions
   ```

4. **database_config.json**:
   ```json
   {
     "created": "2025-01-15T10:30:00",
     "version": "1.0",
     "database_type": "hierarchical_taxonomy",
     "clade_name": "NovelGenus",
     "clade_taxonomy": "d__Bacteria;...;g__NovelGenus",
     "clade_rank": "g",
     "clade_rank_name": "genus",
     "orthophyl_source": "/path/to/orthophyl_output",
     "n_genomes": 45,
     "has_hmms": true,
     "has_trees": true
   }
   ```

5. **Update Master Index**
   - Adds new database to `database_index.json`
   - Updates `database_summary.txt`

**Checkpoint**: `database_{taxon_name}.flag`

---

### Phase 4: Results Aggregation

**Purpose**: Collect and organize all results

**Actions**:

1. **Collect Trees**
   ```
   03_results/trees/
   ├── releaf/
   │   ├── Rhizobiaceae_phylogeny.nwk
   │   └── Enterobacteriaceae_phylogeny.nwk
   └── orthophyl/
       ├── NovelGenus_phylogeny.nwk
       └── AnotherTaxon_phylogeny.nwk
   ```

2. **Generate Summary Report**
   ```
   03_results/pipeline_summary.txt
   ```
   
   Example content:
   ```
   ======================================================================
   ORTHOPHYL PIPELINE - SUMMARY REPORT
   ======================================================================
   
   Input File: assemblies.tsv
   Database Directory: databases/
   Output Directory: results/
   
   RESULTS:
   ----------------------------------------------------------------------
   
   ReLeaf Route (matched databases): 3 databases
     - Rhizobiaceae
     - Enterobacteriaceae
     - Pseudomonadaceae
   
   OrthoPhyl Route (novel taxa): 2 taxa
     - NovelGenus (new database created)
     - AnotherTaxon (new database created)
   
   OUTPUT LOCATIONS:
     Trees: results/03_results/trees/
     Logs: results/logs/
     Routing decisions: results/00_routing/
   
   ======================================================================
   All assemblies have been placed in phylogenetic trees!
   ======================================================================
   ```

**Checkpoint**: `aggregation.flag`

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

## Usage Examples

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

## Troubleshooting

### Common Issues and Solutions

#### 1. "No databases found in {database_dir}"

**Cause**: Database directory is empty or improperly formatted

**Solutions**:
```bash
# Option A: Provide orthophyl_runs.tsv for initial creation
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --orthophyl-runs orthophyl_runs.tsv

# Option B: Manually create databases first
python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir databases/
```

---

#### 2. "Routing failed. Check log: routing.log"

**Cause**: assembly_router_multi.cmd_out3.py encountered an error

**Debug**:
```bash
# Check the log
cat output_dir/logs/routing.log

# Common issues:
# - Invalid taxonomy format
# - Missing assembly files
# - Corrupted database configs

# Test routing manually
python assembly_router/assembly_router_multi.cmd_out3.py \
    --assembly test.fna \
    --taxonomy "d__Bacteria;p__Pseudomonadota;..." \
    --database-dir databases/ \
    --output-dir test_routing/
```

---

#### 3. "Genome download failed for {taxon}"

**Cause**: NCBI server issues, network problems, or invalid taxon name

**Solutions**:
```bash
# Check the download log
cat output_dir/logs/download_{taxon}.log

# Test download manually
utils/gather_filter_asms.sh "Escherichia" test_download/ 8

# If taxon name is wrong, check NCBI Taxonomy:
# https://www.ncbi.nlm.nih.gov/taxonomy

# If NCBI is down, use --skip-download and provide genomes manually
mkdir -p output_dir/02_orthophyl_novel/downloads/{taxon}/genomes_to_keep/
cp /path/to/genomes/*.fna output_dir/02_orthophyl_novel/downloads/{taxon}/genomes_to_keep/

python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir output_dir/ \
    --skip-download \
    --resume
```

---

#### 4. "OrthoPhyl failed for {taxon}"

**Cause**: Various issues in phylogenomic pipeline

**Debug**:
```bash
# Check the log
cat output_dir/logs/orthophyl_{taxon}.log

# Common issues:
# - Too few genomes (need at least 4)
# - Annotation failures
# - OrthoFinder errors
# - Insufficient memory

# Test OrthoPhyl manually
./OrthoPhyl.sh \
    -g genomes/ \
    -s test_output/ \
    -t 8 \
    -p iqtree \
    -o CDS
```

---

#### 5. "ReLeaf failed for {database}"

**Cause**: Issues adding genomes to existing phylogeny

**Debug**:
```bash
# Check the log
cat output_dir/logs/releaf_{database}.log

# Common issues:
# - Missing HMM profiles in database
# - Corrupted alignments
# - Tree inference failures

# Verify database structure
ls -R databases/{database}_db/current/orthophyl_run/

# Test ReLeaf manually
./ReLeaf.sh \
    --store databases/{database}_db/current/orthophyl_run \
    --input_genomes test_genomes/ \
    -t 8 \
    --tree_method iqtree \
    --TREE_DATA CDS
```

---

#### 6. "CheckM failed" or "Out of memory"

**Cause**: CheckM requires significant RAM (~40 GB)

**Solutions**:
```bash
# Option A: Use low RAM mode
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --low-ram

# Option B: Skip CheckM entirely (faster, less stringent)
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --use-bbmap

# Option C: Pre-filter genomes manually and use --skip-download
```

---

#### 7. "Database creation failed for {taxon}"

**Cause**: Issues creating database from OrthoPhyl output

**Debug**:
```bash
# Check the log
cat output_dir/logs/database_{taxon}.log

# Verify OrthoPhyl output is complete
ls output_dir/02_orthophyl_novel/orthophyl_runs/{taxon}/FINAL_SPECIES_TREES/

# Test database creation manually
python assembly_router/create_hierarchical_database_v2.py \
    --input test_runs.tsv \
    --output-dir databases/ \
    --update
```

---

#### 8. "Query {assembly_id} NOT found in tree"

**Cause**: Query genome was filtered out during OrthoPhyl analysis

**Possible Reasons**:
- Genome quality too low
- Too divergent from other genomes
- Annotation failed

**Solutions**:
```bash
# Check if genome was annotated
ls output_dir/02_orthophyl_novel/orthophyl_runs/{taxon}/annots_prots/{assembly_id}.faa

# Check OrthoFinder results
grep {assembly_id} output_dir/02_orthophyl_novel/orthophyl_runs/{taxon}/annots_prots.fixed/OrthoFinder/Results_ortho/Orthogroups/Orthogroups.txt

# If genome is too divergent, consider:
# 1. Running separate OrthoPhyl analysis
# 2. Using broader taxonomic group for download
```

---

### Log File Locations

All logs are in `output_dir/logs/`:

| Log File | Content |
|----------|---------|
| `routing.log` | Assembly routing decisions |
| `releaf_{database}.log` | ReLeaf execution for specific database |
| `download_{taxon}.log` | Genome downloading and filtering |
| `orthophyl_{taxon}.log` | OrthoPhyl phylogenomic analysis |
| `database_{taxon}.log` | Database creation |
| `releaf_version_{database}.log` | Database versioning |

---

## Advanced Features

### 1. Dry Run Mode

**Purpose**: Preview what would happen without executing

**Usage**:
```bash
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --dry-run \
    -vv
```

**What It Does**:
- Creates directory structure
- Shows all commands that would run
- Parses input files
- No actual computation
- Useful for:
  - Debugging
  - Planning resource allocation
  - Verifying input formats

**Output Example**:
```
[DRY RUN] Would run assembly routing
[DRY RUN] Would run ReLeaf for Rhizobiaceae
[DRY RUN] Would download genomes for NovelGenus
[DRY RUN] Would run OrthoPhyl for NovelGenus
[DRY RUN] Would create database for NovelGenus
```

---

### 2. Verbose Logging

**Levels**:

**Level 0** (default):
```bash
python orthophyl_pipeline_wrapper.v2.py --input assemblies.tsv ...
```
- Minimal console output
- All details in log files
- Best for production runs

**Level 1** (`-v`):
```bash
python orthophyl_pipeline_wrapper.v2.py --input assemblies.tsv ... -v
```
- Shows stdout from subprocesses
- stderr still goes to log files
- Good for monitoring progress

**Level 2** (`-vv`):
```bash
python orthophyl_pipeline_wrapper.v2.py --input assemblies.tsv ... -vv
```
- Shows stdout and stderr
- Maximum verbosity
- Best for debugging

---

### 3. Custom Genome Sets

**Scenario**: You want to use specific genomes instead of automatic NCBI download

**Steps**:

1. **Create genome directory**:
   ```bash
   mkdir -p results/02_orthophyl_novel/downloads/MyTaxon/genomes_to_keep/
   ```

2. **Add your genomes**:
   ```bash
   cp /my/genomes/*.fna results/02_orthophyl_novel/downloads/MyTaxon/genomes_to_keep/
   ```

3. **Run with --skip-download**:
   ```bash
   python orthophyl_pipeline_wrapper.v2.py \
       --input assemblies.tsv \
       --database-dir databases/ \
       --output-dir results/ \
       --skip-download
   ```

**Use Cases**:
- Custom genome collections
- Pre-filtered genomes
- Local genome databases
- Avoiding NCBI download limits

---

### 4. Partial Runs

**Scenario**: Only want to run specific phases

**Method**: Use checkpoints strategically

**Example: Only run routing**:
```bash
# Run full pipeline
python orthophyl_pipeline_wrapper.v2.py --input assemblies.tsv ...

# Examine routing results
cat output_dir/00_routing/batch_routing_summary.txt

# Stop here if you just wanted routing decisions
```

**Example: Skip routing, start from ReLeaf**:
```bash
# Manually create routing decisions
# Then create checkpoint
touch output_dir/checkpoints/routing.flag

# Run with --resume
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir output_dir/ \
    --resume
```

---

### 5. Parallel Execution

**Built-in Parallelization**:
- ReLeaf databases processed sequentially (each uses `--threads`)
- OrthoPhyl taxa processed sequentially (each uses `--threads`)
- Within each task, tools use multiple threads

**Manual Parallelization**:

Split input file and run multiple instances:

```bash
# Split assemblies.tsv
split -l 10 assemblies.tsv batch_

# Run in parallel (different output dirs)
python orthophyl_pipeline_wrapper.v2.py \
    --input batch_aa \
    --database-dir databases/ \
    --output-dir results_batch1/ \
    --threads 16 &

python orthophyl_pipeline_wrapper.v2.py \
    --input batch_ab \
    --database-dir databases/ \
    --output-dir results_batch2/ \
    --threads 16 &

wait

# Merge results
mkdir -p results_merged/03_results/trees/
cp results_batch*/03_results/trees/*/*.nwk results_merged/03_results/trees/
```

---

### 6. Database Versioning

**Automatic Versioning**: Enabled by default when `add_releaf_version.py` exists

**Version Structure**:
```
database_dir/Rhizobiaceae_db/
├── v1_orthophyl_initial/      # Original
├── v2_releaf_2025-01-15/      # After first ReLeaf run
├── v3_releaf_2025-02-01/      # After second ReLeaf run
└── current → v3_releaf_2025-02-01
```

**List Versions**:
```bash
python assembly_router/add_releaf_version.py \
    --database-dir databases/Rhizobiaceae_db/ \
    --list-versions
```

**Rollback to Previous Version**:
```bash
cd databases/Rhizobiaceae_db/
rm current
ln -s v2_releaf_2025-01-15 current
```

**Manual Version Creation**:
```bash
python assembly_router/add_releaf_version.py \
    --database-dir databases/Rhizobiaceae_db/ \
    --releaf-output /path/to/ReLeaf_dir/ \
    --version-name v4_custom_update
```

---

### 7. Custom Tree Methods and Data Types

**Available in Databases**: Specified in `database_config.json`

```json
{
  "available_tree_methods": ["iqtree", "fasttree"],
  "available_data_types": ["CDS", "protein"]
}
```

**Router Behavior**:
- Prefers `iqtree` and `CDS` if available
- Falls back to other options if needed
- Generates appropriate ReLeaf commands

**Override for New OrthoPhyl Runs**:

Edit `OrthoPhyl.sh` command in wrapper (line 695-696):
```python
cmd = [
    str(self.orthophyl_script),
    '-g', str(input_dir),
    '-s', str(output_dir),
    '-t', str(self.threads),
    '-p', 'fasttree',  # Change to 'fasttree'
    '-o', 'protein'    # Change to 'protein'
]
```

---

### 8. Integration with HPC Systems

**SLURM Example**:

```bash
#!/bin/bash
#SBATCH --job-name=orthophyl_wrapper
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=orthophyl_%j.log

# Load environment
module load conda
conda activate orthophyl

# Run wrapper
python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir /scratch/databases/ \
    --output-dir /scratch/results_${SLURM_JOB_ID}/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads ${SLURM_CPUS_PER_TASK} \
    -v

# Copy results to permanent storage
cp -r /scratch/results_${SLURM_JOB_ID}/03_results/ /home/user/results/
```

**PBS Example**:

```bash
#!/bin/bash
#PBS -N orthophyl_wrapper
#PBS -l nodes=1:ppn=32
#PBS -l mem=128gb
#PBS -l walltime=48:00:00
#PBS -o orthophyl.log
#PBS -e orthophyl.err

cd $PBS_O_WORKDIR

conda activate orthophyl

python orthophyl_pipeline_wrapper.v2.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32 \
    -v
```

---

## Performance Considerations

### Resource Requirements

**Per Assembly (ReLeaf Route)**:
- CPU: 1-8 cores
- RAM: 4-16 GB
- Time: 10-60 minutes
- Disk: 1-5 GB

**Per Taxon (OrthoPhyl Route)**:
- CPU: 8-64 cores
- RAM: 16-128 GB (depends on genome count and CheckM mode)
- Time: 2-24 hours
- Disk: 10-100 GB

**Genome Download**:
- Network: Depends on NCBI
- Disk: 50-500 MB per genome
- Time: 1-10 minutes per genome

### Optimization Tips

1. **Use --use-bbmap for large datasets**
   - 10x faster than CheckM
   - 1/10th the RAM
   - Trade-off: less stringent QC

2. **Pre-download genomes**
   - Avoids NCBI rate limits
   - More reliable for large batches
   - Use `--skip-download`

3. **Increase --threads**
   - Most tools scale well to 32-64 cores
   - Diminishing returns beyond 64

4. **Use SSD for working directory**
   - Significant speedup for I/O-heavy steps
   - Especially important for OrthoFinder

5. **Split large batches**
   - Run multiple instances in parallel
   - Each with subset of assemblies

---

## Citation

If you use this pipeline wrapper, please cite:

**OrthoPhyl**:
- Middlebrook, E.A., et al. (2024). OrthoPhyl: Automated phylogenomic analysis pipeline. GitHub: https://github.com/eamiddlebrook/OrthoPhyl

**Key Dependencies**:
- OrthoFinder: Emms, D.M. and Kelly, S. (2019). Genome Biology, 20:238
- IQ-TREE: Nguyen, L.T., et al. (2015). Molecular Biology and Evolution, 32:268-274
- CheckM: Parks, D.H., et al. (2015). Genome Research, 25:1043-1055
- MAFFT: Katoh, K. and Standley, D.M. (2013). Molecular Biology and Evolution, 30:772-780

---

## Support and Contact

**Issues**: https://github.com/eamiddlebrook/OrthoPhyl/issues

**Documentation**: https://github.com/eamiddlebrook/OrthoPhyl

**Email**: [Contact information]

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

**End of README.v2.md**

*Last Updated: 2025-01-15*
*Version: 2.0*
