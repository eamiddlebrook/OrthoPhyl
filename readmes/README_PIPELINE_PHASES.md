# Pipeline Phases - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---

This document provides detailed information about each phase of the OrthoPhyl Pipeline Wrapper execution.

## Table of Contents

1. [Phase 1: Initialization](#phase-1-initialization)
2. [Phase 2: Assembly Routing](#phase-2-assembly-routing)
3. [Phase 3a: ReLeaf Route (Matched Databases)](#phase-3a-releaf-route-matched-databases)
4. [Phase 3b: OrthoPhyl Route (Novel Taxa)](#phase-3b-orthophyl-route-novel-taxa)
5. [Phase 4: Results Aggregation](#phase-4-results-aggregation)

---

## Phase 1: Initialization

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
   - `assembly_router/assembly_router.py`
   - `assembly_router/create_hierarchical_database.py`
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

## Phase 2: Assembly Routing

**Purpose**: Determine the appropriate pipeline for each assembly

**Script**: `assembly_router.py`

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

## Phase 3a: ReLeaf Route (Matched Databases)

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

## Phase 3b: OrthoPhyl Route (Novel Taxa)

**Purpose**: Create comprehensive phylogenies for novel taxa

**Process**: Multi-stage workflow for each novel taxon

### Stage 1: Download Genomes

**Script**: `gather_filter_asms.sh`

**Purpose**: Download and filter high-quality genomes from NCBI

**Steps**:

1. **Download from NCBI**
   ```bash
   utils/gather_filter_asms.sh \
       {taxon_name} \
       {output_dir} \
       {threads} \
       [--lowmem]        # Low RAM mode (halves DIAMOND RAM)
       [--use-bbmap]     # Skip CheckM2
   ```

2. **NCBI Datasets API**
   - Downloads all genomes for specified taxon
   - Handles RefSeq and GenBank assemblies
   - Retry logic for server failures

3. **Quality Control** (CheckM2 or bbmap):
   
   **Option A: CheckM2** (default, more stringent):
   ```bash
   checkm2 predict \
       --lowmem \  # Optional: halves DIAMOND RAM
       -t {threads} \
       --input assemblies/ --output-directory checkM_out/
   ```
   - Assesses completeness and contamination via DIAMOND + pretrained ML models
     (no reference tree or pplacer)
   - Default filters:
     - Completeness ≥ 95%
     - Contamination ≤ 1.0%
     - Duplication ≤ 2% (placeholder under CheckM2 — no marker-copy metric, so
       this column is 0.00 and the filter is a no-op)

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

### Stage 2: Add Query Genomes

**Purpose**: Combine downloaded genomes with query assemblies

```python
# Copy query genomes to filtered set
for assembly in query_assemblies:
    shutil.copy(assembly, genomes_to_keep/)
```

**Result**: Complete genome set for phylogenetic analysis

### Stage 3: Run OrthoPhyl

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

### Stage 4: Create Database Entry

**Script**: `create_hierarchical_database.py`

**Purpose**: Make OrthoPhyl output queryable for future runs

**Process**:

1. **Update orthophyl_runs.tsv**
   ```tsv
   {taxon_name}	{orthophyl_output}	{taxonomy}
   ```

2. **Run Database Creator**
   ```bash
   python create_hierarchical_database.py \
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

## Phase 4: Results Aggregation

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

[← Back to Main README](../README.md)
