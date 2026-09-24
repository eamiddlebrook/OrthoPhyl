# Pipeline Phases - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---

This document provides detailed information about each phase of the OrthoPhyl Pipeline Wrapper execution.

## Table of Contents

1. [Phase 1: Initialization](#phase-1-initialization)
2. [Phase 2: Assembly Routing](#phase-2-assembly-routing)
3. [Phase 3a: ReLeaf Route (Matched Databases)](#phase-3a-releaf-route-matched-databases)
4. [Phase 3b: OrthoPhyl Route (Novel Taxa)](#phase-3b-orthophyl-route-novel-taxa)
5. [Phase 3c: Subclade-Build Route (Lazy Subclades)](#phase-3c-subclade-build-route-lazy-subclades)
6. [Phase 4: Results Aggregation](#phase-4-results-aggregation)

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
       if matched_database.is_subclade and not matched_database.built:
           → OrthoPhyl_subclade_build Route
           - Lazily-registered megatree subclade (--megatree-lazy); no tree yet
           - Build it from its own raw members, then ReLeaf the query onto it
       else:
           → ReLeaf Route
           - Use existing database
           - Fast phylogenetic placement
   else:
       → OrthoPhyl Route
       - Download related genomes
       - Run full phylogenetic analysis
       - Create new database
   ```

   A megatree's backbone and its subclades share one taxonomy string, so a
   query can match several databases at the same specificity. That tie is
   broken by `--placement` (`subclade` default, or `backbone`), then by MASH
   distance among the tied subclades (nearest wins; ties broken by ascending
   `clade_name`) — see `readmes/README_ADVANCED_FEATURES.md` § *Placing a
   query*.

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

**Example Routing Decision (OrthoPhyl_subclade_build, `--megatree-lazy` only)**:
```json
{
  "pipeline": "OrthoPhyl_subclade_build",
  "reason": "Query nearest to subclade Andreesenella_2 (parent Andreesenella), which is registered but not yet built; build its tree then ReLeaf",
  "assembly": "/path/to/genome.fna",
  "assembly_id": "genome_id",
  "query_taxonomy": "d__Bacteria;...;g__Andreesenella;s__",
  "matched_database": "Andreesenella_2",
  "parent_taxon": "Andreesenella",
  "subclade_id": 2,
  "subclade_name": "Andreesenella_2",
  "subclade_taxonomy": "d__Bacteria;...;g__Andreesenella",
  "members_file": "/path/to/databases/Andreesenella_2_db/subclade_members.txt",
  "sketch_file": "/path/to/databases/Andreesenella_2_db/subclade_sketch.msh",
  "source_genome_dir": "/path/to/raw/Andreesenella"
}
```
Note `subclade_taxonomy` is the **subclade's own** recorded taxonomy, not the
query's — the on-demand build (Phase 3c) uses this so a rebuild never
overwrites the subclade's taxonomy with whatever a single query happened to
carry.

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

> **Subclade partitioning (`--max-tree-genomes`, default 150).** When a taxon
> downloads more genomes than the ceiling, the raw set is MASH-partitioned into
> size-bounded subclades (`<Taxon>_1`, `<Taxon>_2`, …) *before* QC, and a tree is
> built only for the subclade(s) that contain a query. This reorders the stages
> below: the download is split into a **`--download-only`** phase (Stage 1a) and a
> per-subclade **`--qc-only`** phase (Stage 1c), with partitioning in between
> (Stage 1b). When the raw count is under the ceiling the flow collapses to the
> classic single-tree path (one QC pass on the whole set). See
> `readmes/README_ADVANCED_FEATURES.md` § *Subclade Partitioning* for full detail.

### Stage 1a: Download Genomes (raw, pre-QC)

**Script**: `gather_filter_asms.sh --download-only`

**Purpose**: Download candidate genomes from NCBI and stage them as raw FASTAs at
`{download_dir}/assemblies_all.TMP/*.fna`, stopping before CheckM2 QC. Query
genomes are staged here too so partitioning sees them as clustering leaves.

**Checkpoint**: `download_{taxon_name}.flag`

### Stage 1b: Cap tree size (only if raw count > `--max-tree-genomes`)

By default the wrapper **diverse-subsamples** the oversized taxon down to one tree
(`python_scripts/subsample_genomes.py`, `--subsample-size` genomes) and continues
through the single-tree path (Stage 1c → 3 → 4).

With **`--megatree`** the wrapper instead partitions and builds a full-coverage
merged tree. `python_scripts/subclade_partition.py` runs `mash triangle -k 17 -s
5000 -E` on the raw set, clusters with average-linkage (UPGMA), and recursively
splits so every subclade ≤ `--subclade-size`. It writes `partition_manifest.json`,
a per-subclade `.msh` sketch, and a `.members.txt` list under
`02_orthophyl_novel/partitions/{taxon}/`. By default **every** subclade is then
QC'd and built (Stage 1c → 3 → 4); with **`--megatree-lazy`**, a subclade holding
no query is instead only *registered* — `create_hierarchical_database.py
--register-only`, `built=false`, sketch/members/source-dir recorded, no tree —
and built later on demand when a query MASH-matches it (Phase 3c, below). A
backbone tree is built from `--backbone-reps` diverse reps per (built) subclade,
and `python_scripts/megatree_graft.py` grafts each subclade tree onto its
backbone reps into one merged tree (`03_results/trees/orthophyl/{taxon}_megatree.nwk`),
flagging high-support topology conflicts to `{taxon}_megatree_conflicts.json`.
The backbone database is written with `is_backbone=true` so the router can
place a query into it directly with `--placement backbone`.

**Checkpoint**: `partition_{taxon_name}` (resume re-reads the manifest, never
re-runs mash, so subclade numbering is stable); `megatree_backbone_{taxon}` and
`megatree_graft_{taxon}` for the backbone build and graft.

### Stage 1c: Quality-Control a Subclade

**Script**: `gather_filter_asms.sh --qc-only`

**Purpose**: Download and filter high-quality genomes from NCBI

The wrapper stages a subclade's raw members into a per-subclade
`assemblies_all.TMP/` and runs QC **only** on those genomes (CheckM2 runs here,
not at download time). For the unpartitioned case this is simply the whole raw
set. The QC steps below are identical whether run as the classic combined
download+QC or as this `--qc-only` phase.

**Steps**:

1. **Download from NCBI**
   ```bash
   utils/gather_filter_asms.sh \
       {taxon_name} \
       {output_dir} \
       {threads} \
       [--lowmem]                        # Low RAM mode (halves DIAMOND RAM)
       [--use-bbmap]                     # Skip CheckM2
       [--query-genomes q1.fna,q2.fna]   # QC the query genomes too (see Stage 2)
       [--must-keep GCF_x,GCF_y|file]    # These MUST pass QC or the run aborts
       [--keep-failing-query]            # Warn (don't abort) when a query fails QC
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

### Stage 2: Query Genomes Through QC

**Purpose**: QC the query/input assemblies *alongside* the downloads, rather than
appending them to the filtered set unchecked.

The wrapper passes the query FASTA paths to the gather script via `--query-genomes`.
They are staged into the download set (`assemblies_all.TMP/`), run through the same
CheckM2 stats + threshold filter, and — if they pass — land in `genomes_to_keep/`.

**Failure handling** (elegant + explicit):
- A **query genome that fails QC aborts the run** by default, reporting which genome
  failed and on which metric (e.g. `completeness=82.0 < MIN_completeness=95`). Pass
  `--keep-failing-query` to downgrade this to a loud warning and force the query in.
- Accessions listed in `--must-keep` (comma-separated, or a file with one per line)
  **always** abort the run if QC drops them — they are required downstream and this is
  not overridable.

Per-genome failure reasons are recorded in `qc_removal_reasons.txt` in the download dir.

**Result**: Complete, QC-verified genome set for phylogenetic analysis. The wrapper
keeps a fallback copy step for the `--skip-download` path (no QC runs there).

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

## Phase 3c: Subclade-Build Route (Lazy Subclades)

**Purpose**: Build a lazily-registered megatree subclade on demand, then
ReLeaf the waiting queries onto it.

**Only reachable when `--megatree-lazy` was used at partition time.** A
subclade with no query at partition time is registered (`built=false`) rather
than built; this phase runs when a *later* query's taxonomy/MASH match routes
it to that unbuilt subclade (the `OrthoPhyl_subclade_build` routing decision
from Phase 2).

**Process** (`_phase_subclade_build` / `_process_subclade_build` in
`orthophyl_pipeline_wrapper.py`), grouped by subclade so one build serves every
query that landed there:

1. **Locate the subclade's raw members** via its recorded `source_genome_dir`
   (written at registration time) — the tree is built from the subclade's OWN
   genomes, **not** the routing query.
2. **QC + run OrthoPhyl** on those raw members (`_build_subclade`,
   `force=True`), using the subclade's own `subclade_taxonomy` — not the
   query's taxonomy, so a rebuild never overwrites the subclade's recorded
   taxonomy with whatever a single query happened to carry.
3. **Promote the database entry**: `create_hierarchical_database.py --force`
   overwrites the `built=false` placeholder with a real built entry
   (`built=true`, real tree, `orthophyl_run` symlink).
4. **ReLeaf** the waiting query assemblies onto the freshly-built tree.

**Checkpoints**: registration and build use *separate* keys —
`register_{subclade_name}` (written by `_register_lazy_subclade` at partition
time) vs `qc_{subclade_name}` / `orthophyl_{subclade_name}` /
`database_{subclade_name}` (written during the on-demand build). Sharing one
key between registration and build would let `--resume` mistake a placeholder
for a finished build.

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
