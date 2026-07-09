# Assembly Router Pipeline

A multi-database assembly routing system for OrthoPhyl/ReLeaf that automatically determines whether a genome assembly should be processed with ReLeaf (using an existing database) or requires a new OrthoPhyl run for novel taxa.

## Overview

The Assembly Router Pipeline consists of two main components:

1. **Database Builder** (`create_hierarchical_database_v2.py`) - Creates queryable databases from OrthoPhyl runs
2. **Assembly Router** (`assembly_router_multi.cmd_out2.py`) - Routes assemblies to appropriate pipeline based on taxonomy

### Key Features

- **Automatic taxonomy matching** - Queries multiple databases to find best match
- **Hierarchical database structure** - Supports databases at different taxonomic ranks
- **Smart tree compatibility detection** - Detects available tree methods (iqtree, fasttree, raxml) and data types (CDS, PROT)
- **Automated genome downloading** - Integration with `gather_filter_asms.sh` for NCBI downloads
- **Batch processing** - Route multiple assemblies at once
- **Incremental updates** - Add new databases without rebuilding existing ones

---

## Installation & Requirements

### Dependencies

```bash
# Python 3.7+
pip install biopython  # If needed for sequence handling

# External tools (for genome downloading)
conda create -n gather_genomes \
    -c bioconda -c conda-forge \
    checkm-genome bbmap entrez-direct ncbi-datasets-cli
```

### File Structure

```
OrthoPhyl/
├── assembly_router/
│   ├── assembly_router_multi.cmd_out2.py
│   ├── create_hierarchical_database_v2.py
│   └── README_assembly_router.md (this file)
├── utils/
│   └── gather_filter_asms.sh
├── OrthoPhyl.sh
└── ReLeaf.sh
```

---

## Part 1: Database Creation

### Step 1: Prepare OrthoPhyl Runs Inventory

Create a TSV file listing your completed OrthoPhyl runs:

**Format:** `clade_name<TAB>orthophyl_dir<TAB>clade_taxonomy`

**Example: `orthophyl_runs.tsv`**
```tsv
Gaiellales	/data/orthophyl/gaiellales_run	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales
Escherichia	/data/orthophyl/ecoli_run	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia
Rhizobiaceae	/data/orthophyl/rhizobiaceae_run	d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae
```

**Important Notes:**
- Use GTDB taxonomy format with rank prefixes: `d__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`
- Each OrthoPhyl directory should contain completed run outputs
- Databases can be at any taxonomic rank (order, family, genus, etc.)

### Step 2: Create Initial Database

```bash
python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir /path/to/database_directory/
```

**What this does:**
1. Validates each OrthoPhyl run (checks for HMMs, alignments, trees)
2. Detects available tree methods (iqtree, fasttree, raxml)
3. Detects available data types (CDS, PROT)
4. Creates a `*_db` directory for each entry
5. Generates a master index (`database_index.json`)

**Output Structure:**
```
database_directory/
├── Gaiellales_db/
│   ├── database_config.json          # Metadata + available tree info
│   ├── orthophyl_run/                # Symlink to OrthoPhyl output
│   ├── phylogeny.nwk                 # Species tree
│   ├── genome_list.txt               # List of genomes
│   └── README.txt                    # Human-readable info
├── Escherichia_db/
│   └── ...
├── database_index.json               # Master index
└── database_summary.txt              # Human-readable summary
```

### Step 3: Update Database with New Entries

Add new lines to your `orthophyl_runs.tsv`, then:

```bash
python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir /path/to/database_directory/ \
    --update
```

**Update mode:**
- Skips existing databases (based on clade name)
- Only processes new entries
- Updates master index with all databases

### Step 4: Force Rebuild (if needed)

```bash
python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir /path/to/database_directory/ \
    --force
```

**Force mode:**
- Deletes and rebuilds all databases
- Use when OrthoPhyl runs have been updated

---

## Part 2: Assembly Routing

### Basic Usage: Single Assembly

```bash
python assembly_router/assembly_router_multi.cmd_out2.py \
    --assembly genome.fna \
    --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \
    --database-dir /path/to/database_directory/ \
    --output-dir routing_results/
```

### With Automatic Genome Downloading

```bash
python assembly_router/assembly_router_multi.cmd_out2.py \
    --assembly genome.fna \
    --taxonomy "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae;g__Rhizobium;s__" \
    --database-dir /path/to/database_directory/ \
    --gather-filter-script utils/gather_filter_asms.sh \
    --output-dir routing_results/ \
    --threads 16
```

### Batch Mode: Multiple Assemblies

**1. Create batch input file:**

**Format:** `assembly_path<TAB>taxonomy<TAB>[optional_id]`

**Example: `assemblies_to_route.tsv`**
```tsv
/data/genomes/GCF_001234567.fna	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae	GCF_001234567
/data/genomes/GCF_009876543.fna	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia	GCF_009876543
```

**2. Run batch routing:**

```bash
python assembly_router/assembly_router_multi.cmd_out2.py \
    --batch assemblies_to_route.tsv \
    --database-dir /path/to/database_directory/ \
    --gather-filter-script utils/gather_filter_asms.sh \
    --output-dir routing_results/ \
    --threads 16
```

---

## Understanding the Routing Logic

### Decision Tree

```
Query Assembly with Taxonomy
    ↓
    ├─→ MATCH FOUND in database?
    │   ├─→ YES: Route to ReLeaf
    │   │   ├─→ Use most specific matching database
    │   │   ├─→ Detect available tree methods & data types
    │   │   └─→ Generate compatible ReLeaf command
    │   │
    │   └─→ NO: Route to OrthoPhyl
    │       ├─→ Identify taxonomic rank for download
    │       ├─→ Generate gather_filter_asms.sh command
    │       ├─→ Generate OrthoPhyl command
    │       └─→ Suggest database creation
```

### Matching Rules

1. **Hierarchical matching** - Matches at any taxonomic level
   - Query: `g__Escherichia;s__coli` matches database for `g__Escherichia`
   - Query: `f__Rhizobiaceae;g__Rhizobium` matches database for `f__Rhizobiaceae`

2. **Most specific match** - Uses deepest matching rank
   - If databases exist for both `o__Enterobacterales` and `g__Escherichia`
   - Query with `g__Escherichia` will use the genus-level database

3. **Novel taxa handling** - No match triggers OrthoPhyl route
   - Suggests appropriate taxonomic level for genome download
   - Falls back to higher ranks if species/genus is undefined

---

## Output Files

### For Each Assembly

The router creates three output files per assembly:

#### 1. JSON Decision (`routing_decision_[ID].json`)
Machine-readable routing decision with all metadata.

```json
{
  "pipeline": "ReLeaf",
  "assembly_id": "GCF_001234567",
  "query_taxonomy": "d__Bacteria;p__Actinomycetota;...",
  "matched_database": "Gaiellales",
  "matched_rank": "order",
  "tree_method": "iqtree",
  "tree_data": "CDS",
  "available_methods": ["iqtree", "fasttree"],
  "available_data_types": ["CDS", "PROT"],
  "command": "..."
}
```

#### 2. Human-Readable Summary (`routing_summary_[ID].txt`)

```
======================================================================
ASSEMBLY ROUTING DECISION (Multi-Database Query)
======================================================================

Assembly ID: GCF_001234567
Assembly: /data/genomes/GCF_001234567.fna
Query Taxonomy: d__Bacteria;p__Actinomycetota;c__Thermoleophilia;...

Pipeline: ReLeaf
Reason: Taxonomy matches Gaiellales at order level

Match Details:
  Database: Gaiellales
  Matched at: order level
  Database directory: /path/to/Gaiellales_db
  Contains: 150 genomes

ReLeaf Configuration:
  Tree method: iqtree
  Data type: CDS
  Available methods: iqtree, fasttree
  Available data types: CDS, PROT

======================================================================
COMMAND(S) TO RUN
======================================================================

./ReLeaf.sh \
    --store /path/to/Gaiellales_db/orthophyl_run \
    --input_genomes routing_results/releaf_input \
    -t 16 \
    --tree_method iqtree \
    --TREE_DATA CDS
```

#### 3. Copied Assembly File

The assembly is copied to the appropriate input directory:
- **ReLeaf route**: `routing_results/releaf_input/`
- **OrthoPhyl route**: `routing_results/orthophyl_new_genomes/`

### Batch Processing Output

Additional summary file: `batch_routing_summary.txt`

```
======================================================================
BATCH ROUTING SUMMARY (Multi-Database)
======================================================================

Available databases: 3
  - Gaiellales (order, 150 genomes)
  - Escherichia (genus, 2840 genomes)
  - Rhizobiaceae (family, 450 genomes)

Total assemblies: 25
  → ReLeaf (matched): 18
  → OrthoPhyl (novel): 7

ReLeaf Routing by Database:
----------------------------------------------------------------------

Gaiellales (8 assemblies):
  - GCF_001234567
  - GCF_002345678
  ...

Escherichia (10 assemblies):
  - GCF_000005845
  ...

OrthoPhyl Assemblies (novel taxa):
----------------------------------------------------------------------
  GCF_999999999              - genus: Novigenus
  GCF_888888888              - family: Novifamily
  ...
```

---

## OrthoPhyl Route: Handling Novel Taxa

When an assembly doesn't match any database, the router generates a 3-step workflow:

### Generated Workflow

```bash
# Step 1: Download and filter related genomes from NCBI
# This will download genomes, run CheckM QC, and filter by quality
utils/gather_filter_asms.sh \
    Rhizobium \
    routing_results/orthophyl_new_genomes \
    16

# Note: The script will create routing_results/orthophyl_new_genomes/genomes_to_keep/
# with high-quality, non-redundant genomes.
# Quality filters (adjust in script if needed):
#   - Completeness >= 95%
#   - Contamination <= 1.0%
#   - Duplication <= 2%
#   - Removes RefSeq/GenBank redundancy
#   - Filters by assembly stats (N50, GC content, length)

# Step 2: Run OrthoPhyl on expanded genome set
# Note: Your query genome (GCF_001234567.fna) should be in routing_results/orthophyl_new_genomes/genomes_to_keep
./OrthoPhyl.sh \
    -g routing_results/orthophyl_new_genomes/genomes_to_keep \
    -o routing_results/orthophyl_output \
    -t 16 \
    --tree_method iqtree \
    --TREE_DATA CDS \
    --use_partitions true

# Step 3: Create new database entry
# Add this line to your orthophyl_runs.tsv:
Rhizobium	routing_results/orthophyl_output	d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae;g__Rhizobium

# Then rebuild the database index:
python create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir /path/to/database_directory/ \
    --update
```

### Taxonomic Rank Selection for Download

The router intelligently selects the appropriate rank for genome download:

1. **Uses most specific defined rank** in the query taxonomy
2. **Falls back to higher ranks** if species/genus is undefined:
   - If species (`s__`) is empty → use genus (`g__`)
   - If genus is also empty → use family (`f__`)

**Example:**
```
Query: "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae;g__Rhizobium;s__"
                                                                                                                  ↑ empty
Download target: Rhizobium (genus level)
```

---

## Advanced Features

### Tree Method & Data Type Compatibility

The router automatically detects what tree files are available in each database:

**Database validation detects:**
- **Tree methods**: iqtree, fasttree, raxml (from filename)
- **Data types**: CDS, PROT (from filename)

**ReLeaf command generation:**
- Prefers: `iqtree` > `raxml` > `fasttree`
- Prefers: `CDS` > `PROT`
- Falls back to available options if preferred not present
- Warns user if non-default options are used

**Example:**
```bash
# Database has: fasttree + CDS, iqtree + PROT
# Router will generate:

./ReLeaf.sh \
    --store /path/to/database/orthophyl_run \
    --input_genomes routing_results/releaf_input \
    -t 16 \
    --tree_method iqtree \
    --TREE_DATA PROT

# Note: Using tree_method=iqtree, TREE_DATA=PROT
# Available in database: methods=['fasttree', 'iqtree'], data_types=['CDS', 'PROT']
```

### Database Configuration File

Each database contains a `database_config.json` with metadata:

```json
{
  "created": "2024-01-15T10:30:00",
  "version": "1.0",
  "database_type": "hierarchical_taxonomy",
  "clade_name": "Gaiellales",
  "clade_taxonomy": "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales",
  "clade_rank": "o",
  "clade_rank_name": "order",
  "orthophyl_source": "/data/orthophyl/gaiellales_run",
  "n_genomes": 150,
  "has_hmms": true,
  "has_trees": true,
  "available_tree_methods": ["iqtree", "fasttree"],
  "available_data_types": ["CDS", "PROT"],
  "available_trees": [
    {
      "filename": "SCO_strict.CDS.iqtree.treefile",
      "method": "iqtree",
      "data_type": "CDS"
    },
    {
      "filename": "SCO_strict.CDS.fasttree.tree",
      "method": "fasttree",
      "data_type": "CDS"
    },
    {
      "filename": "SCO_strict.PROT.iqtree.treefile",
      "method": "iqtree",
      "data_type": "PROT"
    }
  ],
  "validation_details": [
    "Found 450 HMM/alignment files in OG_alignmentsToHMM/hmms_final",
    "Found 3 tree(s): SCO_strict.CDS.iqtree.treefile, ..."
  ]
}
```

---

## Command-Line Reference

### create_hierarchical_database_v2.py

```bash
python create_hierarchical_database_v2.py [OPTIONS]

Required:
  --input FILE              TSV file: clade_name, orthophyl_dir, clade_taxonomy
  --output-dir DIR          Output directory for databases

Optional:
  --update                  Skip existing databases, add new ones only
  --force                   Rebuild all databases (overwrite existing)
```

### assembly_router_multi.cmd_out2.py

```bash
python assembly_router_multi.cmd_out2.py [OPTIONS]

Input (one required):
  --assembly FILE           Single assembly FASTA
  --batch FILE              TSV file: assembly_path, taxonomy, [id]

Required:
  --database-dir DIR        Directory containing *_db databases

For --assembly mode:
  --taxonomy STRING         GTDB taxonomy string (required with --assembly)
  --assembly-id STRING      Assembly identifier (optional)

Optional:
  --output-dir DIR          Output directory (default: routing_output)
  --gather-filter-script    Path to gather_filter_asms.sh
  -t, --threads INT         Number of threads (default: 8)
```

---

## Troubleshooting

### Database Creation Issues

**Problem: "Invalid OrthoPhyl run"**
```
Validation details:
  - WARNING: No HMM profiles found
  - WARNING: No species tree found
```

**Solution:** Ensure OrthoPhyl run completed successfully. Check for:
- `OG_alignmentsToHMM/hmms_final/` directory with `.hmm` files
- `FINAL_SPECIES_TREES/` or `phylo_current/SpeciesTree/` with tree files

**Problem: "Could not determine number of genomes"**

**Solution:** Ensure one of these files exists in the OrthoPhyl directory:
- `genome_list`
- `all_input_list`
- `annots_prots/*.faa` files

### Routing Issues

**Problem: No match found but expected a match**

**Solution:** Check taxonomy format:
- Must use GTDB format with rank prefixes: `d__`, `p__`, `c__`, etc.
- Remove spaces after semicolons: `d__Bacteria;p__Actino...` not `d__Bacteria; p__Actino...`
- Ensure all ranks from domain to query rank are present

**Problem: "gather_filter_asms.sh not found"**

**Solution:** Either:
1. Provide correct path: `--gather-filter-script utils/gather_filter_asms.sh`
2. Or omit the flag - router will generate manual download instructions

### ReLeaf Command Issues

**Problem: ReLeaf fails with "tree method not available"**

**Solution:** Check database config:
```bash
cat /path/to/database_db/database_config.json | grep available_tree_methods
```

The router should only generate commands for available methods. If mismatch occurs:
1. Rebuild database: `--force`
2. Check tree files in original OrthoPhyl run

---

## Best Practices

### Database Organization

1. **One database per taxonomic level**
   - Don't create overlapping databases (e.g., both `Escherichia` and `Enterobacteriaceae` with `Escherichia` genomes)
   - Router will use the most specific match

2. **Regular updates**
   - Add new OrthoPhyl runs with `--update` mode
   - No need to rebuild existing databases

3. **Naming conventions**
   - Use clear, consistent clade names in `orthophyl_runs.tsv`
   - Avoid spaces and special characters (will be converted to underscores)

### Genome Downloading

1. **Use gather_filter_asms.sh for reproducibility**
   - Automatic quality filtering (CheckM)
   - Removes RefSeq/GenBank redundancy
   - Filters by assembly statistics

2. **Quality thresholds** (in `gather_filter_asms.sh`):
   - Completeness ≥ 95%
   - Contamination ≤ 1.0%
   - Duplication ≤ 2%
   - Adjust as needed for your taxonomic group

### Batch Processing

1. **Group assemblies by expected outcome**
   - Process known taxa separately from novel taxa
   - Makes monitoring and troubleshooting easier

2. **Check batch summary before running commands**
   - Review `batch_routing_summary.txt`
   - Verify routing decisions make sense
   - Run a small test batch first

---

## Workflow Example: Complete Pipeline

### Scenario
You have:
- 3 completed OrthoPhyl runs (Gaiellales, Escherichia, Rhizobiaceae)
- 20 new genome assemblies to classify

### Step-by-Step

```bash
# 1. Create database inventory
cat > orthophyl_runs.tsv << 'EOF'
Gaiellales	/data/orthophyl/gaiellales	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales
Escherichia	/data/orthophyl/ecoli	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia
Rhizobiaceae	/data/orthophyl/rhizobiaceae	d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Rhizobiaceae
EOF

# 2. Build databases
python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir databases/ \
    2>&1 | tee database_creation.log

# 3. Review database summary
cat databases/database_summary.txt

# 4. Prepare batch input (assemblies with GTDB taxonomy)
# Get taxonomy from GTDB-Tk or other source
cat > assemblies_batch.tsv << 'EOF'
/data/new_genomes/GCF_001.fna	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae	GCF_001
/data/new_genomes/GCF_002.fna	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia	GCF_002
EOF

# 5. Route assemblies
python assembly_router/assembly_router_multi.cmd_out2.py \
    --batch assemblies_batch.tsv \
    --database-dir databases/ \
    --gather-filter-script utils/gather_filter_asms.sh \
    --output-dir routing_results/ \
    --threads 32 \
    2>&1 | tee routing.log

# 6. Review routing decisions
cat routing_results/batch_routing_summary.txt

# 7. Run ReLeaf for matched assemblies
# Commands are in routing_results/routing_summary_*.txt
for summary in routing_results/routing_summary_GCF_*.txt; do
    # Extract and run ReLeaf command
    awk '/^\.\/ReLeaf.sh/,/^$/' "$summary" | bash
done

# 8. For novel taxa, run the 3-step workflow
# (Download genomes, run OrthoPhyl, create new database)

# 9. Add new databases and update index
# After OrthoPhyl runs complete, add to orthophyl_runs.tsv
echo "Novigenus	routing_results/orthophyl_output	d__Bacteria;..." >> orthophyl_runs.tsv

python assembly_router/create_hierarchical_database_v2.py \
    --input orthophyl_runs.tsv \
    --output-dir databases/ \
    --update
```

---

## Version History

### v2 (Current - assembly_router_multi.cmd_out2.py)
- Multi-database query support
- Automatic tree method/data type detection
- Smart ReLeaf command generation with compatibility checking
- Integration with gather_filter_asms.sh
- Batch processing mode
- Incremental database updates

### v1 (Legacy)
- Single database routing
- Manual database specification
- Basic OrthoPhyl/ReLeaf decision logic

---

## Future Development (v3 - Planned)

Potential features for `assembly_router_multi.cmd_out3.py`:
- Confidence scoring for taxonomy matches
- Support for partial taxonomy strings
- Database quality metrics (genome count, alignment quality)
- Automatic GTDB-Tk integration
- Parallel batch processing
- ReLeaf result aggregation

---

## Citation

If you use this pipeline in your research, please cite:

**OrthoPhyl:**
```
Earl A Middlebrook, Robab Katani, Jeanne M Fair
OrthoPhyl - Streamlining large scale, orthology-based phylogenomic studies of bacteria at broad evolutionary scales
G3 Genes|Genomes|Genetics, 2024; jkae119
https://doi.org/10.1093/g3journal/jkae119
```

---

## Support

For issues, questions, or contributions:
- GitHub: https://github.com/eamiddlebrook/OrthoPhyl
- Open an issue with:
  - Command used
  - Error message
  - Relevant log files
  - Database summary output

---

## License

Same as OrthoPhyl main license.

---

**Last Updated:** January 2025  
**Version:** 2.0  
**Maintainer:** Generated for HGTool / OrthoPhyl project