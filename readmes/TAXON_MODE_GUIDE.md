# OrthoPhyl Taxon Mode - Complete Guide

## Overview

The **Taxon Mode** feature enables automatic assembly gathering and database management based on taxonomic queries. Instead of manually curating assembly lists, you can now:

1. **Create databases from taxon names** - Automatically download all assemblies for a taxon
2. **Update existing databases** - Check NCBI for new assemblies and add them via ReLeaf
3. **Track database provenance** - Know exactly which taxon and assemblies are in each database

## Quick Start

### Create New Database from Taxon

```bash
python orthophyl_pipeline_wrapper.py \
    --taxon "Methylorubrum" \
    --taxon-rank genus \
    --database-dir databases/ \
    --output-dir methylorubrum_run/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32
```

**Note**: The `--gather-script` argument is **required** for taxon create mode. It provides the robust genome download and QC filtering pipeline.

This will:
1. Query NCBI for all *Methylorubrum* assemblies
2. Download assemblies via `gather_filter_asms.sh` (with CheckM2/bbmap QC filtering)
3. Run OrthoPhyl to build phylogeny
4. Create `Methylorubrum_db` with full metadata

### Update Existing Database

```bash
python orthophyl_pipeline_wrapper.py \
    --taxon "Methylorubrum" \
    --update-existing \
    --database-dir databases/ \
    --output-dir methylorubrum_update/ \
    --threads 32
```

This will:
1. Query NCBI for all *Methylorubrum* assemblies
2. Compare with existing database
3. Download only NEW assemblies
4. Run ReLeaf to add them to existing phylogeny
5. Update database metadata

## Architecture

### Components

1. **`utils/taxon_assembly_gatherer.py`** - NCBI query and assembly management
2. **`assembly_router/create_hierarchical_database.py`** - Enhanced database creation with metadata
3. **`orthophyl_pipeline_wrapper.py`** - Integrated workflow orchestration

### Database Metadata Structure

Enhanced `database_config.json` now includes:

```json
{
  "created": "2025-07-21T10:00:00",
  "last_updated": "2025-07-21T11:00:00",
  "clade_name": "Methylorubrum",
  "clade_taxonomy": "d__Bacteria;p__Pseudomonadota;...",
  "n_genomes": 150,
  
  // NEW: Taxon source tracking
  "source_taxon_name": "Methylorubrum",
  "source_taxid": "34007",
  "source_rank": "genus",
  "assembly_accessions": [
    "GCF_000022085.1",
    "GCF_000317555.1",
    ...
  ],
  "n_assemblies_at_creation": 150,
  "quality_filters": {
    "min_completeness": 95.0,
    "max_contamination": 1.0,
    "min_n50": 5000
  }
}
```

## Detailed Usage

### Command-Line Arguments

#### Taxon Mode Arguments

- `--taxon NAME` - Taxon name to query (e.g., "Methylorubrum", "Escherichia coli")
- `--taxon-rank RANK` - Taxonomic rank: species, genus, family, order, class, phylum (optional, auto-detected)
- `--update-existing` - Update mode: add new assemblies to existing database

#### Standard Arguments

- `--database-dir DIR` - Directory containing `*_db` databases
- `--output-dir DIR` - Output directory for this run (default: `<database-dir>/.pipeline_runs/<taxon>_<timestamp>`, where `<taxon>` is the `--taxon` name or `TaxID<num>` for a numeric TaxID)
- `--threads N` - Number of CPU threads (default: 8)
- `--dry-run` - Preview what would happen without executing
- `--resume` - Resume from last checkpoint
- `-v, -vv` - Verbose output (1 or 2 levels)

### Workflow Details

#### Create Mode Workflow

```
1. Query NCBI
   ├─ Search for taxon by name
   ├─ Retrieve taxonomy information
   └─ Get list of all assemblies

2. Download Assemblies
   ├─ Download genome files
   ├─ Calculate statistics (CheckM2 or bbmap)
   └─ Filter by quality thresholds

3. Run OrthoPhyl
   ├─ Annotate proteins
   ├─ Run OrthoFinder
   ├─ Build alignments
   ├─ Infer phylogeny (IQ-TREE)
   └─ Generate HMM profiles

4. Create Database
   ├─ Validate OrthoPhyl outputs
   ├─ Create database directory structure
   ├─ Populate metadata with taxon info
   └─ Add to database index
```

#### Update Mode Workflow

```
1. Check Existing Database
   ├─ Find database matching taxon name
   ├─ Load current assembly list
   └─ Verify database is valid

2. Query NCBI for Updates
   ├─ Get all current assemblies for taxon
   ├─ Compare with existing assembly list
   └─ Identify NEW assemblies only

3. Download New Assemblies
   ├─ Download only new genomes
   └─ Apply same quality filters

4. Run ReLeaf
   ├─ Add new assemblies to existing database
   ├─ Update phylogeny incrementally
   └─ Preserve existing HMM profiles

5. Update Metadata
   ├─ Add new assembly accessions to list
   ├─ Update last_updated timestamp
   └─ Update n_genomes count
```

## TaxonAssemblyGatherer API

The core module for NCBI interaction:

```python
from utils.taxon_assembly_gatherer import TaxonAssemblyGatherer

# Initialize
gatherer = TaxonAssemblyGatherer(
    taxon_name="Methylorubrum",
    rank="genus",  # optional
    output_dir=Path("output/")
)

# Query NCBI
assemblies = gatherer.query_ncbi()
# Returns: List[Dict] with assembly metadata

# Filter by quality
filtered = gatherer.filter_assemblies(
    assemblies,
    min_completeness=95.0,
    max_contamination=1.0,
    min_n50=5000
)

# Compare with existing database
existing_accessions = ["GCF_000022085.1", ...]
new_assemblies = gatherer.compare_with_existing(
    assemblies,
    existing_accessions
)

# Download assemblies
gatherer.download_assemblies(
    assemblies,
    output_dir=Path("genomes/")
)

# Get taxonomy string
taxonomy = gatherer.get_taxonomy_string()
# Returns: "d__Bacteria;p__Pseudomonadota;..."
```

## Examples

### Example 1: Create Database for a Genus

```bash
# Create Methylorubrum database
python orthophyl_pipeline_wrapper.py \
    --taxon "Methylorubrum" \
    --taxon-rank genus \
    --database-dir /data/databases/ \
    --output-dir /data/runs/methylorubrum_2025/ \
    --threads 64 \
    -v
```

**Output:**
- `/data/databases/Methylorubrum_db/` - New database
- `/data/runs/methylorubrum_2025/orthophyl_run/` - OrthoPhyl outputs
- `/data/runs/methylorubrum_2025/downloaded_assemblies/` - Downloaded genomes

### Example 2: Update Database with New Assemblies

```bash
# Check for new Methylorubrum assemblies (6 months later)
python orthophyl_pipeline_wrapper.py \
    --taxon "Methylorubrum" \
    --update-existing \
    --database-dir /data/databases/ \
    --output-dir /data/runs/methylorubrum_update_2025/ \
    --threads 64
```

**Output:**
- Updates `/data/databases/Methylorubrum_db/database_config.json`
- `/data/runs/methylorubrum_update_2025/new_assemblies/` - New genomes only
- `/data/runs/methylorubrum_update_2025/releaf_update/` - ReLeaf outputs

### Example 3: Dry Run Preview

```bash
# Preview what would happen
python orthophyl_pipeline_wrapper.py \
    --taxon "Escherichia" \
    --taxon-rank genus \
    --database-dir databases/ \
    --output-dir test_run/ \
    --dry-run
```

**Output:**
- Shows what would be queried and downloaded
- No actual downloads or computations
- Useful for planning and validation

### Example 4: Species-Level Database

```bash
# Create database for specific species
python orthophyl_pipeline_wrapper.py \
    --taxon "Escherichia coli" \
    --taxon-rank species \
    --database-dir databases/ \
    --output-dir ecoli_run/ \
    --threads 32
```

## Quality Filtering

The `TaxonAssemblyGatherer` supports quality filtering:

```python
# In your code or via wrapper
filtered = gatherer.filter_assemblies(
    assemblies,
    min_completeness=95.0,    # Minimum % genome completeness
    max_contamination=1.0,     # Maximum % contamination
    min_n50=5000              # Minimum contig N50
)
```

**Default behavior:**
- No filtering applied by default
- All assemblies from NCBI are included
- Can be customized in future versions

## Database Metadata Tracking

### Why Track Metadata?

1. **Reproducibility** - Know exactly which assemblies were used
2. **Updates** - Identify new assemblies since database creation
3. **Quality Control** - Track filtering criteria applied
4. **Provenance** - Link databases back to source taxa

### Metadata Fields

| Field | Description | Example |
|-------|-------------|---------|
| `source_taxon_name` | Taxon used for query | "Methylorubrum" |
| `source_taxid` | NCBI Taxonomy ID | "34007" |
| `source_rank` | Taxonomic rank | "genus" |
| `assembly_accessions` | List of assembly IDs | ["GCF_000022085.1", ...] |
| `n_assemblies_at_creation` | Initial count | 150 |
| `last_updated` | Last update timestamp | "2025-07-21T11:00:00" |
| `quality_filters` | Filters applied | {"min_completeness": 95.0} |

## Troubleshooting

### No Assemblies Found

```
✗ No assemblies found for taxon 'MyTaxon'
```

**Solutions:**
1. Check taxon name spelling
2. Try different rank (genus vs species)
3. Verify taxon exists in NCBI Taxonomy
4. Check NCBI Assembly database has genomes for this taxon

### Database Already Exists

```
✗ Database already exists for taxon 'Methylorubrum'
  Use --update-existing to add new assemblies
```

**Solutions:**
1. Use `--update-existing` to add new assemblies
2. Use different `--output-dir` for a separate run
3. Manually remove old database if rebuilding

### Import Error

```
Failed to import TaxonAssemblyGatherer
```

**Solutions:**
1. Ensure `utils/taxon_assembly_gatherer.py` exists
2. Check Python path includes OrthoPhyl directory
3. Verify all dependencies installed (requests, etc.)

## Integration with Existing Workflows

### Batch Mode (Original)

Still fully supported:

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --threads 32
```

### Taxon Mode (New)

Alternative workflow:

```bash
python orthophyl_pipeline_wrapper.py \
    --taxon "Methylorubrum" \
    --database-dir databases/ \
    --output-dir results/ \
    --threads 32
```

**Note:** `--taxon` and `--input` are mutually exclusive.

## Future Enhancements

Potential improvements for future versions:

1. **Quality filter arguments** - Command-line control of filtering thresholds
2. **Multiple taxon support** - Process multiple taxa in one run
3. **Scheduled updates** - Cron-friendly update checking
4. **Assembly metadata** - Store more assembly information (date, submitter, etc.)
5. **Alternative sources** - Support for non-NCBI assembly sources
6. **Differential updates** - Smart detection of which assemblies need reprocessing

## Testing Recommendations

### Test 1: Small Taxon (Quick Test)

```bash
# Use a small genus for testing (~10-20 assemblies)
python orthophyl_pipeline_wrapper.py \
    --taxon "Methylorubrum" \
    --database-dir test_db/ \
    --output-dir test_run/ \
    --threads 8 \
    --dry-run  # Preview first
```

### Test 2: Update Mode

```bash
# After Test 1, simulate update
python orthophyl_pipeline_wrapper.py \
    --taxon "Methylorubrum" \
    --update-existing \
    --database-dir test_db/ \
    --output-dir test_update/ \
    --threads 8
```

### Test 3: Verify Metadata

```bash
# Check database metadata
cat test_db/Methylorubrum_db/database_config.json | jq .

# Verify fields present:
# - source_taxon_name
# - source_taxid
# - assembly_accessions
# - last_updated
```

## Summary

The Taxon Mode feature provides:

✅ **Automated assembly gathering** from NCBI by taxon name  
✅ **Database creation** with full OrthoPhyl pipeline  
✅ **Incremental updates** via ReLeaf for new assemblies  
✅ **Complete metadata tracking** for reproducibility  
✅ **Backward compatibility** with existing batch mode  

This enables scalable, reproducible phylogenomic database management for any taxonomic group in NCBI.

---

**Created:** July 21, 2026  
**Version:** 1.0  
**Author:** OrthoPhyl Development Team
