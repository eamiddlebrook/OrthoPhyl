# Assembly Router - Intelligent Routing for OrthoPhyl/ReLeaf

Automatically route new genome assemblies to either **ReLeaf** (fast, incremental addition) or **OrthoPhyl** (comprehensive analysis with related genomes) based on taxonomic representation in your reference database.

## Overview

When you have a new genome assembly, this tool answers the question:
- **Is this taxonomy already in my phylogeny?** → Use ReLeaf (fast)
- **Is this a novel taxon?** → Use OrthoPhyl with expanded genome set (comprehensive)

```
┌─────────────────┐
│  New Assembly   │
│  + Taxonomy     │
└────────┬────────┘
         │
         ↓
   ┌─────────────┐
   │   Router    │
   │  (queries   │
   │  database)  │
   └──────┬──────┘
          │
    ┌─────┴─────┐
    ↓           ↓
┌────────┐  ┌────────────┐
│ ReLeaf │  │ OrthoPhyl  │
│ (fast) │  │ (download  │
│        │  │  related   │
│        │  │  genomes)  │
└────────┘  └────────────┘
```

## Quick Start

### 1. Create Taxonomy Database

First, create a reference database from an existing OrthoPhyl run:

```bash
# Prepare genome metadata file
cat > genome_metadata.tsv << EOF
genome_id	taxonomy	level
GCF_000005845.2	Escherichia coli K-12	strain
GCF_000008865.2	Escherichia coli	species
GCF_000482265.1	Escherichia coli O157:H7	strain
GCF_000006945.2	Salmonella enterica	species
GCF_000011885.1	Salmonella enterica Typhimurium	strain
EOF

# Create database
python setup_taxonomy_database.py \
    --orthophyl-dir /path/to/previous_orthophyl_run/ \
    --genome-metadata genome_metadata.tsv \
    --output-dir reference_db/
```

This creates:
```
reference_db/
├── taxonomy_map.tsv      # Taxon → genome mapping
├── phylogeny.nwk         # Reference tree
├── metadata.json         # Database info
├── orthophyl_run/        # Symlink to OrthoPhyl outputs
└── database_summary.txt  # Human-readable summary
```

### 2. Route Single Assembly

```bash
python assembly_router.py \
    --assembly new_genome.fna \
    --taxonomy "Escherichia coli O104:H4" \
    --database-dir reference_db/ \
    --output-dir routing_results/ \
    -t 24
```

**Output:**
```
[INFO] Routing assembly: new_genome
[INFO] Taxonomy: Escherichia coli O104:H4
[INFO] ✓ Taxonomy matched in database (species level)
[INFO]   Matched: Escherichia coli
[INFO]   Representative genomes: 3

ROUTING COMPLETE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
new_genome:
  Pipeline: ReLeaf
  Reason: Taxonomy represented in database at species level
  Command: ./ReLeaf.sh --store reference_db/orthophyl_run ...

[INFO] Use --execute to run the recommended pipeline
```

### 3. Route Multiple Assemblies (Batch Mode)

Create input table:
```bash
cat > assemblies_to_route.tsv << EOF
assembly_path	taxonomy	assembly_id
genome1.fna	Escherichia coli	ASM001
genome2.fna	Klebsiella pneumoniae	ASM002
genome3.fna	Novel_bacterium sp.	ASM003
EOF

# Route all assemblies
python assembly_router.py \
    --batch assemblies_to_route.tsv \
    --database-dir reference_db/ \
    --output-dir batch_routing/ \
    -t 24
```

## Output Files

For each assembly, the router creates:

```
routing_results/
├── routing_decision_ASM001.json      # Machine-readable decision
├── routing_summary_ASM001.txt        # Human-readable summary
├── releaf_input/                     # Prepared inputs for ReLeaf
│   └── ASM001.fna
├── orthophyl_input/                  # Prepared inputs for OrthoPhyl
│   ├── ASM003.fna
│   └── downloaded_genomes/           # Related genomes from NCBI
└── routing_20240115_143022.log       # Detailed log
```

### Example Decision File (JSON)

```json
{
  "pipeline": "ReLeaf",
  "reason": "Taxonomy represented in database at species level",
  "assembly": "genome1.fna",
  "assembly_id": "ASM001",
  "taxonomy": "Escherichia coli O157:H7",
  "matched_taxonomy": "Escherichia coli",
  "match_type": "species",
  "database_genomes": [
    "GCF_000005845.2",
    "GCF_000008865.2",
    "GCF_000482265.1"
  ],
  "command": "./ReLeaf.sh --store reference_db/orthophyl_run ..."
}
```

### Example Summary File (TXT)

```
══════════════════════════════════════════════════════════════════
ASSEMBLY ROUTING DECISION
══════════════════════════════════════════════════════════════════

Assembly ID: ASM001
Assembly: genome1.fna
Taxonomy: Escherichia coli O157:H7
Pipeline: ReLeaf
Reason: Taxonomy represented in database at species level

Matched taxonomy: Escherichia coli
Match type: species
Database genomes: 3

══════════════════════════════════════════════════════════════════
COMMAND TO RUN
══════════════════════════════════════════════════════════════════

./ReLeaf.sh --store reference_db/orthophyl_run --input_genomes \
    routing_results/releaf_input -t 24 --tree_method iqtree \
    --TREE_DATA CDS
```

## Decision Logic

The router uses hierarchical taxonomy matching:

### 1. Exact Match
```
Query: "Escherichia coli K-12"
Database has: "Escherichia coli K-12"
→ EXACT MATCH → ReLeaf
```

### 2. Species-Level Match
```
Query: "Escherichia coli O157:H7"  (strain)
Database has: "Escherichia coli"    (species)
→ SPECIES MATCH → ReLeaf
```

### 3. Genus-Level Match
```
Query: "Escherichia new_species"   (species)
Database has: "Escherichia"         (genus)
→ GENUS MATCH → ReLeaf (with warning)
```

### 4. No Match
```
Query: "Klebsiella pneumoniae"
Database has: "Escherichia", "Salmonella"
→ NO MATCH → OrthoPhyl
           → Downloads related Klebsiella genomes from NCBI
```

## Taxonomy Database Format

### taxonomy_map.tsv

```tsv
taxonomy	level	genomes
Escherichia	genus	GCF_000005845.2,GCF_000008865.2,GCF_000482265.1
Escherichia coli	species	GCF_000005845.2,GCF_000008865.2,GCF_000482265.1,GCF_000742135.1
Escherichia coli K-12	strain	GCF_000005845.2
Salmonella	genus	GCF_000006945.2,GCF_000011885.1
Salmonella enterica	species	GCF_000006945.2,GCF_000011885.1,GCF_000022165.1
```

**Columns:**
- `taxonomy`: Full taxonomic name
- `level`: `genus`, `species`, or `strain`
- `genomes`: Comma-separated genome IDs (matches files in OrthoPhyl run)

### metadata.json

```json
{
  "created": "2024-01-15T14:30:22",
  "version": "1.0",
  "description": "Taxonomy database for assembly routing",
  "orthophyl_dir": "/path/to/orthophyl_run",
  "n_genomes": 15,
  "n_taxa": 8,
  "taxonomic_levels": ["genus", "species", "strain"],
  "database_type": "OrthoPhyl_reference"
}
```

## Advanced Usage

### Custom NCBI Downloads

By default, the router downloads up to 50 related genomes for novel taxa:

```bash
python assembly_router.py \
    --assembly novel_genome.fna \
    --taxonomy "Novel_species candidatus" \
    --database-dir reference_db/ \
    --ncbi-datasets /path/to/datasets \
    --output-dir results/
```

The script will:
1. Search NCBI for related genomes at genus or species level
2. Download assemblies using `ncbi-datasets` CLI
3. Prepare OrthoPhyl command with expanded genome set

### Execute Pipelines Automatically

```bash
# Dry run (show commands only)
python assembly_router.py \
    --assembly genome.fna \
    --taxonomy "Escherichia coli" \
    --database-dir reference_db/ \
    --dry-run

# Execute the recommended pipeline
python assembly_router.py \
    --assembly genome.fna \
    --taxonomy "Escherichia coli" \
    --database-dir reference_db/ \
    --execute
```

### Integration with Upstream Pipelines

If your assemblies come from another pipeline with taxonomy classifications:

```bash
# Your pipeline outputs: assembly.fna + taxonomy.txt

# Extract taxonomy
taxonomy=$(cat taxonomy.txt)

# Route assembly
python assembly_router.py \
    --assembly assembly.fna \
    --taxonomy "$taxonomy" \
    --assembly-id $(basename assembly.fna .fna) \
    --database-dir reference_db/ \
    --output-dir routing/

# Parse decision and execute
decision=$(cat routing/routing_decision_*.json)
pipeline=$(echo $decision | jq -r '.pipeline')
command=$(echo $decision | jq -r '.command')

echo "Routing to: $pipeline"
eval $command
```

## Batch Processing with Parallel Execution

For large-scale routing:

```bash
# 1. Route all assemblies (fast)
python assembly_router.py \
    --batch assemblies.tsv \
    --database-dir reference_db/ \
    --output-dir batch_results/

# 2. Separate by pipeline
grep -l '"pipeline": "ReLeaf"' batch_results/*.json > releaf_list.txt
grep -l '"pipeline": "OrthoPhyl"' batch_results/*.json > orthophyl_list.txt

# 3. Run ReLeaf assemblies in parallel
cat releaf_list.txt | while read decision_file; do
    command=$(jq -r '.command' $decision_file)
    echo $command
done | parallel -j 4

# 4. Run OrthoPhyl groups
cat orthophyl_list.txt | while read decision_file; do
    command=$(jq -r '.command' $decision_file)
    echo $command
done | parallel -j 2  # OrthoPhyl is more resource-intensive
```

## Updating the Database

As you add more genomes to your collection:

```bash
# 1. Run OrthoPhyl with new genome set
./OrthoPhyl.sh -g all_genomes/ -r reference -o new_run/ -t 48

# 2. Update metadata file
cat >> genome_metadata.tsv << EOF
GCF_NEW001	New_species	species
GCF_NEW002	New_species str123	strain
EOF

# 3. Rebuild database
python setup_taxonomy_database.py \
    --orthophyl-dir new_run/ \
    --genome-metadata genome_metadata.tsv \
    --output-dir reference_db/
```

## Troubleshooting

### Issue: "No HMM profiles found"

**Solution:** Ensure your OrthoPhyl run completed successfully and contains:
- `OG_alignmentsToHMM/hmms_final/*.hmm` (if ANI subset was used)
- OR `annots_prots.fixed/OrthoFinder/Results_ortho/MultipleSequenceAlignments/`

### Issue: "NCBI download failed"

**Solution:** Install ncbi-datasets CLI:
```bash
conda install -c conda-forge ncbi-datasets-cli
# Or
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
chmod +x datasets
```

Alternatively, manually download genomes and place in `orthophyl_input/`.

### Issue: "Taxonomy not matching"

**Problem:** Your query is "E. coli" but database has "Escherichia coli"

**Solution:** Use full genus/species names in taxonomy strings. The router does NOT expand abbreviations.

### Issue: "ReLeaf fails with new assembly"

**Possible causes:**
1. Assembly quality issues (too fragmented)
2. Highly divergent from reference
3. Incorrect taxonomy assignment

**Solution:** Check assembly quality with CheckM/QUAST, verify taxonomy, or force OrthoPhyl run.

## Command-Line Reference

### assembly_router.py

```
Required Arguments:
  --assembly PATH          Single assembly FASTA file
  --taxonomy STRING        Taxonomic classification
  --database-dir PATH      Reference taxonomy database

Optional Arguments:
  --assembly-id STRING     Custom assembly identifier
  --output-dir PATH        Output directory (default: routing_output)
  -t, --threads INT        Number of threads (default: 8)
  --ncbi-datasets PATH     Path to ncbi-datasets tool
  --execute                Execute recommended pipeline
  --dry-run                Show decision without executing

Batch Mode:
  --batch PATH             TSV file with multiple assemblies

Utilities:
  --create-example-db DIR  Create example database structure
```

### setup_taxonomy_database.py

```
Required Arguments:
  --orthophyl-dir PATH        OrthoPhyl output directory
  --genome-metadata PATH      TSV: genome_id, taxonomy, level
  --output-dir PATH           Output database directory
```

## Best Practices

1. **Keep database updated**: Rebuild after significant genome additions

2. **Use consistent taxonomy**: Match NCBI or GTDB taxonomy format

3. **Quality control**: Run CheckM on assemblies before routing

4. **Batch processing**: Route all assemblies first, then execute in groups

5. **Monitor resources**: OrthoPhyl downloads can be large (GB per species)

6. **Version control**: Keep metadata files in version control

7. **Documentation**: Record which assemblies were added when

## Performance

| Operation | Time | Resources |
|-----------|------|----------|
| Route single assembly | < 1 second | Minimal |
| Route 100 assemblies | < 10 seconds | Minimal |
| ReLeaf addition | 5-30 min | Depends on genome set size |
| OrthoPhyl new analysis | 2-24 hours | High CPU/memory |
| NCBI download (50 genomes) | 10-60 min | Network dependent |

## Citation

If you use this routing system, please cite:

- **OrthoPhyl**: Middlebrook EA, Katani R, Fair JM. (2024) OrthoPhyl - Streamlining large scale, orthology-based phylogenomic studies of bacteria at broad evolutionary scales. *G3 Genes|Genomes|Genetics*, jkae119.

- Your **HGT detection pipeline** (when published)

## License

This tool is provided as part of the HGTool package.

## Support

For issues or questions:
1. Check this README
2. Review log files in output directory
3. Open an issue on GitHub
4. Contact: [your contact info]

---

**Happy routing! 🧬**