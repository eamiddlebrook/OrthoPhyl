# Troubleshooting Guide - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---


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

---

[← Back to Main README](../README.md)
