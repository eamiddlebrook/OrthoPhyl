# Troubleshooting Guide - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---


### Common Issues and Solutions

#### 1. "No databases found in {database_dir}"

**Cause**: Database directory is empty or improperly formatted

**Solutions**:
```bash
# Option A: Provide orthophyl_runs.tsv for initial creation
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --orthophyl-runs orthophyl_runs.tsv

# Option B: Manually create databases first
python assembly_router/OP_database_tool.py \
    --input orthophyl_runs.tsv \
    --output-dir databases/
```

---

#### 2. "Routing failed. Check log: routing.log"

**Cause**: assembly_router.py encountered an error

**Debug**:
```bash
# Check the log
cat output_dir/logs/routing.log

# Common issues:
# - Invalid taxonomy format
# - Missing assembly files
# - Corrupted database configs

# Test routing manually
python assembly_router/assembly_router.py \
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

python orthophyl_pipeline_wrapper.py \
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
    -s databases/{database}_db/current/orthophyl_run \
    -g test_genomes/ \
    -t 8 \
    -p iqtree \
    -o CDS
```

---

#### 6. "CheckM2 failed" or "Out of memory"

**Cause**: CheckM2 (`checkm2 predict`) uses DIAMOND + pretrained ML models, so it
needs far less RAM than the legacy CheckM1 (which ran pplacer against a large
reference tree). Memory pressure is usually driven by DIAMOND on large batches.

**Solutions**:
```bash
# Option A: Use low RAM mode (passes CheckM2 --lowmem, halving DIAMOND RAM)
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir databases/ \
    --output-dir results/ \
    --gather-script utils/gather_filter_asms.sh \
    --low-ram

# Option B: Skip CheckM2 entirely (fastest, lowest RAM, less stringent)
python orthophyl_pipeline_wrapper.py \
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
python assembly_router/OP_database_tool.py \
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

#### 9. Router loads 0 databases (worked fine before, moved host/container)

**Cause**: `database_index.json` and/or individual `database_config.json`
files store absolute paths from the environment the database was built in
(e.g. `/databases/...` from inside a container bind mount). On the host, or
under a different bind, those absolute paths don't exist, so every entry is
silently dropped. This affects any database built before path-portability was
added, or built with a bind-mount path that doesn't match where you're now
reading it from.

**Solutions**:
```bash
# Preview what would change (no writes)
python assembly_router/OP_database_tool.py \
    --migrate --output-dir databases/ --dry-run

# If the old mount prefix isn't inferred correctly, give it explicitly
python assembly_router/OP_database_tool.py \
    --migrate --output-dir databases/ --strip-prefix /databases --dry-run

# Run it for real once the dry-run output looks right
python assembly_router/OP_database_tool.py \
    --migrate --output-dir databases/
```

`--migrate` rewrites every stored path (index entries, `sketch_file`,
`members_file`, `source_genome_dir`, `orthophyl_source`, the `orthophyl_run`
symlink target, `orthophyl_runs.tsv`) relative to each database's own
directory, and moves megatree subclade/backbone databases into the nested
`<Taxon>/<Taxon>_N_db` layout if they aren't already there. It reports (but
does not delete) any dangling link it finds — see the next entry. It's safe
to re-run; a second pass over an already-migrated tree is a no-op.

---

#### 10. Dangling `orthophyl_run` symlink inside a database

**Cause**: a database's `orthophyl_run` link (or `current/orthophyl_run`)
points into a `.pipeline_runs/<run>/` directory that has since been deleted.
OrthoPhyl databases don't copy their HMMs/alignments/trees into the database
directory — they're read through this symlink — so deleting an old pipeline
run directory breaks every database that still points at it.

**Debug**:
```bash
# Find dangling links under a database tree
find databases/ -xtype l -name orthophyl_run

# OP_database_tool.py --migrate also reports these (without deleting them)
python assembly_router/OP_database_tool.py --migrate --output-dir databases/ --dry-run
```

**Solutions**:
- If the original `.pipeline_runs/<run>/` directory still exists elsewhere, re-point the symlink at it.
- Otherwise the database is unrecoverable; remove it:
```bash
python assembly_router/OP_database_tool.py --remove {database_name} --output-dir databases/
```
- Going forward, don't delete `.pipeline_runs/<run>/` directories that still have databases pointing at them.

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
