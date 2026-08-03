# Advanced Features - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---


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

### 9. Taxon Mode - Automated Database Creation

**Purpose**: Automatically download genomes for a taxon and create a new database

**Create Mode** (new database from taxon):

```bash
python orthophyl_pipeline_wrapper.v2.py \
    --taxon "Methylorubrum" \
    --taxon-rank genus \
    --database-dir databases/ \
    --output-dir methylorubrum_run/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32
```

**Required Arguments for Taxon Create Mode**:
- `--taxon`: Taxon name (e.g., "Methylorubrum", "Escherichia coli")
- `--database-dir`: Where to create the new database
- `--output-dir`: Working directory for this run
- `--gather-script`: **REQUIRED** - Path to genome download script (e.g., `utils/gather_filter_asms.sh`)

**Optional Arguments**:
- `--taxon-rank`: Taxonomic rank (species, genus, family, etc.) - auto-detected if not provided
- `--threads`: Number of CPU threads
- `--low-ram`: Use CheckM reduced tree mode
- `--use-bbmap`: Skip CheckM, use bbmap for faster QC

**What Happens**:
1. Queries NCBI for all assemblies matching the taxon
2. Downloads genomes using `gather_filter_asms.sh`
3. Applies quality filters (CheckM: completeness ≥95%, contamination ≤1%)
4. Runs full OrthoPhyl pipeline on filtered genomes
5. Creates new database in `database_dir/{taxon}_db/`
6. Database is immediately available for future ReLeaf runs

**Update Mode** (add new assemblies to existing database):

```bash
python orthophyl_pipeline_wrapper.v2.py \
    --taxon "Methylorubrum" \
    --update-existing \
    --database-dir databases/ \
    --output-dir methylorubrum_update/ \
    --threads 32
```

**What Happens**:
1. Finds existing `Methylorubrum_db` in database directory
2. Queries NCBI for current assemblies
3. Identifies NEW assemblies not in database
4. Downloads only new assemblies
5. Runs ReLeaf to add them to existing phylogeny
6. Updates database metadata and version

**Example Workflow**:

```bash
# Year 1: Create initial database
python orthophyl_pipeline_wrapper.v2.py \
    --taxon "Escherichia" \
    --taxon-rank genus \
    --database-dir /data/databases/ \
    --output-dir /data/runs/escherichia_2025/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 64

# Result: /data/databases/Escherichia_db/ created with 500 genomes

# Year 2: Update with new assemblies
python orthophyl_pipeline_wrapper.v2.py \
    --taxon "Escherichia" \
    --update-existing \
    --database-dir /data/databases/ \
    --output-dir /data/runs/escherichia_update_2026/ \
    --threads 64

# Result: 50 new genomes added via ReLeaf, database updated to 550 genomes
```

**See Also**: `TAXON_MODE_GUIDE.md` for complete documentation

---

---

[← Back to Main README](../README.md)
