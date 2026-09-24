# Advanced Features - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---


### 1. Dry Run Mode

**Purpose**: Preview what would happen without executing

**Usage**:
```bash
python orthophyl_pipeline_wrapper.py \
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
python orthophyl_pipeline_wrapper.py --input assemblies.tsv ...
```
- Minimal console output
- All details in log files
- Best for production runs

**Level 1** (`-v`):
```bash
python orthophyl_pipeline_wrapper.py --input assemblies.tsv ... -v
```
- Shows stdout from subprocesses
- stderr still goes to log files
- Good for monitoring progress

**Level 2** (`-vv`):
```bash
python orthophyl_pipeline_wrapper.py --input assemblies.tsv ... -vv
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
   python orthophyl_pipeline_wrapper.py \
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
python orthophyl_pipeline_wrapper.py --input assemblies.tsv ...

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
python orthophyl_pipeline_wrapper.py \
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
python orthophyl_pipeline_wrapper.py \
    --input batch_aa \
    --database-dir databases/ \
    --output-dir results_batch1/ \
    --threads 16 &

python orthophyl_pipeline_wrapper.py \
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
python orthophyl_pipeline_wrapper.py \
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

python orthophyl_pipeline_wrapper.py \
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
python orthophyl_pipeline_wrapper.py \
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
- `--gather-script`: **REQUIRED** - Path to genome download script (e.g., `utils/gather_filter_asms.sh`)

**Optional Arguments**:
- `--output-dir`: Working directory for this run (default: `<database-dir>/.pipeline_runs/<taxon>_<timestamp>`, where `<taxon>` is the `--taxon` name or `TaxID<num>` for a numeric TaxID)
- `--taxon-rank`: Taxonomic rank (species, genus, family, etc.) - auto-detected if not provided
- `--threads`: Number of CPU threads
- `--low-ram`: Pass CheckM2 `--lowmem` (halves DIAMOND RAM)
- `--use-bbmap`: Skip CheckM2, use bbmap for faster QC
- `--must-keep`: Accessions that MUST survive QC or the run aborts with a per-metric
  report. Supply either a comma-separated list (`GCF_000...,GCF_001...`) or a path to a
  file with one accession per line.
- `--keep-failing-query`: Let query/input genomes that fail QC through with a loud
  warning instead of aborting (default: a query genome failing QC aborts the run).

**What Happens**:
1. Queries NCBI for all assemblies matching the taxon
2. Downloads genomes using `gather_filter_asms.sh`
3. Applies quality filters (CheckM2: completeness ≥95%, contamination ≤1%)
4. Runs full OrthoPhyl pipeline on filtered genomes
5. Creates new database in `database_dir/{taxon}_db/`
6. Database is immediately available for future ReLeaf runs

**Update Mode** (add new assemblies to existing database):

```bash
python orthophyl_pipeline_wrapper.py \
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
python orthophyl_pipeline_wrapper.py \
    --taxon "Escherichia" \
    --taxon-rank genus \
    --database-dir /data/databases/ \
    --output-dir /data/runs/escherichia_2025/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 64

# Result: /data/databases/Escherichia_db/ created with 500 genomes

# Year 2: Update with new assemblies
python orthophyl_pipeline_wrapper.py \
    --taxon "Escherichia" \
    --update-existing \
    --database-dir /data/databases/ \
    --output-dir /data/runs/escherichia_update_2026/ \
    --threads 64

# Result: 50 new genomes added via ReLeaf, database updated to 550 genomes
```

**See Also**: `TAXON_MODE_GUIDE.md` for complete documentation

---

### 10. Large-Taxon Handling (`--max-tree-genomes`, `--subsample-size`)

Some taxa (large genera especially) have far more assemblies than OrthoPhyl can
put into a single tree in reasonable time/RAM. The `ANI_shortlist` (`-n`) only
shrinks the OrthoFinder step — the full genome set still lands in the final
tree. To cap tree size, the pipeline wrapper **diverse-subsamples an oversized
taxon down to one manageable tree** by default.

#### Default: diverse subsampling

When a taxon downloads more raw genomes than `--max-tree-genomes` (default
**2000**), the wrapper builds **one** tree from a maximally-diverse MASH subsample
of `--subsample-size` genomes (default **500**):

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir /data/databases/ \
    --output-dir /data/runs/ \
    --gather-script utils/gather_filter_asms.sh \
    --max-tree-genomes 2000 \
    --subsample-size 500 \
    --threads 64
```

- **Method** (`python_scripts/subsample_genomes.py`): sketch every genome once
  with `mash sketch` (linear — no all-vs-all matrix), then greedily pick the most
  diverse subset via farthest-point (max-min) sampling. Each pick is a single
  `mash dist <candidate> combined.msh` call, so the whole selection is O(n·N) time
  and O(n) memory. This scales to very large taxa (tens of thousands of genomes)
  that would OOM the O(n²) partitioner.
- **Query / must-keep genomes are always retained**: they seed the greedy pick, so
  the subset is guaranteed to include them.
- Selection is deterministic (sorted input, lexical tie-break), so `--resume`
  re-uses the same subset.
- QC (CheckM2) runs on the subsampled set, so it also skips genomes that will
  never enter the tree.

The per-subclade **partitioning / megatree** path below is opt-in (enabled with
`--megatree`) and is what `--max-total-genomes` guards.

#### Opt-in: full-coverage megatree (`--megatree`)

Where diverse subsampling keeps *one* tree by discarding most genomes, `--megatree`
keeps **every** genome: it partitions the oversized taxon into size-bounded
subclades, builds a full tree per subclade, builds a small backbone tree from a few
diverse representatives of each subclade, and grafts each subclade's full tree onto
its representatives in the backbone — producing **one merged tree containing every
genome**.

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir /data/databases/ \
    --output-dir /data/runs/ \
    --gather-script utils/gather_filter_asms.sh \
    --megatree \
    --max-tree-genomes 2000 \
    --subclade-size 150 \
    --backbone-reps 5 \
    --conflict-min-support 90 \
    --max-total-genomes 5000 \
    --threads 64
```

`--megatree` is only engaged when a taxon's raw count exceeds `--max-tree-genomes`;
under that ceiling it collapses to the normal single-tree build. It is **not**
mutually exclusive with subsampling — subsampling is simply the default when
`--megatree` is absent.

**How it works** (novel-taxon / OrthoPhyl route):

1. **Download raw, pre-QC.** The wrapper downloads the full candidate genome set
   with `gather_filter_asms.sh --download-only` — it stops *before* the expensive
   CheckM2 QC pass.
2. **Enforce the guardrail.** The partitioner builds a dense `N×N` MASH matrix, so
   `--max-total-genomes` (default 5000) is enforced here: a raw set larger than the
   ceiling is refused rather than OOM-killing the node.
3. **Partition.** `python_scripts/subclade_partition.py` runs `mash triangle`
   (all-vs-all, the same `-k 17 -s 5000` parameters OrthoPhyl uses), clusters with
   average-linkage (UPGMA), and recursively splits the tree so every subclade holds
   ≤ `--subclade-size` genomes (default **150**). Subclades are numbered
   deterministically: `Andreesenella_1`, `Andreesenella_2`, … A combined MASH sketch
   (`.msh`) and a member list are written per subclade.
4. **Build every subclade** (or lazily register it — see below). For each subclade
   the wrapper QCs its raw members (`gather_filter_asms.sh --qc-only`, CheckM2 runs
   here) and runs OrthoPhyl, producing a full per-subclade tree.
5. **Build the backbone.** From each subclade, `min(subclade_size, --backbone-reps)`
   (default **5**) diverse representatives are picked (that subclade's own MASH
   greedy max-min, seeded by any query genomes so they anchor the backbone). The
   pooled reps are run through OrthoPhyl once to produce a backbone tree.
6. **Graft.** `python_scripts/megatree_graft.py` replaces each subclade's
   representative clade in the backbone with that subclade's full tree
   (MRCA-replace with a monophyly check), writing the merged tree to
   `03_results/trees/orthophyl/<taxon>_megatree.nwk`.
7. **Flag conflicts (do not resolve).** Where a subclade tree and the backbone
   disagree on a bipartition that is strongly supported (≥ `--conflict-min-support`,
   default **90**) on both sides, the disagreement is recorded to
   `<taxon>_megatree_conflicts.json`. Topology reconciliation is deliberate future
   work; this pass only flags.

The taxon database is created from the **backbone** OrthoPhyl run so ReLeaf has a
coherent HMM set. **`--taxon` create mode** takes the same path (there is no query,
so all subclades and their reps are built).

**Notes and caveats:**

- `--subclade-size` is a ceiling compared against the *raw* (pre-QC) count, so a
  subclade will usually end up somewhat smaller after QC.
- Non-monophyletic representatives (reps interleaved with other subclades in the
  backbone) are grafted best-effort — foreign leaves are preserved and the subclade
  is flagged `monophyletic: false` in the conflict report.
- Partitioning is deterministic (sorted input + UPGMA + size-desc numbering), so
  subclade names are stable across runs — required for `--resume`.
- **`--max-total-genomes` (default 5000)** guards this path only: it bounds the
  dense `N×N` MASH distance matrix (O(n²) in memory, ~20 GB at n=50k). The
  **default subsample path does not build this matrix and is unaffected** — it
  handles arbitrarily large taxa.

#### Opt-in: lazy build-on-demand (`--megatree-lazy`)

By default `--megatree` builds a full tree for **every** subclade immediately,
even ones no query landed in. `--megatree-lazy` defers that: a subclade with no
query at partition time is only *registered* — `built=false`, its MASH sketch +
member list + source genome directory recorded — instead of QC'd and built. This
lets an oversized taxon be partitioned once and have its subclades built
incrementally, only as queries actually route to them.

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir /data/databases/ \
    --output-dir /data/runs/ \
    --gather-script utils/gather_filter_asms.sh \
    --megatree --megatree-lazy \
    --max-tree-genomes 2000 \
    --subclade-size 150 \
    --threads 64
```

A later run whose query taxonomy-matches the registered (but unbuilt) subclade
triggers an on-demand build: the wrapper QCs and runs OrthoPhyl on *that
subclade's own raw members* (not the query), promotes the database entry from
`built=false` to `built=true`, then ReLeafs the waiting query assemblies onto
the freshly-built tree. This is the `OrthoPhyl_subclade_build` routing decision
(Phase 3C) — see the pipeline-phases diagram in the main README.

#### Placing a query: dense subclade or sparse backbone (`--placement`)

A megatree's backbone and every one of its subclades are written with the
**same parent-level taxonomy string** (they differ only by name/rank
metadata), so a query can taxonomy-match several of them at once. `--placement`
picks how that tie is broken:

- **`subclade` (default)** — routes to the most MASH-similar dense subclade,
  for the best local phylogenetic resolution. The query is sketched
  (`mash sketch -k 17 -s 5000`, matching `subclade_partition.py`'s params) and
  compared against each tied subclade's sketch; ties are broken deterministically
  by ascending subclade name. If no subclade sketch is usable, this falls back
  to the backbone with a warning rather than failing to route.
- **`backbone`** — routes straight to the megatree's sparse overview tree, no
  MASH call needed. Useful when you want a quick broad-context placement rather
  than committing to one dense subclade.

```bash
python orthophyl_pipeline_wrapper.py \
    --input assemblies.tsv \
    --database-dir /data/databases/ \
    --output-dir /data/runs/ \
    --placement backbone \
    --threads 64
```

The backbone database is distinguished from a dense subclade sharing the same
taxonomy by an `is_backbone` flag in its `database_config.json` (written via
`create_hierarchical_database.py --is-backbone`).

**Note:** `MASH_K`/`MASH_S` in `assembly_router.py` must stay in lockstep with
the sketch parameters in `subclade_partition.py` and `script_lib/functions.sh`
— comparing sketches built with different `-k`/`-s` produces meaningless
distances.

---

---

[← Back to Main README](../README.md)
