# Apptainer / Singularity Containers - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---

### Container Overview

OrthoPhyl is available as an Apptainer/Singularity container that includes all dependencies pre-installed. This is ideal for:
- HPC environments without conda/mamba
- Reproducible analyses
- Systems where you lack admin privileges
- Ensuring consistent software versions

`apptainer` and `singularity` are drop-in replacements for each other in every command on this page -- use whichever is installed on your system.

### Getting the Container

Pull the pre-built container from Sylabs Cloud Library:

```bash
# Create directory for containers (optional but recommended)
mkdir -p ~/singularity_images
cd ~/singularity_images

# Pull the latest OrthoPhyl container
apptainer pull library://earlyevol/default/orthophyl

# Rename to version-specific name
mv orthophyl_latest.sif OrthoPhyl.v3.3.sif
```

**Note**: The container is ~2.7 GB and includes all dependencies. No building or admin privileges required!

### Building the Container Yourself

The recipe is `Singularity.OP.v3.3.recipe` at the repo root. It stages the CheckM2 DIAMOND reference database from a pre-downloaded host copy (to avoid a ~2.9 GB network fetch on every build) via a **required** build-arg:

```bash
apptainer build \
    --build-arg checkm2_db_path=/path/to/CheckM2_database/uniref100.KO.1.dmnd \
    OrthoPhyl.v3.3.sif \
    Singularity.OP.v3.3.recipe
```

### Container Structure

The container includes:
- **Base environment**: OrthoPhyl, OrthoFinder, IQ-TREE, MAFFT, trimAl, FastTree, HMMER, etc.
- **gather_genomes environment**: CheckM2, bbmap, NCBI datasets CLI, entrez-direct (the CheckM2 DIAMOND database is baked in and `CHECKM2DB` is exported)
- **Additional tools**: ASTRAL, catfasta2phyml, Alignment Assessment
- **OrthoPhyl code**: Cloned from GitHub at `/opt/gits/OrthoPhyl/`

### `exec` vs `run`

The container's default `%runscript` just drops you into an interactive `/bin/bash` shell with every entrypoint already on `PATH` -- it does **not** forward arguments to a pipeline script. That means:

```bash
# This does NOT run OrthoPhyl.sh with these args -- it opens a shell and
# prints a banner, ignoring everything after the image name:
apptainer run OrthoPhyl.v3.3.sif -T TESTER -s output_dir -t 4   # WRONG

# Always exec the entrypoint directly instead:
apptainer exec OrthoPhyl.v3.3.sif \
    bash /opt/gits/OrthoPhyl/OrthoPhyl.sh -T TESTER -s output_dir -t 4
```

`apptainer run` (no args) is still useful on its own, to drop into the shell:

```bash
apptainer run OrthoPhyl.v3.3.sif        # drops into /bin/bash
```

### Critical: Database Directory Mounting

**IMPORTANT**: The database directory MUST be bind-mounted to be accessible inside the container.

```bash
# Correct: Bind mount the database directory at an IDENTICAL path (host path
# == container path). Databases store internal paths (symlinks, sketch/
# members files) relative to the database dir, so they resolve correctly
# under any bind mapping that's consistent between build and later use --
# an identity bind is the simplest way to guarantee that.
apptainer exec \
    --bind /path/to/databases:/path/to/databases \
    OrthoPhyl.v3.3.sif \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
        --input assemblies.tsv \
        --database-dir /path/to/databases \
        --output-dir results/ \
        --threads 32

# Also wrong, in a subtler way: binding to a DIFFERENT container-side path
# (e.g. --bind /path/to/databases:/databases) works for the run that creates
# it, but any absolute paths baked in at that point won't resolve later from
# the host or from a container invoked with a different bind. Prefer an
# identity bind so host and container paths always match.

# Wrong: Database directory not mounted (will fail!)
apptainer exec OrthoPhyl.v3.3.sif \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
        --database-dir /path/to/databases \  # Not accessible!
        ...
```

### Moving Databases Between Host and Container

As long as every bind used to build and use a database is an **identity bind**
(host path == container path), the database's internal paths are already
relative and portable: `sketch_file`/`members_file`/`metadata_file`,
`source_genome_dir`, `orthophyl_source`, and the `orthophyl_run` symlink are
all stored relative to the database's own directory, and `orthophyl_runs.tsv`
entries are relative to `--database-dir`. A database built inside the
container can be used directly from the host, and vice versa, with no
conversion step.

**Databases built before this (or built with a mismatched, non-identity
bind) may still have absolute, environment-specific paths baked in.** You'll
know this is the problem if the router loads 0 databases from a directory
that clearly has `*_db` subdirectories, or an `orthophyl_run` symlink points
at a path that doesn't exist. Repair them with the migration tool:

```bash
# Preview every move and path rewrite without changing anything:
python assembly_router/OP_database_tool.py --migrate --output-dir /path/to/databases --dry-run

# If the preview looks right, run it for real:
python assembly_router/OP_database_tool.py --migrate --output-dir /path/to/databases
```

`--migrate` infers the old mount prefix automatically from
`database_index.json` (override with `--strip-prefix` if needed), nests
megatree backbone/subclade DBs under their parent-taxon directory
(`<Taxon>/<Taxon>_db`, `<Taxon>/<Taxon>_N_db`), rewrites every path field
relative to the database's new location, relativizes `orthophyl_runs.tsv`,
and regenerates the index. It is idempotent -- running it again is a no-op --
and it **reports** dangling symlink targets rather than deleting anything; a
DB whose `orthophyl_run` points at a deleted `.pipeline_runs/<run>/` dir needs
that run's output restored, or should be removed with
`OP_database_tool.py --remove`.

### Basic Usage Examples

#### Example 1: Batch Mode with Existing Databases

```bash
# Prepare input file (outside container)
cat > assemblies.tsv << EOF
/data/genomes/genome1.fna	d__Bacteria;p__Pseudomonadota;...	genome1
/data/genomes/genome2.fna	d__Bacteria;p__Actinomycetota;...	genome2
EOF

# Run wrapper in container
apptainer exec \
    --bind /data:/data \
    --bind /scratch:/scratch \
    OrthoPhyl.v3.3.sif \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
        --input /data/assemblies.tsv \
        --database-dir /data/databases \
        --output-dir /scratch/results \
        --threads 32
```

#### Example 2: Taxon Mode (Create New Database)

```bash
# Create database from taxon name
apptainer exec \
    --bind /data:/data \
    --bind /scratch:/scratch \
    OrthoPhyl.v3.3.sif \
    bash -c "
        source /opt/conda/etc/profile.d/conda.sh && \
        conda activate gather_genomes && \
        python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
            --taxon 'Methylorubrum' \
            --taxon-rank genus \
            --database-dir /data/databases \
            --output-dir /scratch/methylorubrum_run \
            --gather-script /opt/gits/OrthoPhyl/utils/gather_filter_asms.sh \
            --threads 32
    "
```

**Note**: Taxon mode requires activating the `gather_genomes` conda environment for NCBI tools.

#### Example 3: Interactive Shell

```bash
# Enter container interactively
apptainer shell \
    --bind /data:/data \
    --bind /scratch:/scratch \
    OrthoPhyl.v3.3.sif

# Inside container:
Apptainer> cd /data/my_project
Apptainer> python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py --help
Apptainer> # Run your analysis...
```

#### Example 4: Direct OrthoPhyl.sh Execution

```bash
# Run OrthoPhyl directly (not via wrapper)
apptainer exec \
    --bind /data:/data \
    OrthoPhyl.v3.3.sif \
    bash /opt/gits/OrthoPhyl/OrthoPhyl.sh \
        -g /data/genomes/ \
        -s /data/output/ \
        -t 32 \
        -p iqtree \
        -o CDS
```

### HPC SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=orthophyl_container
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=orthophyl_%j.log

# Set up paths
CONTAINER=/shared/containers/OrthoPhyl.v3.3.sif
DATABASES=/shared/databases
SCRATCH=/scratch/$SLURM_JOB_ID
INPUT=/home/user/assemblies.tsv

# Create scratch directory
mkdir -p $SCRATCH

# Run OrthoPhyl wrapper in container
apptainer exec \
    --bind $DATABASES:$DATABASES \
    --bind $SCRATCH:$SCRATCH \
    --bind $(dirname $INPUT):$(dirname $INPUT) \
    $CONTAINER \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
        --input $INPUT \
        --database-dir $DATABASES \
        --output-dir $SCRATCH/results \
        --threads $SLURM_CPUS_PER_TASK \
        -v

# Copy results to permanent storage
cp -r $SCRATCH/results/03_results /home/user/final_results/

# Cleanup
rm -rf $SCRATCH
```

### Environment Variables in Container

The container sets these environment variables:

```bash
PATH=/opt/conda/bin:/opt/gits/OrthoPhyl:/opt/gits/OrthoPhyl/utils:/opt/gits/OrthoPhyl/assembly_router:$PATH
Path_to_gits=/opt/gits
OrthoPhyl=bash /opt/gits/OrthoPhyl/OrthoPhyl.sh
orthophyl_pipeline_wrapper=python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py
assembly_router=python /opt/gits/OrthoPhyl/assembly_router/assembly_router.py
CHECKM2DB=/opt/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd
ASTRAL_cmd=/opt/gits/ASTRAL/Astral/astral.5.7.8.jar
catfasta2phyml_cmd=/opt/gits/catfasta2phyml/catfasta2phyml.pl
Alignment_Assessment=/opt/gits/Alignment_Assessment/Alignment_Assessment_v2.py
```

You can use these shortcuts (they're `export`ed by `%environment`, so they're
only usable from inside a shell, e.g. via `apptainer shell` or
`apptainer exec ... bash -c '$OrthoPhyl -h'`):

```bash
apptainer shell OrthoPhyl.v3.3.sif
Apptainer> $OrthoPhyl -h
Apptainer> $assembly_router --help
```

### Bind Mount Best Practices

1. **Always bind mount**:
   - Database directory
   - Input genome directory
   - Output directory
   - Any reference data

2. **Use absolute paths** in bind mounts:
   ```bash
   --bind /absolute/path:/absolute/path
   ```

3. **Bind parent directories** if paths are complex:
   ```bash
   --bind /data:/data  # Covers /data/genomes, /data/databases, etc.
   ```

4. **Check what's automatically bound**:
   - `$HOME` (usually)
   - `/tmp` (usually)
   - Current working directory (usually)

5. **Verify mounts** before long runs:
   ```bash
   apptainer exec --bind /data:/data OrthoPhyl.v3.3.sif ls /data
   ```

6. **Point `APPTAINER_TMPDIR`/`SINGULARITY_TMPDIR` at local disk, not NFS**, on
   hosts that lack unprivileged user namespaces/FUSE (every `exec`/`run` call
   prints "Converting SIF file to temporary sandbox..." on these hosts,
   re-extracting the whole multi-GB image each time). NFS-backed tmpdirs can
   cause slow extractions and isolated subprocess failures mid-run -- see
   Troubleshooting below. If local disk isn't available either, `/dev/shm`
   (tmpfs, RAM-backed) is a usable fallback on nodes with enough free memory
   for the image's uncompressed size.

### Running Tests in the Container

OrthoPhyl includes a comprehensive test suite that can be run inside the
container. The test suite is container-aware and handles temporary
directories intelligently.

#### Quick Start

```bash
# Unit tests (fast, ~2 minutes)
apptainer exec \
    --pwd /opt/gits/OrthoPhyl \
    OrthoPhyl.v3.3.sif \
    pytest tests/unit/ -v

# Integration tests (slow, ~30 minutes, requires setup)
apptainer exec \
    --bind /scratch:/scratch \
    --pwd /opt/gits/OrthoPhyl \
    --env TMPDIR=/scratch/pytest_tmp \
    --env ORTHOPHYL_RUN_INTEGRATION=1 \
    OrthoPhyl.v3.3.sif \
    pytest -m integration -v
```

#### Test Categories

**Unit Tests** (`tests/unit/`):
- Fast (< 5 minutes total)
- No external dependencies
- Test individual components
- Safe to run anywhere

**Integration Tests** (`tests/integration/`):
- Slow (30+ minutes)
- Require full OrthoPhyl environment
- Test complete workflows
- Require `ORTHOPHYL_RUN_INTEGRATION=1` to enable (or `-m integration`)

#### Temporary Directory Setup

Containers often have read-only or limited `/tmp`. The test suite handles
this automatically, in priority order:

1. `PYTEST_TMP_DIR` (explicit override)
2. `TMPDIR` / `TMP` / `TEMP` (standard)
3. `.pytest_tmp/` in current directory (fallback)

```bash
# Create writable temp directory
mkdir -p /scratch/pytest_tmp

# Set TMPDIR
apptainer exec \
    --bind /scratch:/scratch \
    --env TMPDIR=/scratch/pytest_tmp \
    OrthoPhyl.v3.3.sif \
    pytest tests/unit/
```

#### Working Directory

Tests expect to run from the OrthoPhyl repository root:

```bash
# Correct: Set working directory
apptainer exec --pwd /opt/gits/OrthoPhyl OrthoPhyl.v3.3.sif pytest

# Wrong: Run from different directory
apptainer exec OrthoPhyl.v3.3.sif pytest  # May fail!
```

#### Test Data Access

Test data is in `TESTER/` directory. Ensure it's accessible:

```bash
# If using external OrthoPhyl repo (not container's copy)
apptainer exec \
    --bind /home/user/OrthoPhyl:/work \
    --pwd /work \
    OrthoPhyl.v3.3.sif \
    pytest tests/unit/
```

#### Container Detection

The test suite automatically detects when running in a container by checking:
- `/.singularity.d/` directory
- `SINGULARITY_CONTAINER` / `APPTAINER_CONTAINER` environment variable
- `/run/.containerenv` file (Podman/Apptainer)
- `/proc/1/cgroup` contents (generic container detection)

When detected, you'll see:
```
================================================================================
CONTAINER DETECTED: Running in container-aware mode
================================================================================
```

#### More Testing Examples

```bash
# Specific test file
apptainer exec --pwd /opt/gits/OrthoPhyl OrthoPhyl.v3.3.sif \
    pytest tests/unit/test_wrapper_batch.py -v

# Parallel execution (requires pytest-xdist in the container)
apptainer exec --pwd /opt/gits/OrthoPhyl OrthoPhyl.v3.3.sif \
    pytest tests/unit/ -n 4 -v

# Exclude slow tests
apptainer exec --pwd /opt/gits/OrthoPhyl OrthoPhyl.v3.3.sif \
    pytest -m "not slow"
```

#### HPC SLURM Test Job

```bash
#!/bin/bash
#SBATCH --job-name=orthophyl_tests
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=test_%j.log

# Setup
CONTAINER=/shared/containers/OrthoPhyl.v3.3.sif
TMPDIR=/scratch/$SLURM_JOB_ID/pytest_tmp
mkdir -p $TMPDIR

# Run tests
apptainer exec \
    --bind /scratch:/scratch \
    --pwd /opt/gits/OrthoPhyl \
    --env TMPDIR=$TMPDIR \
    --env ORTHOPHYL_RUN_INTEGRATION=1 \
    $CONTAINER \
    pytest -m integration -v --tb=short

# Cleanup
rm -rf /scratch/$SLURM_JOB_ID
```

#### Test Selection Reference

```bash
pytest tests/unit/                                            # unit only
pytest -m integration                                          # integration only
pytest tests/unit/test_router.py::TestMultiDatabaseRouter       # a class
pytest tests/unit/test_router.py::TestMultiDatabaseRouter::test_load_databases  # one test
pytest -m "not slow"                                             # skip slow tests
pytest -v                                                        # verbose
pytest -s                                                        # show prints
pytest -x                                                        # stop on first failure
```

#### Testing Best Practices

1. **Always set TMPDIR** explicitly in containers
2. **Use --pwd** to set working directory
3. **Bind mount scratch space** for temp files
4. **Run unit tests first** to verify setup
5. **Use -v flag** for better error messages
6. **Clean up temp directories** after tests
7. **Check available disk space** before integration tests

See also `tests/README.md` for the general (non-container) testing guide and
`pytest.ini` for configuration.

### Troubleshooting

#### Issue: "Database not found"

**Cause**: Database directory not bind-mounted

**Solution**:
```bash
# Add bind mount (identity bind -- see "Critical: Database Directory
# Mounting" above for why host and container paths should match)
apptainer exec --bind /path/to/databases:/path/to/databases ...
```

#### Issue: router loads 0 databases from a directory that clearly has `*_db` dirs in it

**Cause**: the databases were built with absolute, environment-specific paths
baked in -- usually from a pre-portability-fix version of OrthoPhyl, or from
a non-identity bind mount (`--bind /host/dbs:/databases`).

**Solution**: run `OP_database_tool.py --migrate` (see "Moving Databases
Between Host and Container" above).

#### Issue: a database's `orthophyl_run` symlink (or a subclade's sketch/
members file) is broken -- points at a path that doesn't exist

**Cause**: same root cause as above (environment-mismatched absolute path),
or the underlying `.pipeline_runs/<run>/` output was deleted -- databases
still depend on that run's HMMs/alignments; they aren't copied into the DB
itself.

**Solution**: try `--migrate` first. If the link is dangling because the run
output is genuinely gone, either restore it or remove the now-broken database
with `OP_database_tool.py --remove`. To repair a single link manually:
```bash
# From inside the environment where the real orthophyl output is visible:
cd /path/to/databases/<Taxon>_db
ln -sfn "$(realpath --relative-to=. /real/path/to/orthophyl_run)" orthophyl_run
```
Going forward, always bind-mount the database directory at an identical
host/container path to avoid this.

#### Issue: "Permission denied" writing output

**Cause**: Output directory not writable or not bind-mounted

**Solution**:
```bash
# Ensure output directory exists and is bind-mounted
mkdir -p /scratch/output
apptainer exec --bind /scratch:/scratch ...
```

#### Issue: "Command not found"

**Cause**: Using relative paths or wrong environment

**Solution**:
```bash
# Use absolute paths to scripts
python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py

# Or for taxon mode, activate gather_genomes environment
bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate gather_genomes && ..."
```

#### Issue: "NCBI datasets not found" (taxon mode)

**Cause**: Not using gather_genomes conda environment

**Solution**:
```bash
apptainer exec OrthoPhyl.v3.3.sif \
    bash -c "
        source /opt/conda/etc/profile.d/conda.sh && \
        conda activate gather_genomes && \
        python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py --taxon ...
    "
```

#### Issue: isolated subprocess failures mid-run (e.g. "mash: error while loading shared libraries: libgsl.so.25: cannot open shared object file") on a host with no missing library, no OOM, and no resource contention

**Cause**: `APPTAINER_TMPDIR`/`SINGULARITY_TMPDIR` (not the pipeline's own
`TMPDIR` -- see above) points at NFS or other shared/networked storage. On
hosts where Apptainer can't do its normal fast squashfs mount (no
unprivileged user namespaces / FUSE), *every* `exec`/`run` call re-extracts
the whole SIF into a fresh temporary sandbox first -- you'll see `INFO:
Converting SIF file to temporary sandbox...` on each invocation. This image's
sandbox is large (conda envs + the baked-in CheckM2 database, likely 15-20 GB
uncompressed). Extracting that onto NFS is slow, and NFS client-side
attribute/data caching can make an individual file (e.g. a `.so` a subprocess
`dlopen()`s moments after "successful" extraction) appear briefly missing or
stale -- producing an isolated single-call failure rather than every call
failing identically. The same misconfiguration can also starve the node of
memory during the (slow, NFS-bound) extraction step itself, surfacing as an
unrelated-looking OOM-killed `unsquashfs` on a different run.

**Diagnosis**: confirm every `exec`/`run` call prints "Converting SIF file to
temporary sandbox..." (true fast-mount hosts only do this once, or never),
then check where extraction lands:
```bash
echo "APPTAINER_TMPDIR=$APPTAINER_TMPDIR SINGULARITY_TMPDIR=$SINGULARITY_TMPDIR TMPDIR=$TMPDIR"
mount | grep -E "$(dirname "${APPTAINER_TMPDIR:-${SINGULARITY_TMPDIR:-${TMPDIR:-/tmp}}}")"
```
If that mount is `nfs`/`nfs4`, or `df -T` on the directory shows a network
filesystem type, this is almost certainly the cause.

**Solution**: point `APPTAINER_TMPDIR` (and/or `SINGULARITY_TMPDIR`, for the
`singularity` binary) at **local node disk**, not NFS or shared scratch:
```bash
export APPTAINER_TMPDIR=/local/scratch/$USER   # whatever local disk exists on the node
mkdir -p "$APPTAINER_TMPDIR"
apptainer exec --bind /path/to/databases:/databases OrthoPhyl.v3.3.sif \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py ...
```
Note this is a *separate* setting from the pipeline/pytest `TMPDIR` covered
earlier in this doc -- `APPTAINER_TMPDIR`/`SINGULARITY_TMPDIR` governs where
Apptainer/Singularity itself unpacks the container image, before your command
inside it ever runs.

#### "No space left on device" in `/tmp` during tests

**Solution**: Use a location with more space:
```bash
export TMPDIR=/scratch/large_tmp
mkdir -p $TMPDIR
```

---

[← Back to Main README](../README.md)
