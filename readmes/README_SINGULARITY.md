# Singularity Containers - OrthoPhyl Wrapper v2

[← Back to Main README](../README.md)

---


### Container Overview

OrthoPhyl is available as a Singularity container that includes all dependencies pre-installed. This is ideal for:
- HPC environments without conda/mamba
- Reproducible analyses
- Systems where you lack admin privileges
- Ensuring consistent software versions

### Getting the Container

Pull the pre-built container from Sylabs Cloud Library:

```bash
# Create directory for containers (optional but recommended)
mkdir -p ~/singularity_images
cd ~/singularity_images

# Pull the latest OrthoPhyl container
singularity pull library://earlyevol/default/orthophyl

# Rename to version-specific name
mv orthophyl_latest.sif OrthoPhyl.v3.1.0.sif
```

**Note**: The container is ~2.7 GB and includes all dependencies. No building or admin privileges required!

### Container Structure

The container includes:
- **Base environment**: OrthoPhyl, OrthoFinder, IQ-TREE, MAFFT, trimAl, FastTree, HMMER, etc.
- **gather_genomes environment**: CheckM, bbmap, NCBI datasets CLI, entrez-direct
- **Additional tools**: ASTRAL, catfasta2phyml, Alignment Assessment
- **OrthoPhyl code**: Cloned from GitHub at `/opt/gits/OrthoPhyl/`

### Critical: Database Directory Mounting

**IMPORTANT**: The database directory MUST be bind-mounted to be accessible inside the container.

```bash
# Correct: Bind mount the database directory
singularity exec \
    --bind /path/to/databases:/databases \
    OrthoPhyl.v3.1.0.sif \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
        --input assemblies.tsv \
        --database-dir /databases \
        --output-dir results/ \
        --threads 32

# Wrong: Database directory not mounted (will fail!)
singularity exec OrthoPhyl.v3.1.0.sif \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
        --database-dir /path/to/databases \  # Not accessible!
        ...
```

### Basic Usage Examples

#### Example 1: Batch Mode with Existing Databases

```bash
# Prepare input file (outside container)
cat > assemblies.tsv << EOF
/data/genomes/genome1.fna	d__Bacteria;p__Pseudomonadota;...	genome1
/data/genomes/genome2.fna	d__Bacteria;p__Actinomycetota;...	genome2
EOF

# Run wrapper in container
singularity exec \
    --bind /data:/data \
    --bind /scratch:/scratch \
    OrthoPhyl.v3.1.0.sif \
    python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py \
        --input /data/assemblies.tsv \
        --database-dir /data/databases \
        --output-dir /scratch/results \
        --threads 32
```

#### Example 2: Taxon Mode (Create New Database)

```bash
# Create database from taxon name
singularity exec \
    --bind /data:/data \
    --bind /scratch:/scratch \
    OrthoPhyl.v3.1.0.sif \
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
singularity shell \
    --bind /data:/data \
    --bind /scratch:/scratch \
    OrthoPhyl.v3.1.0.sif

# Inside container:
Singularity> cd /data/my_project
Singularity> python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py --help
Singularity> # Run your analysis...
```

#### Example 4: Direct OrthoPhyl.sh Execution

```bash
# Run OrthoPhyl directly (not via wrapper)
singularity exec \
    --bind /data:/data \
    OrthoPhyl.v3.1.0.sif \
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
CONTAINER=/shared/containers/OrthoPhyl.v3.1.0.sif
DATABASES=/shared/databases
SCRATCH=/scratch/$SLURM_JOB_ID
INPUT=/home/user/assemblies.tsv

# Create scratch directory
mkdir -p $SCRATCH

# Run OrthoPhyl wrapper in container
singularity exec \
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
assembly_router=python /opt/gits/OrthoPhyl/assembly_router/assembly_router.py
ASTRAL_cmd=/opt/gits/ASTRAL/Astral/astral.5.7.8.jar
```

You can use these shortcuts:

```bash
# Use environment variable shortcuts
singularity exec OrthoPhyl.v3.1.0.sif $OrthoPhyl -h
singularity exec OrthoPhyl.v3.1.0.sif $assembly_router --help
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
   singularity exec --bind /data:/data OrthoPhyl.v3.1.0.sif ls /data
   ```

### Running Tests in Container

See the dedicated section below for comprehensive testing documentation.

**Quick test**:

```bash
# Unit tests (fast)
singularity exec \
    --pwd /opt/gits/OrthoPhyl \
    OrthoPhyl.v3.1.0.sif \
    pytest tests/unit/ -v

# Integration tests (requires setup)
singularity exec \
    --bind /scratch:/scratch \
    --pwd /opt/gits/OrthoPhyl \
    --env TMPDIR=/scratch/pytest_tmp \
    --env ORTHOPHYL_RUN_INTEGRATION=1 \
    OrthoPhyl.v3.1.0.sif \
    pytest -m integration -v
```

### Troubleshooting Container Issues

#### Issue: "Database not found"

**Cause**: Database directory not bind-mounted

**Solution**:
```bash
# Add bind mount
singularity exec --bind /path/to/databases:/databases ...
```

#### Issue: "Permission denied" writing output

**Cause**: Output directory not writable or not bind-mounted

**Solution**:
```bash
# Ensure output directory exists and is bind-mounted
mkdir -p /scratch/output
singularity exec --bind /scratch:/scratch ...
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
singularity exec OrthoPhyl.v3.1.0.sif \
    bash -c "
        source /opt/conda/etc/profile.d/conda.sh && \
        conda activate gather_genomes && \
        python /opt/gits/OrthoPhyl/orthophyl_pipeline_wrapper.py --taxon ...
    "
```

---

## Running Tests in Containers

### Overview

OrthoPhyl includes a comprehensive test suite that can be run inside Singularity containers. The test suite is container-aware and handles temporary directories intelligently.

### Quick Start

```bash
# Unit tests (fast, ~2 minutes)
singularity exec \
    --pwd /opt/gits/OrthoPhyl \
    OrthoPhyl.v3.1.0.sif \
    pytest tests/unit/ -v

# Integration tests (slow, ~30 minutes, requires setup)
singularity exec \
    --bind /scratch:/scratch \
    --pwd /opt/gits/OrthoPhyl \
    --env TMPDIR=/scratch/pytest_tmp \
    --env ORTHOPHYL_RUN_INTEGRATION=1 \
    OrthoPhyl.v3.1.0.sif \
    pytest -m integration -v
```

### Test Categories

**Unit Tests** (`tests/unit/`):
- Fast (< 5 minutes total)
- No external dependencies
- Test individual components
- Safe to run anywhere

**Integration Tests** (`tests/integration/`):
- Slow (30+ minutes)
- Require full OrthoPhyl environment
- Test complete workflows
- Require `ORTHOPHYL_RUN_INTEGRATION=1` to enable

### Container-Specific Considerations

#### 1. Temporary Directory Setup

Containers often have read-only or limited `/tmp`. The test suite handles this automatically:

**Priority order**:
1. `PYTEST_TMP_DIR` (explicit override)
2. `TMPDIR` / `TMP` / `TEMP` (standard)
3. `.pytest_tmp/` in current directory (fallback)

**Recommended approach**:

```bash
# Create writable temp directory
mkdir -p /scratch/pytest_tmp

# Set TMPDIR
singularity exec \
    --bind /scratch:/scratch \
    --env TMPDIR=/scratch/pytest_tmp \
    OrthoPhyl.v3.1.0.sif \
    pytest tests/unit/
```

#### 2. Working Directory

Tests expect to run from the OrthoPhyl repository root:

```bash
# Correct: Set working directory
singularity exec --pwd /opt/gits/OrthoPhyl OrthoPhyl.v3.1.0.sif pytest

# Wrong: Run from different directory
singularity exec OrthoPhyl.v3.1.0.sif pytest  # May fail!
```

#### 3. Test Data Access

Test data is in `TESTER/` directory. Ensure it's accessible:

```bash
# If using external OrthoPhyl repo (not container's copy)
singularity exec \
    --bind /home/user/OrthoPhyl:/work \
    --pwd /work \
    OrthoPhyl.v3.1.0.sif \
    pytest tests/unit/
```

### Complete Testing Examples

#### Example 1: Basic Unit Tests

```bash
singularity exec \
    --pwd /opt/gits/OrthoPhyl \
    OrthoPhyl.v3.1.0.sif \
    pytest tests/unit/ -v --tb=short
```

#### Example 2: Integration Tests with Scratch Space

```bash
# Create temp directory
mkdir -p /scratch/$USER/pytest_tmp

# Run integration tests
singularity exec \
    --bind /scratch:/scratch \
    --pwd /opt/gits/OrthoPhyl \
    --env TMPDIR=/scratch/$USER/pytest_tmp \
    --env ORTHOPHYL_RUN_INTEGRATION=1 \
    OrthoPhyl.v3.1.0.sif \
    pytest -m integration -v

# Cleanup
rm -rf /scratch/$USER/pytest_tmp
```

#### Example 3: Specific Test File

```bash
singularity exec \
    --pwd /opt/gits/OrthoPhyl \
    OrthoPhyl.v3.1.0.sif \
    pytest tests/unit/test_wrapper_batch.py -v
```

#### Example 4: HPC SLURM Job

```bash
#!/bin/bash
#SBATCH --job-name=orthophyl_tests
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=test_%j.log

# Setup
CONTAINER=/shared/containers/OrthoPhyl.v3.1.0.sif
TMPDIR=/scratch/$SLURM_JOB_ID/pytest_tmp
mkdir -p $TMPDIR

# Run tests
singularity exec \
    --bind /scratch:/scratch \
    --pwd /opt/gits/OrthoPhyl \
    --env TMPDIR=$TMPDIR \
    --env ORTHOPHYL_RUN_INTEGRATION=1 \
    $CONTAINER \
    pytest -m integration -v --tb=short

# Cleanup
rm -rf /scratch/$SLURM_JOB_ID
```

#### Example 5: Parallel Test Execution

```bash
# Install pytest-xdist in container (if not already present)
# Then run tests in parallel

singularity exec \
    --pwd /opt/gits/OrthoPhyl \
    OrthoPhyl.v3.1.0.sif \
    pytest tests/unit/ -n 4 -v  # 4 parallel workers
```

### Test Selection

```bash
# Run only unit tests
pytest tests/unit/

# Run only integration tests
pytest -m integration

# Run specific test class
pytest tests/unit/test_router.py::TestMultiDatabaseRouter

# Run specific test method
pytest tests/unit/test_router.py::TestMultiDatabaseRouter::test_load_databases

# Exclude slow tests
pytest -m "not slow"

# Run with verbose output
pytest -v

# Show print statements
pytest -s

# Stop on first failure
pytest -x
```

### Container Detection

The test suite automatically detects when running in a container:

```
CONTAINER DETECTED: Running in container-aware mode
Temp directory: /scratch/pytest_tmp
```

Detection checks for:
- `/.singularity.d/` directory
- `SINGULARITY_CONTAINER` environment variable
- `/.dockerenv` file
- `/run/.containerenv` file

### Troubleshooting Test Issues

#### "Permission denied" in /tmp

**Solution**:
```bash
export TMPDIR=/scratch/pytest_tmp
mkdir -p $TMPDIR
```

#### "No space left on device"

**Solution**: Use larger scratch space:
```bash
export TMPDIR=/scratch/large_space
```

#### "TESTER/genomes not found"

**Solution**: Ensure correct working directory:
```bash
singularity exec --pwd /opt/gits/OrthoPhyl ...
```

#### Tests hang or timeout

**Solution**: Check if integration tests are running (they're slow):
```bash
# Run only fast unit tests
pytest tests/unit/ -v
```

### Best Practices

1. **Always set TMPDIR** explicitly in containers
2. **Use --pwd** to set working directory
3. **Bind mount scratch space** for temp files
4. **Run unit tests first** to verify setup
5. **Use -v flag** for better error messages
6. **Clean up temp directories** after tests
7. **Check available disk space** before integration tests

### Performance Tips

1. **Use node-local storage** (`/scratch` on HPC)
2. **Run tests in parallel** with `-n` flag
3. **Skip slow tests** during development: `-m "not slow"`
4. **Use session-scoped fixtures** (already implemented)

### Additional Resources

- **Detailed container testing guide**: `tests/CONTAINER_TESTING.md`
- **General testing guide**: `tests/README.md`
- **Test configuration**: `pytest.ini`

---

---

[← Back to Main README](../README.md)
