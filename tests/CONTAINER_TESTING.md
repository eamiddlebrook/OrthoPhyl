# Running Tests in Singularity/Docker Containers

This guide explains how to run OrthoPhyl tests inside Singularity or Docker containers, which have special considerations for temporary directories and file permissions.

## Quick Start

### Singularity Container

```bash
# Basic unit tests (no special setup needed)
singularity exec OrthoPhyl.v3.1.0.sif pytest tests/unit/

# Integration tests with writable temp directory
singularity exec \
  --bind /scratch:/scratch \
  --env TMPDIR=/scratch/pytest_tmp \
  --env ORTHOPHYL_RUN_INTEGRATION=1 \
  OrthoPhyl.v3.1.0.sif \
  pytest -m integration

# Alternative: Use explicit PYTEST_TMP_DIR
singularity exec \
  --bind /scratch:/scratch \
  --env PYTEST_TMP_DIR=/scratch/pytest_tmp \
  --env ORTHOPHYL_RUN_INTEGRATION=1 \
  OrthoPhyl.v3.1.0.sif \
  pytest -m integration
```

### Docker Container

```bash
# Basic unit tests
docker run --rm -v $(pwd):/work -w /work orthophyl:latest pytest tests/unit/

# Integration tests with volume mount
docker run --rm \
  -v $(pwd):/work \
  -v /tmp/pytest_tmp:/tmp/pytest_tmp \
  -e TMPDIR=/tmp/pytest_tmp \
  -e ORTHOPHYL_RUN_INTEGRATION=1 \
  -w /work \
  orthophyl:latest \
  pytest -m integration
```

## Container-Specific Issues

### 1. Temporary Directory Handling

**Problem:** Containers often have read-only `/tmp` or bind-mounted `/tmp` that behaves differently than expected.

**Solution:** The test suite now automatically detects containers and handles temporary directories intelligently:

1. **PYTEST_TMP_DIR** (highest priority): Explicit override
   ```bash
   export PYTEST_TMP_DIR=/scratch/my_tests
   ```

2. **TMPDIR/TMP/TEMP** (standard): Respects standard environment variables
   ```bash
   export TMPDIR=/scratch/tmp
   ```

3. **Fallback**: If `/tmp` is not writable, creates `.pytest_tmp/` in current directory
   ```bash
   # Automatically used if above options fail
   # Creates: /work/.pytest_tmp/
   ```

### 2. Container Detection

The test suite automatically detects when running in a container by checking:
- `/.singularity.d` directory (Singularity)
- `SINGULARITY_CONTAINER` environment variable
- `/.dockerenv` file (Docker)
- `/run/.containerenv` file (Podman)
- `/proc/1/cgroup` contents (generic container detection)

When detected, you'll see:
```
================================================================================
CONTAINER DETECTED: Running in container-aware mode
================================================================================
```

### 3. File Permissions

**Singularity:** By default, runs as your user, so permissions usually work.

**Docker:** May run as root, causing permission issues with output files.

**Solution:** Use `--user` flag:
```bash
docker run --rm --user $(id -u):$(id -g) -v $(pwd):/work -w /work ...
```

### 4. Path Binding/Mounting

**Singularity:** Automatically binds `$HOME`, `/tmp`, and current directory. For other paths:
```bash
singularity exec --bind /scratch:/scratch --bind /data:/data ...
```

**Docker:** Must explicitly mount all needed paths:
```bash
docker run -v $(pwd):/work -v /data:/data -w /work ...
```

## Environment Variables

### Test Control
- `ORTHOPHYL_RUN_INTEGRATION=1`: Enable integration tests (disabled by default)
- `PYTEST_TMP_DIR=/path`: Override temporary directory location
- `TMPDIR=/path`: Standard temp directory (fallback to PYTEST_TMP_DIR)

### Container-Specific
- `SINGULARITY_CONTAINER`: Set by Singularity (read-only)
- `SINGULARITY_BIND`: Comma-separated bind paths

## Common Scenarios

### HPC with Singularity (SLURM)

```bash
#!/bin/bash
#SBATCH --job-name=orthophyl_tests
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Use node-local scratch for temp files
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p $TMPDIR

# Run tests
singularity exec \
  --bind $TMPDIR:$TMPDIR \
  --env TMPDIR=$TMPDIR \
  --env ORTHOPHYL_RUN_INTEGRATION=1 \
  /path/to/OrthoPhyl.sif \
  pytest -m integration -v

# Cleanup
rm -rf $TMPDIR
```

### CI/CD with Docker

```yaml
# .gitlab-ci.yml or similar
test:
  image: orthophyl:latest
  variables:
    PYTEST_TMP_DIR: /tmp/pytest
    ORTHOPHYL_RUN_INTEGRATION: "1"
  before_script:
    - mkdir -p /tmp/pytest
  script:
    - pytest -m integration -v
  artifacts:
    when: always
    paths:
      - .pytest_tmp/  # Fallback location if /tmp fails
```

### Local Development with Singularity

```bash
# Create a persistent temp directory
mkdir -p ~/orthophyl_test_tmp

# Run tests with bind mount
singularity exec \
  --bind ~/orthophyl_test_tmp:/tmp \
  --pwd /home/earlm/gits/OrthoPhyl \
  OrthoPhyl.sif \
  pytest tests/unit/ -v

# Integration tests (longer running)
singularity exec \
  --bind ~/orthophyl_test_tmp:/tmp \
  --pwd /home/earlm/gits/OrthoPhyl \
  --env ORTHOPHYL_RUN_INTEGRATION=1 \
  OrthoPhyl.sif \
  pytest -m integration -v
```

## Troubleshooting

### "Permission denied" writing to /tmp

**Symptom:**
```
PermissionError: [Errno 13] Permission denied: '/tmp/pytest-...'
```

**Solution:**
```bash
# Option 1: Set TMPDIR to writable location
export TMPDIR=/scratch/tmp
mkdir -p $TMPDIR

# Option 2: Use PYTEST_TMP_DIR
export PYTEST_TMP_DIR=/scratch/pytest
mkdir -p $PYTEST_TMP_DIR

# Option 3: Let pytest use CWD fallback (automatic)
# Creates .pytest_tmp/ in current directory
```

### "No space left on device" in /tmp

**Symptom:**
```
OSError: [Errno 28] No space left on device
```

**Solution:** Use a location with more space:
```bash
export TMPDIR=/scratch/large_tmp
mkdir -p $TMPDIR
```

### Tests can't find input files

**Symptom:**
```
FileNotFoundError: TESTER/genomes_fasttest not found
```

**Solution:** Ensure repository is properly mounted:
```bash
# Singularity: Use --pwd to set working directory
singularity exec --pwd /path/to/OrthoPhyl OrthoPhyl.sif pytest

# Docker: Use -w flag
docker run -v $(pwd):/work -w /work orthophyl:latest pytest
```

### Container not detected (should be but isn't)

**Symptom:** Tests fail with container-specific errors but no "CONTAINER DETECTED" message.

**Solution:** Manually set indicator:
```bash
export SINGULARITY_CONTAINER=1  # For Singularity
# or
touch /.dockerenv  # For Docker (requires root)
```

## Test Selection

```bash
# Run only unit tests (fast, no container issues)
pytest tests/unit/

# Run only integration tests (requires full environment)
pytest -m integration

# Run specific test file
pytest tests/unit/test_router.py

# Run specific test
pytest tests/unit/test_router.py::TestRouter::test_basic_routing

# Exclude slow tests
pytest -m "not slow"

# Verbose output
pytest -v

# Show print statements
pytest -s
```

## Performance Tips

1. **Use node-local storage** for temp files (e.g., `/scratch` on HPC)
2. **Bind mount** only necessary paths to reduce overhead
3. **Session-scoped fixtures** run expensive operations once (already implemented)
4. **Parallel execution** with pytest-xdist:
   ```bash
   pytest -n 4  # Run 4 tests in parallel
   ```

## Best Practices

1. **Always set TMPDIR** explicitly in container environments
2. **Use absolute paths** for bind mounts
3. **Check available space** in temp directory before long tests
4. **Clean up** temp directories after tests (especially in CI)
5. **Use --tb=short** for cleaner error messages (already in pytest.ini)

## Getting Help

If tests fail in containers:

1. Check the "CONTAINER DETECTED" message appears
2. Verify temp directory is writable: `touch $TMPDIR/test && rm $TMPDIR/test`
3. Check disk space: `df -h $TMPDIR`
4. Run with verbose output: `pytest -v -s`
5. Check the test logs in `.pytest_tmp/` or `$PYTEST_TMP_DIR/`

For issues specific to OrthoPhyl tests, see `tests/README.md`.
