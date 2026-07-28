# OrthoPhyl wrapper/router test suite

Unit tests for the pipeline-wrapper / assembly-router layer. See
[`../WRAPPER_UNIT_TEST_PLAN.md`](../WRAPPER_UNIT_TEST_PLAN.md) for the full plan and the
list of known bugs these tests target.

## Setup

The test dependencies (`pytest`, `pytest-cov`, `pytest-mock`) ship with the OrthoPhyl
conda environment (`orthophyl_env.2.2.1.yml`), so if you have that activated you already
have everything:

```bash
conda activate OrthoPhyl
pytest tests/
```

If you prefer an isolated environment (or don't use conda), build a venv from the
requirements file instead:

```bash
python3 -m venv .venv-test
.venv-test/bin/pip install -r tests/requirements-test.txt
.venv-test/bin/python -m pytest tests/
```

## Running

```bash
# whole suite
pytest tests/

# one module
pytest tests/unit/test_router.py -v

# skip slow / network-dependent tests (none yet, but the marker is reserved)
pytest tests/ -m "not slow"

# coverage
pytest tests/ --cov=. --cov-report=term-missing
```

(Prefix with `.venv-test/bin/python -m ` if you used the venv route above.)

## Layout

```
tests/
├── conftest.py                 # module loader (dotted filenames) + shared fixtures
├── requirements-test.txt
├── unit/
│   ├── test_gtdb_taxonomy.py   # §1  GTDBTaxonomy (both copies, parametrized)
│   ├── test_router.py          # §2  MultiDatabaseRouter routing decisions
│   ├── test_db_creator.py      # §3  create_hierarchical_database_v2
│   └── test_wrapper_batch.py   # §4  wrapper orchestration (subprocess mocked)
└── fixtures/                   # static fixture data (added as suites grow)
```

## How dotted module names are imported

`assembly_router_multi.cmd_out3.py` and `orthophyl_pipeline_wrapper.v2.py` are not legal
Python module identifiers, so they cannot be `import`ed normally. `conftest.py` provides
a `load_module(path)` helper (via `importlib`) and exposes each target as a
session-scoped fixture (`router_module`, `wrapper_module`, `db_creator_module`, …).

## Regression tests for fixed bugs (B1–B4)

Bugs B1–B4 from [`../WRAPPER_UNIT_TEST_PLAN.md`](../WRAPPER_UNIT_TEST_PLAN.md) are fixed;
each has a guarding regression test:

| Test | Bug (fixed) |
|------|-------------|
| `test_wrapper_batch.py::TestReleafPhase::test_releaf_actually_invoked_when_not_dry_run` | **B1** — misindented `return` in `_run_releaf` meant `ReLeaf.sh` was never executed |
| `test_wrapper_batch.py::TestReleafPhase::test_version_creation_runs_when_db_found` | **B2/B3** — unconditional `return`s made `_create_releaf_version` dead code |
| `test_taxon_assembly_gatherer.py::TestWrapperContract` | **B4** — wrapper↔gatherer API mismatch (ctor kwarg, missing methods, wrong attr) |

The `xfail(strict=True)` pattern is still the recommended approach for *future* bugs
found before their fix: mark the intended-behavior test xfail so the suite stays green
today, then it XPASSes (a hard failure under `strict=True`) once fixed, prompting removal
of the marker.

Still open (documented in the plan, not yet fixed): the `.fna`-only genome count in
`_run_orthophyl`, and the hardcoded `added_genomes = 0` in
`add_releaf_version.count_genomes_in_releaf`.

---

## Integration tests

**Location:** `tests/integration/`

Integration tests run the **real** `OrthoPhyl.sh` and `ReLeaf.sh` scripts on small test
datasets, asserting correct outputs and topology. See
[`../INTEGRATION_TEST_PLAN.md`](../INTEGRATION_TEST_PLAN.md) for the full plan.

### Why integration tests?

Unit tests (above) mock all subprocess calls, so they catch logic bugs but **cannot detect**:
- Wrong CLI flags passed to real tools (the mock accepts any argv).
- False-success scenarios where a tool fails but the wrapper doesn't detect it.
- Actual tool integration breakage (e.g., iqtree output format changes).

Integration tests fill that gap by running the real bioinformatics pipeline and asserting
on the actual scientific outputs (trees, alignments, taxon presence, **topology
correctness via Robinson-Foulds distance**).

### Running integration tests

Integration tests are **opt-in** (too slow for rapid iteration):

```bash
# Requires the full OrthoPhyl conda environment
conda activate OrthoPhyl

# Run integration tests explicitly
pytest -m integration -v

# Or via environment variable (for CI)
ORTHOPHYL_RUN_INTEGRATION=1 pytest tests/integration/
```

**Why opt-in?** OrthoPhyl takes ~5–15 minutes on the fasttest dataset. Unit tests run in
seconds; integration tests are for nightly CI or pre-release validation.

### Test data

- **Base dataset:** `TESTER/genomes_fasttest` (6 orchid chloroplast genomes, ~150 kb each)
- **Add-assembly dataset:** `TESTER/genomes_fasttest_addasm` (2 additional genomes for ReLeaf)
- **Reference trees:** `TESTER/REFERENCE_TESTER_TREES/` (gold-standard topologies for comparison)

### What's tested

| Suite | Tests |
|-------|-------|
| `test_orthophyl_fasttest.py` | OrthoPhyl.sh end-to-end: exit 0, trees generated, HMMs/alignments present, all taxa in tree, **topology matches reference (RF=0)** |
| `test_releaf_fasttest.py` | ReLeaf.sh add-assembly: exit 0, new alignments/trees, added taxa present, original taxa retained, 8 total taxa |

The **topology comparison** (Robinson-Foulds distance = 0) is the strongest assertion — it
catches not just "a tree was produced" but "the *correct* tree was produced."

### Tool requirements

Integration tests are automatically skipped if required tools are missing:
- `iqtree`, `orthofinder`, `mafft`, `hmmbuild`, `hmmsearch`, `trimal`

Activate the OrthoPhyl conda environment to get all tools:
```bash
conda activate OrthoPhyl
```
