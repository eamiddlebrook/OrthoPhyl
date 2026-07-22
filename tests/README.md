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
