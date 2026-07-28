# OrthoPhyl & ReLeaf — Integration Test Plan

## Purpose

This plan covers **end-to-end integration tests** that invoke the real `OrthoPhyl.sh` and
`ReLeaf.sh` scripts on small test datasets, asserting correct outputs and topology. It is
**complementary to** [`WRAPPER_UNIT_TEST_PLAN.md`](WRAPPER_UNIT_TEST_PLAN.md), which tests
the orchestration/routing layer with all subprocess calls mocked.

**Why both?** Unit tests with mocked `subprocess.run` catch logic bugs (routing decisions,
path wiring, checkpoint handling) but cannot detect:
- Wrong CLI flags passed to the real tools (the mock accepts any argv).
- False-success scenarios where a tool fails but the wrapper doesn't detect it.
- Actual tool integration breakage (e.g., iqtree output format changes).

Integration tests fill that gap by running the real bioinformatics pipeline and asserting
on the actual scientific outputs (trees, alignments, taxon presence, topology correctness).

---

## Scope

| Component | What's tested | Data |
|-----------|---------------|------|
| `OrthoPhyl.sh` | Full pipeline: annotation → OrthoFinder → alignment → tree inference | `TESTER/genomes_fasttest` (6 orchid chloroplasts) |
| `ReLeaf.sh` | Add-assembly workflow: HMM search → alignment → tree update | `TESTER/genomes_fasttest_addasm` (2 additional genomes) |
| Topology correctness | Robinson–Foulds distance vs reference trees | `TESTER/REFERENCE_TESTER_TREES/` |

**Out of scope** (covered by unit tests):
- Wrapper orchestration logic (`orthophyl_pipeline_wrapper.v2.py`).
- Router decision-making (`assembly_router_multi.cmd_out3.py`).
- Database creation (`create_hierarchical_database_v2.py`).

---

## Test data

### Base dataset: `TESTER/genomes_fasttest`
6 orchid chloroplast genomes (~150 kb each):
- `AB893950.1_Dendrobium_moniliforme_chloroplast.fasta`
- `KU551270.1_Neottia_fugongensis_voucher_Jin_X.H._11375_(PE).fasta`
- `NC_026568.1_Cattleya_crispata.fasta`
- `NC_035154.1_Dendrobium_moniliforme_chloroplast.fasta`
- `NC_042686.1_Gastrochilus_calceolaris.fasta`
- `NC_050919.1_Apostasia_ramifera.fasta`

Pre-annotated CDS/protein files in `TESTER/annots_{nucls,prots}_fasttest/`.

### Add-assembly dataset: `TESTER/genomes_fasttest_addasm`
2 additional genomes for ReLeaf:
- `MK836106.1_Holcoglossum_tsii_isolate_LDKAe67.fasta`
- `MK936427.1_Spiranthes_sinensis.fasta`

Pre-annotated in `TESTER/annots_{nucls,prots}_fasttest_addasm/`.

### Reference trees: `TESTER/REFERENCE_TESTER_TREES/`
Gold-standard topologies from a known-good run:
- `fastTree.SCO_strict.codon_aln.tree` — 6-taxon strict SCO tree (FastTree, CDS)
- `fastTree.SCO_6.codon_aln.tree` — 6-taxon relaxed (≥6 taxa per OG)
- `OG_SCO_strict.astral.tree` — ASTRAL coalescent tree (strict SCO)
- `OG_SCO_6.astral.tree` — ASTRAL (relaxed)

**Note:** Reference trees are named for the tool that produced them (`fastTree`, `astral`),
but topology comparison is tool-agnostic (unrooted RF distance). The current pipeline uses
`iqtree` by default; we compare topologies, not branch lengths or tree methods.

---

## Test framework

**Tool:** `pytest` with `ete3` for tree comparison.

**Markers:**
- `@pytest.mark.integration` — all tests in this suite.
- `@pytest.mark.slow` — long-running (OrthoPhyl takes ~5–15 min on fasttest data).

**Opt-in execution:** Integration tests are **not** run by default (too slow for rapid
iteration). Invoke explicitly:
```bash
pytest -m integration
```

Or via environment variable (for CI):
```bash
ORTHOPHYL_RUN_INTEGRATION=1 pytest tests/integration/
```

**Tool availability guards:** Tests are skipped if required tools (`iqtree`, `orthofinder`,
`mafft`, `hmmer`, etc.) are not on `PATH`. This allows the unit test suite to run in
minimal environments while integration tests require the full OrthoPhyl conda environment.

---

## Test suite structure

```
tests/
├── unit/                                   # existing, unchanged
└── integration/
    ├── conftest.py                         # session fixtures + tool guards
    ├── test_orthophyl_fasttest.py          # OrthoPhyl end-to-end
    ├── test_releaf_fasttest.py             # ReLeaf add-assembly
    └── helpers/
        └── tree_compare.py                 # ete3 wrappers: RF distance, topology checks
```

---

## 1. Session fixtures (`tests/integration/conftest.py`)

### `orthophyl_run` (session-scoped)
Runs `OrthoPhyl.sh` **once** at the start of the test session into a shared `tmp_path`,
reused by all tests. Mirrors the invocation from `test.sh`:

```bash
./OrthoPhyl.sh \
    -s <tmp_output_dir> \
    -a ./TESTER/annots_nucls_fasttest,./TESTER/annots_prots_fasttest \
    -o BOTH \
    -R full \
    -g ./TESTER/genomes_fasttest \
    -t 4 \
    -c control_file.user \
    -n 100
```

Returns the output directory path. Fails the entire session if OrthoPhyl exits non-zero.

### Tool availability guards
```python
@pytest.fixture(scope="session", autouse=True)
def require_tools():
    """Skip all integration tests if required tools are missing."""
    required = ["iqtree", "orthofinder", "mafft", "hmmbuild", "hmmsearch"]
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"Integration tests require: {', '.join(missing)}")
```

### Opt-in guard
```python
def pytest_configure(config):
    if not (config.getoption("-m") == "integration" or 
            os.getenv("ORTHOPHYL_RUN_INTEGRATION")):
        config.option.markexpr = "not integration"
```

---

## 2. Tree comparison helpers (`tests/integration/helpers/tree_compare.py`)

Wraps `ete3` for topology assertions:

```python
from ete3 import Tree

def load_tree(path: str) -> Tree:
    """Load a Newick tree, handling common format variations."""
    return Tree(path, format=1)  # flexible format

def normalize_labels(tree: Tree) -> Tree:
    """Strip whitespace/quotes from leaf names for robust comparison."""
    for leaf in tree.iter_leaves():
        leaf.name = leaf.name.strip().strip('"').strip("'")
    return tree

def rf_distance(tree1: Tree, tree2: Tree, unrooted: bool = True) -> tuple[int, int, int, int, int]:
    """
    Compute Robinson-Foulds distance between two trees.
    
    Returns: (rf, max_rf, common_leaves, tree1_leaves, tree2_leaves)
    """
    return tree1.robinson_foulds(tree2, unrooted_trees=unrooted)

def same_topology(tree1: Tree, tree2: Tree, tolerance: int = 0) -> bool:
    """
    Assert trees have identical topology (RF distance ≤ tolerance).
    
    Args:
        tolerance: Allow small RF differences (default 0 = exact match).
    """
    rf, max_rf, common, _, _ = rf_distance(tree1, tree2)
    return rf <= tolerance

def assert_taxa_present(tree: Tree, expected_taxa: list[str]):
    """Assert all expected taxa appear in the tree."""
    leaf_names = {leaf.name for leaf in tree.iter_leaves()}
    missing = set(expected_taxa) - leaf_names
    assert not missing, f"Missing taxa in tree: {missing}"
```

---

## 3. OrthoPhyl integration tests (`test_orthophyl_fasttest.py`)

### 3.1 `test_orthophyl_completes_successfully`
```python
def test_orthophyl_completes_successfully(orthophyl_run):
    """OrthoPhyl.sh exits 0 and produces expected directory structure."""
    assert orthophyl_run.exists()
    assert (orthophyl_run / "genome_list").exists()
    assert (orthophyl_run / "phylo_current").exists()
```

### 3.2 `test_species_trees_generated`
```python
def test_species_trees_generated(orthophyl_run):
    """Final species trees exist for both SCO_strict and SCO_3."""
    tree_dir = orthophyl_run / "phylo_current" / "SpeciesTree"
    
    # Top-level tree copies
    assert (tree_dir / "iqtree.SCO_strict.CDS.tree").exists()
    assert (tree_dir / "iqtree.SCO_3.CDS.tree").exists()
    
    # Detailed iqtree outputs
    assert (tree_dir / "iqtree" / "iqtree.SCO_strict.CDS.treefile").exists()
    assert (tree_dir / "iqtree" / "iqtree.SCO_strict.CDS.iqtree").exists()
```

### 3.3 `test_hmms_and_alignments_present`
```python
def test_hmms_and_alignments_present(orthophyl_run):
    """HMMs and trimmed alignments are generated."""
    # HMMs (multiple possible locations, check both)
    hmm_candidates = [
        orthophyl_run / "OG_alignmentsToHMM" / "hmms_final",
        orthophyl_run / "phylo_current" / "OG_alignmentsToHMM" / "hmms_final",
    ]
    assert any(d.exists() and list(d.glob("*.hmm")) for d in hmm_candidates)
    
    # Trimmed alignments
    assert (orthophyl_run / "phylo_current" / "AlignmentsProts.trm.nm").exists()
    assert (orthophyl_run / "phylo_current" / "AlignmentsCDS.trm.nm").exists()
```

### 3.4 `test_all_input_taxa_in_tree`
```python
def test_all_input_taxa_in_tree(orthophyl_run):
    """All 6 input genomes appear in the final species tree."""
    tree_path = orthophyl_run / "phylo_current" / "SpeciesTree" / "iqtree.SCO_strict.CDS.tree"
    tree = load_tree(str(tree_path))
    
    expected_taxa = [
        "AB893950.1_Dendrobium_moniliforme_chloroplast",
        "KU551270.1_Neottia_fugongensis_voucher_Jin_X.H._11375_PE",
        "NC_026568.1_Cattleya_crispata",
        "NC_035154.1_Dendrobium_moniliforme_chloroplast",
        "NC_042686.1_Gastrochilus_calceolaris",
        "NC_050919.1_Apostasia_ramifera",
    ]
    assert_taxa_present(tree, expected_taxa)
```

### 3.5 `test_topology_matches_reference` ⭐ **Key test**
```python
def test_topology_matches_reference(orthophyl_run):
    """Produced tree has exact topology (RF=0) vs reference tree."""
    produced = load_tree(
        str(orthophyl_run / "phylo_current" / "SpeciesTree" / "iqtree.SCO_strict.CDS.tree")
    )
    reference = load_tree("TESTER/REFERENCE_TESTER_TREES/fastTree.SCO_strict.codon_aln.tree")
    
    produced = normalize_labels(produced)
    reference = normalize_labels(reference)
    
    assert same_topology(produced, reference, tolerance=0), \
        "Tree topology differs from reference (RF > 0)"
```

**Rationale:** This is the strongest assertion — it catches not just "a tree was produced"
but "the *correct* tree was produced." If the pipeline silently uses wrong parameters,
drops taxa, or has a tool integration bug, the topology will diverge.

---

## 4. ReLeaf integration tests (`test_releaf_fasttest.py`)

### 4.1 `test_releaf_completes_successfully`
```python
def test_releaf_completes_successfully(orthophyl_run, tmp_path):
    """ReLeaf.sh exits 0 when adding assemblies to an existing OrthoPhyl run."""
    cmd = [
        "./ReLeaf.sh",
        "-g", str(Path("TESTER/genomes_fasttest_addasm").resolve()),
        "-a", "TESTER/annots_nucls_fasttest_addasm,TESTER/annots_prots_fasttest_addasm",
        "-s", str(orthophyl_run),
        "-t", "4",
        "-p", "iqtree",
        "-o", "BOTH",
    ]
    result = subprocess.run(cmd, cwd=Path.cwd(), capture_output=True, text=True)
    
    assert result.returncode == 0, \
        f"ReLeaf.sh failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
```

### 4.2 `test_releaf_outputs_present`
```python
def test_releaf_outputs_present(orthophyl_run):
    """ReLeaf produces new alignments and trees."""
    releaf_dir = orthophyl_run / "ReLeaf_dir"
    
    # New alignments
    assert (releaf_dir / "new_prot_alignments.trm.nm").exists()
    assert (releaf_dir / "new_CDS_alignments.trm.nm").exists()
    
    # New trees
    tree_dir = releaf_dir / "new_trees"
    assert tree_dir.exists()
    assert list(tree_dir.glob("*.tree*"))  # at least one tree file
```

### 4.3 `test_added_taxa_in_releaf_tree`
```python
def test_added_taxa_in_releaf_tree(orthophyl_run):
    """Both newly-added genomes appear in the ReLeaf output tree."""
    # ReLeaf typically produces phylogeny_with_new_genomes.nwk or similar
    tree_candidates = [
        orthophyl_run / "ReLeaf_dir" / "phylogeny_with_new_genomes.nwk",
        orthophyl_run / "ReLeaf_dir" / "new_trees" / "phylogeny.nwk",
    ]
    tree_path = next((p for p in tree_candidates if p.exists()), None)
    assert tree_path, f"No ReLeaf output tree found in {tree_candidates}"
    
    tree = load_tree(str(tree_path))
    
    added_taxa = [
        "MK836106.1_Holcoglossum_tsii_isolate_LDKAe67",
        "MK936427.1_Spiranthes_sinensis",
    ]
    assert_taxa_present(tree, added_taxa)
```

### 4.4 `test_original_taxa_retained`
```python
def test_original_taxa_retained(orthophyl_run):
    """Original 6 taxa are still present after ReLeaf (no dropouts)."""
    tree_path = orthophyl_run / "ReLeaf_dir" / "phylogeny_with_new_genomes.nwk"
    tree = load_tree(str(tree_path))
    
    original_taxa = [
        "AB893950.1_Dendrobium_moniliforme_chloroplast",
        "KU551270.1_Neottia_fugongensis_voucher_Jin_X.H._11375_PE",
        "NC_026568.1_Cattleya_crispata",
        "NC_035154.1_Dendrobium_moniliforme_chloroplast",
        "NC_042686.1_Gastrochilus_calceolaris",
        "NC_050919.1_Apostasia_ramifera",
    ]
    assert_taxa_present(tree, original_taxa)
```

### 4.5 `test_releaf_tree_has_8_taxa`
```python
def test_releaf_tree_has_8_taxa(orthophyl_run):
    """Final ReLeaf tree contains exactly 8 taxa (6 original + 2 added)."""
    tree_path = orthophyl_run / "ReLeaf_dir" / "phylogeny_with_new_genomes.nwk"
    tree = load_tree(str(tree_path))
    
    leaf_count = len(list(tree.iter_leaves()))
    assert leaf_count == 8, f"Expected 8 taxa, found {leaf_count}"
```

---

## 5. Running the tests

### Local development (opt-in)
```bash
# Activate OrthoPhyl conda environment (provides all tools)
conda activate OrthoPhyl

# Run integration tests only
pytest -m integration -v

# Run with coverage (integration tests don't contribute much to code coverage,
# but can verify the scripts execute without crashing)
pytest -m integration --cov=. --cov-report=term-missing
```

### CI (GitHub Actions / GitLab CI)
Add a separate stage to `.github/workflows/tests.yml`:

```yaml
integration-tests:
  runs-on: ubuntu-latest
  needs: unit-tests  # run after fast unit tests pass
  if: github.event_name == 'schedule' || contains(github.event.head_commit.message, '[integration]')
  steps:
    - uses: actions/checkout@v3
    - name: Set up Conda
      uses: conda-incubator/setup-miniconda@v2
      with:
        environment-file: orthophyl_env.2.2.1.yml
        activate-environment: OrthoPhyl
    - name: Run integration tests
      run: |
        conda activate OrthoPhyl
        pytest -m integration -v --tb=short
      env:
        ORTHOPHYL_RUN_INTEGRATION: "1"
```

**Trigger strategy:**
- **Nightly cron** — catch tool version drift / environment changes.
- **Manual opt-in** — commit message contains `[integration]`.
- **Not on every PR** — too slow for rapid feedback.

---

## 6. Maintenance & evolution

### When to update reference trees
If the pipeline changes in a way that intentionally alters topology (e.g., switching from
FastTree to IQ-TREE with a different model), regenerate reference trees:

```bash
# Run the pipeline on fasttest data
./test.sh 4 100

# Copy the output trees to the reference dir
cp TESTER/FULLTEST_OUT.*/phylo_current/SpeciesTree/iqtree.SCO_strict.CDS.tree \
   TESTER/REFERENCE_TESTER_TREES/iqtree.SCO_strict.codon_aln.tree

# Update test expectations if tree filenames change
```

Document the regeneration in a commit message and verify the new reference trees are
scientifically correct (manual inspection / comparison to known phylogenies).

### Handling nondeterminism
If a test fails due to solver nondeterminism (e.g., IQ-TREE finds a slightly different
tree on different runs), options:
1. **Increase tolerance** — change `tolerance=0` to `tolerance=2` in `same_topology()`.
2. **Fix the seed** — add `--seed 1234` to IQ-TREE invocations (already done in
   `control_file.user`).
3. **Swap the dataset** — use a dataset with stronger phylogenetic signal where the
   optimal tree is unambiguous.

The user has indicated they'll swap data if fragility arises, so start with `tolerance=0`
and adjust only if needed.

### Adding more test datasets
To test additional clades or edge cases (e.g., highly divergent taxa, missing data):
1. Add a new `TESTER/genomes_<name>/` directory with genomes + annotations.
2. Create a new test module `test_orthophyl_<name>.py` with a separate session fixture.
3. Generate reference trees for that dataset and add topology assertions.

---

## 7. Cross-reference with unit tests

| What | Unit tests (mocked) | Integration tests (real tools) |
|------|---------------------|--------------------------------|
| Routing decisions | ✅ `test_router.py` | ❌ (out of scope) |
| Wrapper orchestration | ✅ `test_wrapper_batch.py` | ❌ (out of scope) |
| CLI flag shape | ⚠️ (mock accepts any argv) | ✅ `test_releaf_completes_successfully` |
| Tool integration | ❌ (tools never run) | ✅ all tests |
| Topology correctness | ❌ (no real trees) | ✅ `test_topology_matches_reference` |
| False-success detection | ⚠️ (can mock returncode) | ✅ (real exit codes) |
| Execution speed | Fast (~seconds) | Slow (~5–15 min) |
| CI frequency | Every commit | Nightly / on-demand |

**Complementary, not redundant.** Both layers are necessary for comprehensive coverage.

---

## 8. Known limitations & future work

1. **No wrapper end-to-end test yet** — the integration suite currently tests
   `OrthoPhyl.sh` and `ReLeaf.sh` directly, not via `orthophyl_pipeline_wrapper.v2.py`.
   A future test could invoke the wrapper in non-dry-run mode on a fixture database,
   asserting `pipeline_status.json` reflects real success/failure (this would have caught
   the false-success bug B1–B3 at the integration level).

2. **Single dataset** — only fasttest (orchid chloroplasts) is covered. Additional
   datasets (bacteria, larger genomes, edge cases) would strengthen confidence.

3. **No performance regression tracking** — tests assert correctness, not speed. Consider
   adding timing assertions or benchmarks if runtime regressions become a concern.

4. **Manual reference tree validation** — the reference trees are assumed correct. A
   one-time manual review against published phylogenies would add confidence.

---

## Summary

This integration test suite provides **end-to-end validation** of the OrthoPhyl and ReLeaf
scientific pipelines, catching bugs that unit tests (with mocked subprocess calls) cannot
detect. By comparing produced trees to reference topologies via Robinson-Foulds distance,
the tests assert not just that the pipeline runs, but that it produces **scientifically
correct results**.

Combined with the unit test suite in `WRAPPER_UNIT_TEST_PLAN.md`, this provides
comprehensive coverage of both orchestration logic and scientific correctness.
