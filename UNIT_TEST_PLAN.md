# OrthoPhyl Unit Testing Plan

## Overview
This document outlines a comprehensive unit testing strategy for the OrthoPhyl v2.2.1 pipeline and its ReLeaf extension. The testing plan is organized by component type and covers both Python scripts and Bash functions.

---

## 1. Python Scripts Testing

### 1.1 ANI_genome_picking.py

**Purpose**: Clusters genomes by Average Nucleotide Identity (ANI) and selects a diverse subset for OrthoFinder analysis.

**Test Suite**: `tests/test_ani_genome_picking.py`

#### Unit Tests:
1. **test_parse_ani_output**
   - Input: Sample fastANI output file
   - Expected: Correctly populated full_dict with ANI values
   - Edge cases: Empty file, malformed lines, special characters in genome names

2. **test_divergent_genome_handling**
   - Input: ANI values including comparisons that failed (set to 50)
   - Expected: Proper handling of divergent genomes (>25% identity cutoff)
   - Verify: Count of too-divergent comparisons is accurate

3. **test_clustering_algorithm**
   - Input: Full ANI distance matrix
   - Expected: Correct hierarchical clustering based on maximum ANI pairs
   - Verify: Number of clusters equals requested n_clusters

4. **test_cluster_merge**
   - Input: Two genomes with known ANI distances to others
   - Expected: Merged cluster has averaged distances to remaining genomes
   - Verify: Distance calculation formula: (sum/2) + (100-max_val)/2

5. **test_representative_selection**
   - Input: Clustered genome groups
   - Expected: One representative selected per cluster (highest average ANI)
   - Verify: Output file "Species_shortlist" contains correct number of genomes

6. **test_special_character_handling**
   - Input: Genome names with parentheses, quotes, tabs
   - Expected: Characters properly stripped/normalized
   - Verify: No parsing errors, consistent genome identification

#### Integration Tests:
1. **test_full_workflow**
   - Input: Complete fastANI output + genome list + cluster count
   - Expected: Valid Species_shortlist file with expected genome count
   - Performance: Run time <1 min for 100 genomes

2. **test_with_tester_data**
   - Input: TEST or TESTER_chloroplast dataset
   - Expected: Matches expected shortlist from known-good runs

---

### 1.2 OG_sco_filter.py

**Purpose**: Filters orthogroups to identify strict single-copy orthologs (SCOs) based on taxon representation threshold.

**Test Suite**: `tests/test_og_sco_filter.py`

#### Unit Tests:
1. **test_parse_orthogroup_counts**
   - Input: OrthoFinder Orthogroups.GeneCount.tsv file
   - Expected: Correctly parsed orthogroup counts per taxon
   - Edge cases: Header handling, empty fields, trailing commas

2. **test_single_copy_detection**
   - Input: Orthogroup with 0 or 1 copy per taxon
   - Expected: Identified as single-copy ortholog
   - Negative cases: Orthogroups with paralogs (>1 copy) rejected

3. **test_taxon_representation_threshold**
   - Input: Various homologue counts (e.g., 10, 50, 100)
   - Expected: Only OGs present in ≥threshold taxa are kept
   - Verify: Correct filtering based on total column count

4. **test_paralog_exclusion**
   - Input: Orthogroup with mixed 0/1/2+ counts
   - Expected: Any OG with >1 copy in any taxon is excluded
   - Verify: one_zero list [0,1] logic works correctly

5. **test_output_file_format**
   - Input: Filtered orthogroup list
   - Expected: Output file "SCO_{threshold}" with one OG per line
   - Verify: No whitespace issues, correct line endings

#### Integration Tests:
1. **test_strict_sco_filtering**
   - Input: Real OrthoFinder output with homologues=total_taxa
   - Expected: Only strict SCOs (present in 100% taxa) returned

2. **test_relaxed_sco_filtering**
   - Input: Real OrthoFinder output with homologues=0.3*total_taxa
   - Expected: Relaxed SCOs (present in ≥30% taxa) returned

---

### 1.3 Newick2FastTreeConstraints.py

**Purpose**: Converts a Newick tree into FastTree constraint format, filtering by support values.

**Test Suite**: `tests/test_newick_constraints.py`

#### Unit Tests:
1. **test_tree_parsing**
   - Input: Valid Newick format tree string
   - Expected: ete3 Tree object created successfully
   - Edge cases: Malformed Newick, missing branch lengths

2. **test_support_value_filtering**
   - Input: Tree with mixed support values (0.5-1.0)
   - Expected: Nodes with support <0.95 are deleted
   - Verify: Tree topology changes correctly

3. **test_binary_constraint_generation**
   - Input: Tree with N leaves and M internal nodes
   - Expected: M binary constraint strings (0s and 1s)
   - Verify: Each leaf has correct membership in each clade

4. **test_postorder_traversal**
   - Input: Sample tree structure
   - Expected: Constraints generated in postorder
   - Verify: Child clades processed before parent clades

5. **test_fasta_output_format**
   - Input: Constraint dictionary
   - Expected: FASTA format with >name headers
   - Verify: Binary strings match tree topology

6. **test_leaf_name_handling**
   - Input: Leaves with special characters
   - Expected: Names preserved in output
   - Edge cases: Spaces, punctuation, long names

#### Integration Tests:
1. **test_full_conversion**
   - Input: Real tree from TESTER dataset
   - Expected: Valid constraint file usable by FastTree
   - Verify: Number of constraints equals number of well-supported nodes

2. **test_roundtrip_consistency**
   - Input: Tree → constraints → FastTree with constraints
   - Expected: Resulting tree respects original topology

---

## 2. Bash Function Testing

### 2.1 Core Pipeline Functions (functions.sh)

**Test Suite**: `tests/test_bash_functions.bats` (using BATS framework)

#### 2.1.1 SET_UP_DIR_STRUCTURE

**Tests:**
1. **test_directory_creation**
   - Verify all required directories are created
   - Check: AlignmentsProts, AlignmentsCDS, logs, etc.

2. **test_idempotency**
   - Run function twice with setup.complete present
   - Expected: No duplicate work, proper detection

3. **test_min_frac_orthos_calculation**
   - Input: Various genome counts and fractions
   - Expected: Correct integer calculation of threshold

4. **test_genome_list_aggregation**
   - Input: Genomes + pre-annotated samples
   - Expected: Correct all_input_list file

#### 2.1.2 CLEAN_N_COPY_GENOMES

**Tests:**
1. **test_special_character_removal**
   - Input: Filenames with parentheses, spaces, etc.
   - Expected: Clean filenames in output directory

2. **test_symlink_creation**
   - Verify soft links created correctly
   - Check no data duplication

3. **test_empty_directory_handling**
   - Input: Empty genome directory
   - Expected: Graceful error or warning

#### 2.1.3 PRODIGAL_PREDICT

**Tests:**
1. **test_annotation_output**
   - Input: Sample bacterial genome
   - Expected: CDS and protein FASTA files generated

2. **test_parallel_execution**
   - Verify parallel processes respect thread limit
   - Check: No race conditions in file writing

3. **test_partial_failure_handling**
   - Input: Mix of valid and corrupted FASTA files
   - Expected: Valid genomes processed, errors logged

#### 2.1.4 CHECK_ANNOTS

**Tests:**
1. **test_duplicate_removal**
   - Input: Proteins with 100% identity
   - Expected: Duplicates removed with dedupe.sh

2. **test_sequence_validation**
   - Check for: valid FASTA format, protein alphabet
   - Expected: Invalid sequences flagged

#### 2.1.5 ORTHO_RUN

**Tests:**
1. **test_orthofinder_execution**
   - Mock OrthoFinder call
   - Verify: Correct parameters passed

2. **test_output_parsing**
   - Expected: Orthogroups.txt and related files exist
   - Verify: Parsing doesn't fail on empty OGs

#### 2.1.6 ANI_ORTHOFINDER_TO_ALL_SEQS

**Tests:**
1. **test_hmm_building**
   - Input: Protein alignment
   - Expected: Valid HMM file created with hmmbuild

2. **test_hmmsearch_execution**
   - Input: HMM + target proteomes
   - Expected: Hits identified and filtered by e-value

3. **test_iterative_expansion**
   - Verify: Round 1 and Round 2 HMMs created
   - Check: New sequences added to orthogroups

4. **test_parallel_hmm_searches**
   - Performance: Thread utilization efficient
   - Correctness: No output file corruption

#### 2.1.7 REALIGN_ORTHOGROUP_PROTS

**Tests:**
1. **test_mafft_alignment**
   - Input: Unaligned orthogroup sequences
   - Expected: Aligned FASTA with gaps

2. **test_alignment_quality_check**
   - Verify: No sequences dropped unexpectedly
   - Check: Alignment length consistency

#### 2.1.8 GET_OG_NAMES

**Tests:**
1. **test_orthogroup_listing**
   - Expected: Complete list of OG identifiers
   - Verify: No duplicates, proper formatting

#### 2.1.9 TRIM

**Tests:**
1. **test_trimal_backtranslation**
   - Input: Protein alignment + CDS sequences
   - Expected: Codon-aligned CDS with gaps in correct positions

2. **test_gap_threshold_filtering**
   - Input: Various -gt values
   - Expected: Columns removed per threshold

3. **test_segfault_handling**
   - Input: Alignment with all sites removed
   - Expected: Graceful failure, logged error

#### 2.1.10 GET_TRIMMED_COLS

**Tests:**
1. **test_column_tracking**
   - Verify: Accurate count of removed alignment columns
   - Output: Correct summary statistics

#### 2.1.11 ALIGNMENT_STATS

**Tests:**
1. **test_stat_calculation**
   - Input: Trimmed alignment
   - Expected: Length, gap%, GC content, etc.

2. **test_r_script_integration**
   - Verify: R scripts execute without error
   - Check: Output files parseable

#### 2.1.12 SCO_MIN_ALIGN

**Tests:**
1. **test_sco_filtering**
   - Input: Various min_frac_orthos values
   - Expected: Correct OGs kept/excluded

2. **test_concatenation**
   - Tool: catfasta2phyml
   - Expected: Supermatrix with proper partitioning

3. **test_matrix_stats**
   - Verify: Missing data calculations accurate
   - Check: Per-taxon coverage reporting

#### 2.1.13 TREE_BUILD

**Tests:**
1. **test_fasttree_execution**
   - Input: Concatenated alignment
   - Expected: Newick tree file

2. **test_raxml_execution**
   - Verify: RAxML parameters correct
   - Check: Bootstrap support values present

3. **test_iqtree_execution**
   - Test: Model selection, partition merging
   - Verify: Best scheme file generated

4. **test_astral_execution**
   - Input: Gene trees
   - Expected: Species tree with local posterior probabilities

5. **test_constraint_tree_building**
   - Input: Previous tree + new samples
   - Expected: Constraint respected in new tree

#### 2.1.14 WRAP_UP

**Tests:**
1. **test_file_organization**
   - Verify: Final trees moved to FINAL_SPECIES_TREES
   - Check: Summary files generated

2. **test_completion_marker**
   - Expected: Success indicators written

---

### 2.2 ReLeaf Functions (functions_addem.sh)

**Test Suite**: `tests/test_releaf_functions.bats`

#### 2.2.1 SET_UP_ADDASM_DIRS

**Tests:**
1. **test_directory_setup**
   - Verify: ReLeaf-specific directories created
   - Check: No conflict with original OP run

#### 2.2.2 GET_OLD_TREE_INFO

**Tests:**
1. **test_hmm_directory_detection**
   - Input: Completed OP storage directory
   - Expected: HMM directory path identified

2. **test_tree_parsing**
   - Expected: Extract tree methods and omics types used
   - Verify: Parameters match original run

3. **test_model_extraction**
   - Input: IQTree best_scheme file
   - Expected: Partition models extracted for reuse

#### 2.2.3 ANNOTATIONS

**Tests:**
1. **test_new_sample_annotation**
   - Input: New genome assemblies
   - Expected: CDS and protein files generated

2. **test_pre_annotated_handling**
   - Input: User-provided CDS/protein files
   - Expected: Files validated and integrated

#### 2.2.4 ADD_2_ALIGNMENTS

**Tests:**
1. **test_hmmsearch_on_new_samples**
   - Input: Old HMMs + new proteomes
   - Expected: Orthogroup membership assigned

2. **test_alignment_guide_usage**
   - Input: Old alignment as guide
   - Expected: New sequences aligned consistently

3. **test_sequence_addition**
   - Verify: New sequences appended to existing alignments
   - Check: No disruption to existing sequences

#### 2.2.5 ADD_2_IQTREE

**Tests:**
1. **test_model_reuse**
   - Input: Original IQTree partition models
   - Expected: Same models applied to expanded alignment

2. **test_constraint_tree_application**
   - Input: Original tree as constraint
   - Expected: New tree extends original topology

#### 2.2.6 RUN_ML_SUPERMATRIX_WORKFLOW

**Tests:**
1. **test_full_releaf_ml_pipeline**
   - Input: New samples + old run directory
   - Expected: Updated tree with all samples

2. **test_runtime_comparison**
   - Verify: ReLeaf significantly faster than full rerun
   - Performance: <20% of original runtime for 10% new samples

#### 2.2.7 ASTRAL_WORKFLOW

**Tests:**
1. **test_gene_tree_generation**
   - Input: Expanded gene alignments
   - Expected: Updated gene trees

2. **test_astral_tree_building**
   - Note: Currently not using constraint trees
   - Expected: Species tree from all gene trees

---

## 3. Integration Testing

### 3.1 Full Pipeline Tests

**Test Suite**: `tests/integration/test_orthophyl_pipeline.sh`

#### Tests:
1. **test_tester_fasttest**
   - Runtime: ~8 minutes on 3 cores
   - Expected: 4 species trees in FINAL_SPECIES_TREES
   - Data: TESTER/genomes_fasttest

2. **test_tester_chloroplast**
   - Runtime: ~8 minutes on 3 cores
   - Expected: 4 species trees
   - Data: TESTER/genomes_chloroplast

3. **test_full_tester**
   - Runtime: ~20 hours on 20 cores
   - Expected: Complete workflow with all tree methods
   - Data: TESTER/genomes

4. **test_pre_annotated_input**
   - Input: Mix of genomes + CDS/protein files
   - Expected: Integrated analysis

5. **test_releaf_addition**
   - Step 1: Run TESTER_fasttest
   - Step 2: Add samples with ReLeaf
   - Expected: Updated trees with new samples

### 3.2 Parameter Variation Tests

**Test Suite**: `tests/integration/test_parameter_variations.sh`

#### Tests:
1. **test_rigor_levels**
   - Run with: fast, medium, full
   - Verify: Appropriate tools called

2. **test_tree_methods**
   - Separately test: fasttree, raxml, iqtree, astral
   - Verify: Each produces valid output

3. **test_omics_types**
   - Run with: CDS, PROT, BOTH
   - Verify: Correct molecular data used

4. **test_ani_thresholds**
   - Vary: max genomes for OrthoFinder (10, 20, 50)
   - Verify: Correct subset selected

5. **test_min_frac_orthos**
   - Values: 0.3, 0.5, 0.7, 1.0
   - Verify: SCO filtering thresholds correct

### 3.3 Error Handling Tests

**Test Suite**: `tests/integration/test_error_handling.sh`

#### Tests:
1. **test_empty_genome_directory**
2. **test_corrupted_fasta_files**
3. **test_insufficient_disk_space**
4. **test_missing_dependencies**
5. **test_keyboard_interrupt**
6. **test_partial_run_recovery**

---

## 4. Test Data Organization

### 4.1 Directory Structure
```
tests/
├── unit/
│   ├── python/
│   │   ├── test_ani_genome_picking.py
│   │   ├── test_og_sco_filter.py
│   │   └── test_newick_constraints.py
│   └── bash/
│       ├── test_bash_functions.bats
│       └── test_releaf_functions.bats
├── integration/
│   ├── test_orthophyl_pipeline.sh
│   ├── test_parameter_variations.sh
│   └── test_error_handling.sh
├── fixtures/
│   ├── sample_ani_output.txt
│   ├── sample_orthogroups.tsv
│   ├── sample_tree.nwk
│   ├── sample_genomes/
│   └── expected_outputs/
└── test_utils.sh
```

### 4.2 Test Fixtures
- **Minimal test genome set**: 5-10 small bacterial genomes (<1MB each)
- **Synthetic ANI matrix**: Known cluster structure
- **Reference orthogroup counts**: Pre-computed for validation
- **Expected trees**: Known-good outputs for regression testing

---

## 5. Testing Framework Setup

### 5.1 Python Testing
**Framework**: pytest + pytest-cov

**Installation**:
```bash
pip install pytest pytest-cov pytest-mock
```

**Run commands**:
```bash
# All Python tests
pytest tests/unit/python/

# With coverage
pytest --cov=python_scripts --cov-report=html tests/unit/python/

# Specific test
pytest tests/unit/python/test_ani_genome_picking.py::test_clustering_algorithm
```

### 5.2 Bash Testing
**Framework**: BATS (Bash Automated Testing System)

**Installation**:
```bash
git clone https://github.com/bats-core/bats-core.git
cd bats-core
./install.sh /usr/local
```

**Run commands**:
```bash
# All bash tests
bats tests/unit/bash/*.bats

# Specific test file
bats tests/unit/bash/test_bash_functions.bats
```

### 5.3 Integration Testing
**Framework**: Custom bash scripts

**Run commands**:
```bash
# Quick integration test
bash tests/integration/test_orthophyl_pipeline.sh fasttest

# Full integration suite
bash tests/integration/run_all_integration_tests.sh
```

---

## 6. Continuous Integration Setup

### 6.1 GitHub Actions Workflow
**File**: `.github/workflows/test.yml`

**Stages**:
1. **Lint**: shellcheck, pylint
2. **Unit tests**: Python + Bash
3. **Fast integration**: TESTER_fasttest only
4. **Nightly full test**: Complete TESTER suite

### 6.2 Test Environments
- **Python**: 3.7, 3.8, 3.9, 3.10
- **Bash**: 4.4, 5.0, 5.1
- **OS**: Ubuntu 20.04, 22.04

---

## 7. Coverage Goals

### 7.1 Target Coverage
- **Python scripts**: >90% line coverage
- **Bash functions**: >80% function coverage
- **Integration paths**: All major parameter combinations

### 7.2 Priority Areas
1. **Critical path**: Annotation → OrthoFinder → Alignment → Tree
2. **Error handling**: File I/O, dependency failures
3. **Edge cases**: Empty inputs, extreme parameters
4. **ReLeaf specific**: HMM search, constraint trees

---

## 8. Test Execution Order

### 8.1 Development Workflow
1. Run relevant unit tests during development
2. Run fast integration test before commit
3. Full integration test before PR

### 8.2 CI/CD Workflow
1. Lint and static analysis (1 min)
2. Python unit tests (5 min)
3. Bash unit tests (10 min)
4. Fast integration test (10 min)
5. Nightly: Full integration suite (1-2 hours)

---

## 9. Mocking Strategy

### 9.1 External Dependencies
**Mock**:
- OrthoFinder (use pre-computed results)
- FastTree/RAxML/IQTree (use small test cases or cached trees)
- HMMER (use truncated databases)
- MAFFT (use pre-aligned fixtures)

**Reason**: Reduce test time from hours to minutes

### 9.2 File System
- Use temporary directories for all tests
- Clean up after each test
- Use fixtures for input data (read-only)

---

## 10. Performance Testing

### 10.1 Benchmarks
**Test Suite**: `tests/performance/benchmark_*.sh`

**Metrics**:
1. **Scalability**: Runtime vs. number of genomes (10, 50, 100, 500)
2. **Thread efficiency**: Speedup vs. number of cores
3. **Memory usage**: Peak RAM per stage
4. **Disk I/O**: Read/write patterns

### 10.2 Regression Detection
- Track runtime of standard tests
- Alert if >10% slower than baseline
- Compare memory usage between versions

---

## 11. Documentation Tests

### 11.1 README Examples
**Test**: All code examples in README.md are executable and produce expected output

**Test Suite**: `tests/doc/test_readme_examples.sh`

### 11.2 Help Messages
**Test**: All scripts produce valid help output with `-h` flag

---

## 12. Maintenance and Evolution

### 12.1 Test Maintenance
- Update tests when pipeline changes
- Add regression tests for every bug fix
- Retire obsolete tests

### 12.2 Test Review
- Quarterly review of test coverage
- Annual review of integration test suite
- Remove redundant tests

---

## 13. Success Criteria

A test suite is successful if:
1. ✅ All critical pipeline paths have tests
2. ✅ Tests catch regressions before production
3. ✅ Tests run fast enough to not block development (<15 min for unit tests)
4. ✅ Test failures provide actionable error messages
5. ✅ Coverage metrics meet targets
6. ✅ Integration tests validate real-world usage

---

## Appendix A: Quick Start Testing Guide

### Running Your First Test
```bash
# 1. Install test dependencies
pip install pytest pytest-cov
git clone https://github.com/bats-core/bats-core.git
cd bats-core && ./install.sh /usr/local && cd ..

# 2. Run Python unit tests
pytest tests/unit/python/test_ani_genome_picking.py -v

# 3. Run Bash unit tests
bats tests/unit/bash/test_bash_functions.bats

# 4. Run quick integration test
bash tests/integration/test_orthophyl_pipeline.sh fasttest

# 5. Check coverage
pytest --cov=python_scripts --cov-report=term tests/unit/python/
```

### Writing Your First Test
See `tests/examples/example_test.py` and `tests/examples/example_test.bats` for templates.

---

## Appendix B: Known Testing Challenges

1. **Long runtimes**: Full pipeline takes 20+ hours
   - Solution: Use smaller test datasets, mock expensive steps

2. **External dependencies**: Many bioinformatics tools required
   - Solution: Docker/Singularity containers for testing

3. **Stochastic outputs**: Some tools have randomness
   - Solution: Test for properties (valid tree, correct format) not exact matches

4. **Large data**: Real genomes are 2-5 MB each
   - Solution: Synthetic or minimal genomes for unit tests

5. **Platform differences**: Conda, OS variations
   - Solution: Matrix testing in CI

---

## Appendix C: Test Implementation Priority

### Phase 1 (Weeks 1-2): Foundation
1. Python unit tests for all 3 scripts
2. Basic bash function tests (5 most critical)
3. One fast integration test (TESTER_fasttest)

### Phase 2 (Weeks 3-4): Expansion
1. Complete bash function test coverage
2. ReLeaf-specific tests
3. Parameter variation integration tests

### Phase 3 (Weeks 5-6): Polish
1. Error handling tests
2. Performance benchmarks
3. CI/CD setup

### Phase 4 (Ongoing): Maintenance
1. Add tests for bug fixes
2. Update tests for new features
3. Monitor and improve coverage
