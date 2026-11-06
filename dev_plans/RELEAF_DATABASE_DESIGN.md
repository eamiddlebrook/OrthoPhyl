# ReLeaf Database Design - Analysis & Implementation

## How ReLeaf Works

### Input Requirements
ReLeaf needs access to a **complete OrthoPhyl run directory** (`--store` argument) containing:

1. **HMM profiles** (for finding orthologs in new genomes):
   - Primary: `OG_alignmentsToHMM/hmms_final/*.hmm`
   - Fallback: `annots_prots.fixed/OrthoFinder/Results_ortho/MultipleSequenceAlignments/*.fa`

2. **Old alignments** (to add new sequences to):
   - Proteins: `phylo_current/AlignmentsProts/OG*.fa`
   - CDS: `phylo_current/AlignmentsCDS/OG*.fa`

3. **Trimmed column info** (to maintain alignment consistency):
   - `phylo_current/trimmed_columns/`

4. **Existing phylogenies** (to constrain new trees):
   - `FINAL_SPECIES_TREES/iqtree.SCO_strict.CDS.tree`
   - `FINAL_SPECIES_TREES/iqtree.SCO_strict.PROT.tree`
   - `FINAL_SPECIES_TREES/fastTree.*.tree` (if used)

5. **IQTree partition schemes** (if used):
   - `phylo_current/SpeciesTree/iqtree/iqtree.SCO_strict.CDS.best_scheme`

6. **SCO (Single Copy Ortholog) sets**:
   - `phylo_current/SCO_strict`
   - `phylo_current/SCO_[min_frac]` (if relaxed mode used)

### ReLeaf Workflow

```
Input: New genome assemblies
    ↓
1. Annotate genomes (Prodigal)
    ↓
2. HMM search against pre-computed OG profiles
   → Identifies orthologs in new genomes
    ↓
3. Add new sequences to existing alignments (MAFFT --add --keeplength)
   → Maintains alignment structure
    ↓
4. Apply old trimming scheme to new alignments
   → Ensures columns match original
    ↓
5. Concatenate alignments using original SCO sets
    ↓
6. Build constrained ML tree (IQTree -g or FastTree with constraints)
   → New sequences added to existing tree topology
    ↓
Output: Updated phylogeny with new genomes placed
```

---

## Database Design Problem

### Current Situation
- `create_hierarchical_database_v2.py` creates a symlink: `database_name_db/orthophyl_run → /path/to/original/orthophyl/output`
- ReLeaf expects `--store /path/to/orthophyl/run` which should contain all the files listed above
- **Problem**: After running ReLeaf, the output goes to `$store/ReLeaf_dir/`, but this creates a **version control issue**

### Version Control Challenge

```
Scenario:
1. Initial OrthoPhyl run with 100 genomes → Create database v1
2. ReLeaf adds 10 new genomes → outputs to ReLeaf_dir/
3. Want to use this as new database → Need database v2
4. ReLeaf adds 5 more genomes → Need database v3
5. Want to run ReLeaf on v1 AND v3 simultaneously
```

**Requirements:**
- Each ReLeaf run creates a new "version" of the database
- Each version should be independently usable
- Need to avoid file duplication (use symlinks)
- Must be clear which version is being used

---

## Proposed Solution: Database Versions

### Directory Structure

```
database_dir/
├── Rhizobiaceae_db/
│   ├── v1_orthophyl_initial/           # Original OrthoPhyl run
│   │   ├── database_config.json
│   │   ├── version_info.json
│   │   ├── orthophyl_run@ → /original/path
│   │   ├── phylogeny.nwk
│   │   └── genome_list.txt
│   │
│   ├── v2_releaf_2024-01-15/           # After first ReLeaf
│   │   ├── database_config.json
│   │   ├── version_info.json
│   │   ├── base_version@ → ../v1_orthophyl_initial
│   │   ├── releaf_output@ → /path/to/ReLeaf_dir
│   │   ├── orthophyl_run_composite/    # Synthesized structure
│   │   │   ├── OG_alignmentsToHMM/
│   │   │   │   └── hmms_final@ → /original/path/.../hmms_final
│   │   │   ├── phylo_current/
│   │   │   │   ├── AlignmentsProts@ → /ReLeaf_dir/new_prot_alignments.trm.nm
│   │   │   │   ├── AlignmentsCDS@ → /ReLeaf_dir/new_CDS_alignments.trm.nm
│   │   │   │   ├── trimmed_columns@ → /original/.../trimmed_columns
│   │   │   │   └── SpeciesTree@ → /original/.../SpeciesTree
│   │   │   └── FINAL_SPECIES_TREES/
│   │   │       ├── iqtree.SCO_strict.CDS.tree@ → /ReLeaf_dir/new_trees/...
│   │   │       └── ...
│   │   ├── phylogeny.nwk@ → releaf_output/new_trees/iqtree.SCO_strict.CDS.tree
│   │   └── genome_list.txt              # Updated list
│   │
│   ├── v3_releaf_2024-01-20/           # After second ReLeaf
│   │   └── ...
│   │
│   ├── current@ → v2_releaf_2024-01-15  # Symlink to latest version
│   └── database_config.json@ → current/database_config.json
```

### Version Info JSON

```json
{
  "version": "v2_releaf_2024-01-15",
  "version_number": 2,
  "version_type": "releaf",
  "created": "2024-01-15T10:30:00",
  "parent_version": "v1_orthophyl_initial",
  "base_genomes": 100,
  "added_genomes": 10,
  "total_genomes": 110,
  "added_genome_ids": ["GCF_001", "GCF_002", ...],
  "releaf_run": "/path/to/store/ReLeaf_dir",
  "tree_methods": ["iqtree"],
  "data_types": ["CDS", "PROT"],
  "notes": "Added 10 Rhizobium genomes from NCBI"
}
```

---

## Implementation: create_hierarchical_database_v3.py

### Key Features

1. **Version-aware database creation**
   - Initial: Creates v1 from OrthoPhyl run
   - Updates: Creates new version from ReLeaf output

2. **Composite directory synthesis**
   - Creates `orthophyl_run_composite/` with correct structure
   - Uses symlinks to avoid duplication
   - Points to most recent versions of each file type

3. **Version tracking**
   - Each version has `version_info.json`
   - `current` symlink points to latest
   - Can reference specific versions

4. **Backward compatibility**
   - Old scripts can use `database_name_db/orthophyl_run`
   - Points to composite directory of current version

### Usage

```bash
# Initial creation from OrthoPhyl runs
python create_hierarchical_database_v3.py \
    --input orthophyl_runs.tsv \
    --output-dir databases/ \
    --mode initial

# Add ReLeaf run as new version
python create_hierarchical_database_v3.py \
    --database-dir databases/Rhizobiaceae_db/ \
    --releaf-output /path/to/store/ReLeaf_dir/ \
    --mode releaf-update

# List versions
python create_hierarchical_database_v3.py \
    --database-dir databases/Rhizobiaceae_db/ \
    --mode list-versions

# Set active version
python create_hierarchical_database_v3.py \
    --database-dir databases/Rhizobiaceae_db/ \
    --set-current v1_orthophyl_initial
```

---

## ReLeaf Output Structure

### What ReLeaf Creates

```
$store/ReLeaf_dir/
├── new_prot_alignments.trm.nm/    # Updated protein alignments
│   └── OG*.new_aligned.fa
├── new_CDS_alignments.trm.nm/     # Updated CDS alignments
│   └── OG*.new_aligned.fa
├── new_trees/                      # New phylogenies
│   ├── iqtree.SCO_strict.CDS.addasm.treefile
│   └── iqtree.SCO_strict.PROT.addasm.treefile
├── genome_list                     # File listing all genomes
├── annots_prots/                   # New genome annotations
├── annots_CDS/
└── logs/
```

### What to Preserve for Next ReLeaf

1. **HMM profiles** → Keep from original (unchanged)
2. **Alignments** → Use NEW versions from ReLeaf_dir
3. **Trimmed columns** → Keep from original (unchanged)
4. **Trees** → Use NEW trees from ReLeaf_dir
5. **Partition schemes** → May need update (check IQTree output)
6. **SCO sets** → Keep from original (should be same)

---

## Composite Directory Construction

### For v2 (after first ReLeaf)

```python
def create_composite_orthophyl_dir(base_version_dir, releaf_output_dir, composite_dir):
    """
    Create a composite OrthoPhyl structure that ReLeaf can use.
    
    Combines:
    - Unchanged files from base version (HMMs, trimmed_columns, SCO sets)
    - Updated files from ReLeaf output (alignments, trees)
    """
    
    # Structure needed by ReLeaf
    composite_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. HMMs - from base version (unchanged)
    (composite_dir / "OG_alignmentsToHMM" / "hmms_final").mkdir(parents=True)
    symlink(
        base_version_dir / "orthophyl_run" / "OG_alignmentsToHMM" / "hmms_final",
        composite_dir / "OG_alignmentsToHMM" / "hmms_final"
    )
    
    # 2. Alignments - from ReLeaf output (UPDATED)
    (composite_dir / "phylo_current").mkdir(parents=True)
    symlink(
        releaf_output_dir / "new_prot_alignments.trm.nm",
        composite_dir / "phylo_current" / "AlignmentsProts.trm"
    )
    symlink(
        releaf_output_dir / "new_CDS_alignments.trm.nm",
        composite_dir / "phylo_current" / "AlignmentsCDS.trm"
    )
    
    # 3. Trimmed columns - from base version (unchanged)
    symlink(
        base_version_dir / "orthophyl_run" / "phylo_current" / "trimmed_columns",
        composite_dir / "phylo_current" / "trimmed_columns"
    )
    
    # 4. Trees - from ReLeaf output (UPDATED)
    (composite_dir / "FINAL_SPECIES_TREES").mkdir(parents=True)
    for tree_file in releaf_output_dir.glob("new_trees/*.tree*"):
        # Rename to match expected pattern
        new_name = tree_file.name.replace(".addasm", "")
        symlink(tree_file, composite_dir / "FINAL_SPECIES_TREES" / new_name)
    
    # 5. IQTree info - from base version (may need update logic)
    (composite_dir / "phylo_current" / "SpeciesTree" / "iqtree").mkdir(parents=True)
    symlink(
        base_version_dir / "orthophyl_run" / "phylo_current" / "SpeciesTree" / "iqtree",
        composite_dir / "phylo_current" / "SpeciesTree" / "iqtree"
    )
    
    # 6. SCO sets - from base version (unchanged)
    for sco_file in (base_version_dir / "orthophyl_run" / "phylo_current").glob("SCO_*"):
        symlink(sco_file, composite_dir / "phylo_current" / sco_file.name)
```

---

## Workflow Integration

### Updated Wrapper Workflow

```
1. Route assemblies
    ↓
2a. ReLeaf route (matched database)
    → Run ReLeaf with database version (default: current)
    → Create new database version from ReLeaf output
    → Update 'current' symlink
    
2b. OrthoPhyl route (novel taxa)
    → Download genomes + add queries
    → Run OrthoPhyl
    → Create v1 database
    
3. Aggregate results
```

### Command Updates

```bash
# ReLeaf with specific version
./ReLeaf.sh \
    --store databases/Rhizobiaceae_db/v1_orthophyl_initial/orthophyl_run \
    --input_genomes input/

# ReLeaf with current version (default)
./ReLeaf.sh \
    --store databases/Rhizobiaceae_db/current/orthophyl_run_composite \
    --input_genomes input/

# After ReLeaf, create new version
python create_hierarchical_database_v3.py \
    --database-dir databases/Rhizobiaceae_db/ \
    --releaf-output /path/to/ReLeaf_dir/ \
    --mode releaf-update
```

---

## Benefits

1. **Version control**: Clear provenance of each database state
2. **Reproducibility**: Can go back to any version
3. **Efficiency**: Symlinks avoid file duplication
4. **Flexibility**: Run ReLeaf on different versions simultaneously
5. **Automation**: Wrapper can automatically create new versions
6. **Compatibility**: Works with existing ReLeaf code

---

## Next Steps

1. Implement `create_hierarchical_database_v3.py`
2. Add version management functions
3. Update wrapper to:
   - Create new version after each ReLeaf run
   - Track which version was used for routing
4. Add validation to ensure composite directories have all required files
5. Implement version listing/comparison tools
```