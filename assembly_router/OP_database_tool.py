#!/usr/bin/env python3
"""
Create Hierarchical Taxonomy Database from OrthoPhyl Runs (v2 - with update support)

This script creates a queryable database from one or more OrthoPhyl output directories.
Supports incremental updates - only processes new entries when run again.

Input Format (TSV):
    clade_name<TAB>orthophyl_dir<TAB>clade_taxonomy
    
    Example:
    Gaiellales	/path/to/gaiellales_run	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales
    Escherichia	/path/to/ecoli_run	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia

Usage:
    # Initial creation
    python OP_database_tool.py \
        --input orthophyl_runs.tsv \
        --output-dir taxonomy_databases/
    
    # Update with new entries (skips existing)
    python OP_database_tool.py \
        --input orthophyl_runs.tsv \
        --output-dir taxonomy_databases/ \
        --update
    
    # Force rebuild all
    python OP_database_tool.py \
        --input orthophyl_runs.tsv \
        --output-dir taxonomy_databases/ \
        --force

    # Preview removing one or more clade databases (dry run, deletes nothing)
    python OP_database_tool.py \
        --output-dir taxonomy_databases/ \
        --remove Gaiellales Escherichia

    # Actually remove them and refresh the master index
    python OP_database_tool.py \
        --output-dir taxonomy_databases/ \
        --remove Gaiellales Escherichia \
        --yes

Author: Generated for HGTool
Date: 2024
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import logging
from datetime import datetime
import shutil

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def relative_symlink(link: Path, target: Path) -> None:
    """Create `link` -> `target` using a relative path.

    A relative symlink stays correct under any bind-mount mapping (e.g.
    Singularity `--bind host:/container`), as long as `link` and `target`
    are on the same mount -- unlike `target.resolve()`, which bakes in
    whatever path the current process sees (the container's view), breaking
    on the host or under a different container bind later.
    """
    target = target.resolve()
    rel = os.path.relpath(target, link.parent.resolve())
    if os.path.commonpath([target, link.parent.resolve()]) == os.sep:
        logger.warning(
            f"  ⚠ {link.name} and its target share no common path above '/' "
            f"-- this symlink will only resolve where both paths are visible "
            f"the same way (e.g. inside the same container bind)."
        )
    link.symlink_to(rel)


def db_dir_for(output_dir: Path, clade_name: str, parent_taxon: Optional[str] = None) -> Path:
    """Compute the on-disk directory for a database, nested or flat.

    A megatree subclade/backbone DB (parent_taxon set) lives at
    <output_dir>/<parent_taxon>/<safe_name>_db so a taxon's backbone and all
    its subclades group together instead of sitting as unrelated-looking
    siblings next to every other top-level taxon DB. Everything else stays
    flat at <output_dir>/<safe_name>_db, unchanged from before nesting existed.
    """
    safe_name = clade_name.replace(' ', '_').replace('/', '_')
    if parent_taxon:
        safe_parent = parent_taxon.replace(' ', '_').replace('/', '_')
        return output_dir / safe_parent / f"{safe_name}_db"
    return output_dir / f"{safe_name}_db"


def iter_database_dirs(output_dir: Path):
    """Yield every *_db directory under output_dir, flat or one level nested.

    Looks at output_dir/*_db (flat databases) plus output_dir/<X>/*_db for
    any <X> that isn't hidden and doesn't itself end in "_db" (so a nested
    subclade dir is never mistaken for a parent-taxon grouping dir, and
    hidden dirs like .pipeline_runs are never descended into). Matches
    purely on naming/location -- it does NOT require database_config.json
    to exist, mirroring the plain output_dir.glob("*_db") this replaces
    (some callers, e.g. the wrapper's ReLeaf-tree collection, care about a
    DB dir's other contents even when its config is missing; callers that
    DO need the config already tolerate a missing/unreadable one via
    try/except around the open() call).
    """
    if not output_dir.exists():
        return
    yield from sorted(output_dir.glob("*_db"))
    for sub in sorted(output_dir.iterdir()):
        if not sub.is_dir() or sub.name.startswith('.') or sub.name.endswith('_db'):
            continue
        yield from sorted(sub.glob("*_db"))


def find_database_dir(output_dir: Path, clade_name: str) -> Optional[Path]:
    """Locate an existing database's directory, flat or nested.

    Tries the flat path first (the common case, and where non-megatree DBs
    always live), then falls back to a one-level scan for <parent>/<safe>_db.
    """
    safe_name = clade_name.replace(' ', '_').replace('/', '_')
    flat = output_dir / f"{safe_name}_db"
    if flat.exists() and (flat / "database_config.json").exists():
        return flat
    if not output_dir.exists():
        return None
    for sub in sorted(output_dir.iterdir()):
        if not sub.is_dir() or sub.name.startswith('.') or sub.name.endswith('_db'):
            continue
        candidate = sub / f"{safe_name}_db"
        if candidate.exists() and (candidate / "database_config.json").exists():
            return candidate
    return None


class GTDBTaxonomy:
    """Parse GTDB taxonomy strings."""
    
    RANKS = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    RANK_NAMES = {
        'd': 'domain',
        'p': 'phylum', 
        'c': 'class',
        'o': 'order',
        'f': 'family',
        'g': 'genus',
        's': 'species'
    }
    
    def __init__(self, taxonomy_string: str):
        self.raw = taxonomy_string.strip()
        self.levels = self._parse_taxonomy(self.raw)
    
    def _parse_taxonomy(self, taxonomy_string: str) -> Dict[str, str]:
        """Parse taxonomy into dictionary."""
        import re
        levels = {}
        
        parts = taxonomy_string.split(';')
        for part in parts:
            part = part.strip()
            if not part:
                continue
            
            match = re.match(r'^([dpcofgs])__(.*)$', part)
            if match:
                rank = match.group(1)
                name = match.group(2).strip()
                levels[rank] = name if name else None
        
        return levels
    
    def get_rank(self, rank: str) -> Optional[str]:
        """Get taxonomy at specific rank."""
        return self.levels.get(rank)
    
    def get_most_specific_rank(self) -> Optional[str]:
        """Get the most specific (lowest) defined rank."""
        for rank in reversed(self.RANKS):
            if self.levels.get(rank):
                return rank
        return None
    
    def get_rank_name(self, rank: str) -> str:
        """Get full name of rank."""
        return self.RANK_NAMES.get(rank, rank)
    
    def __str__(self):
        return self.raw


def validate_orthophyl_run(orthophyl_dir: Path) -> Dict:
    """
    Validate that OrthoPhyl directory has required outputs.
    Detects available tree methods (iqtree, fasttree) and data types (CDS, PROT).
    
    Returns:
        Dictionary with validation results and metadata
    """
    orthophyl_dir = Path(orthophyl_dir)
    
    if not orthophyl_dir.exists():
        raise FileNotFoundError(f"OrthoPhyl directory not found: {orthophyl_dir}")
    
    validation = {
        'valid': False,
        'has_hmms': False,
        'has_alignments': False,
        'has_trees': False,
        'n_genomes': 0,
        'hmm_dir': None,
        'tree_file': None,
        'warnings': [],
        'validation_details': [],
        # Track available tree methods and data types
        'available_trees': [],  # List of available tree files with metadata
        'tree_methods': set(),  # iqtree, fasttree, raxml
        'tree_data_types': set(),  # CDS, PROT
    }
    
    # Check for HMM profiles (required for ReLeaf)
    hmm_locations = [
        orthophyl_dir / "OG_alignmentsToHMM" / "hmms_final",
        orthophyl_dir / "annots_prots.fixed" / "OrthoFinder" / "Results_ortho" / "MultipleSequenceAlignments",
        orthophyl_dir / "phylo_current" / "AlignmentsProts"  # Also check for alignments here
    ]
    
    for hmm_dir in hmm_locations:
        if hmm_dir.exists():
            hmm_files = list(hmm_dir.glob("*.hmm")) + list(hmm_dir.glob("*.fa")) + list(hmm_dir.glob("*.faa"))
            if hmm_files:
                validation['has_hmms'] = True
                validation['hmm_dir'] = hmm_dir
                validation['validation_details'].append(f"Found {len(hmm_files)} HMM/alignment files in {hmm_dir.name}")
                logger.info(f"  ✓ Found {len(hmm_files)} files in {hmm_dir.relative_to(orthophyl_dir)}")
                break
    
    if not validation['has_hmms']:
        validation['warnings'].append("No HMM profiles or alignments found - ReLeaf may not work")
        validation['validation_details'].append("WARNING: No HMM profiles found")
        logger.warning("  ⚠ No HMM profiles found")
    
    # Check for alignments
    alignment_locations = [
        orthophyl_dir / "phylo_current" / "AlignmentsProts.trm",
        orthophyl_dir / "phylo_current" / "AlignmentsCDS.trm",
        orthophyl_dir / "phylo_current" / "AlignmentsProts",
        orthophyl_dir / "phylo_current" / "AlignmentsCDS"
    ]
    
    for aln_dir in alignment_locations:
        if aln_dir.exists():
            aln_files = list(aln_dir.glob("*.fa")) + list(aln_dir.glob("*.faa"))
            if aln_files:
                validation['has_alignments'] = True
                validation['validation_details'].append(f"Found {len(aln_files)} alignment files in {aln_dir.name}")
                logger.info(f"  ✓ Found {len(aln_files)} alignments in {aln_dir.relative_to(orthophyl_dir)}")
                break
    
    # Check for species trees - detect method and data type from filename
    tree_locations = [
        orthophyl_dir / "FINAL_SPECIES_TREES",
        orthophyl_dir / "phylo_current" / "SpeciesTree"
    ]
    
    for tree_dir in tree_locations:
        if tree_dir.exists():
            tree_files = list(tree_dir.glob("*.tree")) + list(tree_dir.glob("*.nwk")) + list(tree_dir.glob("*.treefile"))
            
            for tree_file in tree_files:
                tree_info = {
                    'file': tree_file,
                    'name': tree_file.name,
                    'method': None,
                    'data_type': None
                }
                
                # Detect tree method from filename
                name_lower = tree_file.name.lower()
                if 'iqtree' in name_lower:
                    tree_info['method'] = 'iqtree'
                    validation['tree_methods'].add('iqtree')
                elif 'fasttree' in name_lower:
                    tree_info['method'] = 'fasttree'
                    validation['tree_methods'].add('fasttree')
                elif 'raxml' in name_lower:
                    tree_info['method'] = 'raxml'
                    validation['tree_methods'].add('raxml')
                
                # Detect data type from filename
                if 'cds' in name_lower:
                    tree_info['data_type'] = 'CDS'
                    validation['tree_data_types'].add('CDS')
                elif 'prot' in name_lower or 'protein' in name_lower:
                    tree_info['data_type'] = 'PROT'
                    validation['tree_data_types'].add('PROT')
                
                validation['available_trees'].append(tree_info)
            
            if validation['available_trees']:
                validation['has_trees'] = True
                
                # Prefer IQ-TREE > RAxML > FastTree, and CDS > PROT
                priority_order = [
                    ('iqtree', 'CDS'),
                    ('iqtree', 'PROT'),
                    ('raxml', 'CDS'),
                    ('raxml', 'PROT'),
                    ('fasttree', 'CDS'),
                    ('fasttree', 'PROT'),
                ]
                
                for method, data_type in priority_order:
                    for tree in validation['available_trees']:
                        if tree['method'] == method and tree['data_type'] == data_type:
                            validation['tree_file'] = tree['file']
                            break
                    if validation['tree_file']:
                        break
                
                # Fallback to first tree if no match
                if not validation['tree_file']:
                    validation['tree_file'] = validation['available_trees'][0]['file']
                
                validation['validation_details'].append(
                    f"Found {len(validation['available_trees'])} tree(s): " +
                    ", ".join([t['name'] for t in validation['available_trees']])
                )
                
                logger.info(f"  ✓ Found {len(validation['available_trees'])} tree file(s)")
                for tree in validation['available_trees']:
                    logger.info(f"    - {tree['name']} (method={tree['method']}, data={tree['data_type']})")
                
                break
    
    if not validation['has_trees']:
        validation['warnings'].append("No species tree found")
        validation['validation_details'].append("WARNING: No species tree found")
        logger.warning("  ⚠ No species tree found")
    
    # Count genomes - try multiple sources
    genome_sources = [
        orthophyl_dir / "genome_list",
        orthophyl_dir / "all_input_list",
        orthophyl_dir / "store" / "genome_list",
        orthophyl_dir / "store" / "all_input_list"
    ]
    
    for genome_list in genome_sources:
        if genome_list.exists():
            with open(genome_list, 'r') as f:
                validation['n_genomes'] = len([l for l in f if l.strip() and not l.startswith('#')])
            if validation['n_genomes'] > 0:
                validation['validation_details'].append(f"Counted genomes from {genome_list.name}")
                logger.info(f"  ✓ Counted {validation['n_genomes']} genomes from {genome_list.name}")
                break
    
    # Fallback: count protein files
    if validation['n_genomes'] == 0:
        prot_dirs = [
            orthophyl_dir / "annots_prots",
            orthophyl_dir / "store" / "annots_prots"
        ]
        for prot_dir in prot_dirs:
            if prot_dir.exists():
                validation['n_genomes'] = len(list(prot_dir.glob("*.faa")))
                if validation['n_genomes'] > 0:
                    validation['validation_details'].append(f"Counted {validation['n_genomes']} .faa files")
                    logger.info(f"  ✓ Counted {validation['n_genomes']} protein files")
                    break
    
    if validation['n_genomes'] == 0:
        validation['warnings'].append("Could not determine number of genomes")
        validation['validation_details'].append("WARNING: Could not count genomes")
        logger.warning("  ⚠ Could not determine genome count")
    
    # Overall validation - accept if we have HMMs/alignments OR trees
    validation['valid'] = (validation['has_hmms'] or validation['has_alignments']) and validation['n_genomes'] > 0
    
    if not validation['valid']:
        logger.error("  ✗ Validation FAILED")
        logger.error(f"    Has HMMs/alignments: {validation['has_hmms'] or validation['has_alignments']}")
        logger.error(f"    Genome count: {validation['n_genomes']}")
    
    return validation


def database_exists(clade_name: str, output_dir: Path) -> Optional[Path]:
    """
    Check if database already exists for this clade.

    Returns:
        Path to existing database, or None
    """
    return find_database_dir(output_dir, clade_name)


def create_database_for_run(
    orthophyl_dir: Path,
    clade_taxonomy: str,
    clade_name: str,
    output_dir: Path,
    force: bool = False,
    subclade_meta: Optional[Dict] = None,
    taxonomy_source: str = "ncbi",
    qc_applied: bool = True,
    register_only: bool = False,
    metadata_file: Optional[str] = None,
) -> Path:
    """
    Create a database directory for a single OrthoPhyl run.

    metadata_file (optional): path to an assembly metadata TSV (e.g.
    all_asm_acc_metadata) to copy into the DB, independent of subclade_meta.

    subclade_meta (optional): metadata for a MASH-partitioned subclade. When
        provided, extra fields (is_subclade, parent_taxon, subclade_id,
        sketch_file, members_file, built, source_genome_dir, is_backbone) are
        written to the config, and any source sketch/members files are copied
        into the DB dir. A plain (non-subclade) DB is created byte-identically
        when this is None.

    taxonomy_source: "ncbi" (default) or "user_supplied" -- provenance of
        clade_taxonomy, surfaced for the local-genome-ingest mode where the
        taxon name/lineage is not NCBI-assigned.

    qc_applied: whether CheckM2 QC was run on the genomes before tree-building.
        Defaults to True; the local-genome-ingest mode's --skip-qc sets it False.

    register_only (optional): register a subclade DB WITHOUT a built tree -- for
        lazy build-on-demand. Skips OrthoPhyl-run validation, sets built=false /
        has_trees=false, writes a placeholder tree, and skips the orthophyl_run
        symlink. Requires subclade_meta.

    Returns:
        Path to created database directory
    """
    logger.info(f"Creating database for: {clade_name}")
    logger.info(f"  Taxonomy: {clade_taxonomy}")
    logger.info(f"  Source: {orthophyl_dir}")

    # Check if already exists
    existing_db = database_exists(clade_name, output_dir)
    if existing_db and not force:
        logger.info(f"  ⊙ Database already exists: {existing_db}")
        logger.info(f"    Use --force to rebuild or --update to skip")
        raise FileExistsError(f"Database already exists: {existing_db}")

    if register_only:
        # Lazy registration: no OrthoPhyl run yet, so nothing to validate. Build a
        #   minimal validation dict so the rest of the function can proceed.
        validation = {
            'valid': True,
            'n_genomes': (subclade_meta or {}).get('n_genomes', 0),
            'has_hmms': False,
            'has_trees': False,
            'tree_file': None,
            'hmm_dir': None,
            'tree_methods': set(),
            'tree_data_types': set(),
            'available_trees': [],
            'validation_details': ['Registered lazily (built=false); tree not yet built'],
            'warnings': [],
        }
    else:
        # Validate OrthoPhyl run
        try:
            validation = validate_orthophyl_run(orthophyl_dir)
        except Exception as e:
            logger.error(f"  ✗ Failed to validate: {e}")
            raise

        if not validation['valid']:
            error_msg = f"Invalid OrthoPhyl run: {orthophyl_dir}\n"
            error_msg += "Validation details:\n"
            for detail in validation['validation_details']:
                error_msg += f"  - {detail}\n"
            raise ValueError(error_msg)

        if validation['warnings']:
            for warning in validation['warnings']:
                logger.warning(f"  ⚠ {warning}")

    # Parse taxonomy to determine clade rank
    tax = GTDBTaxonomy(clade_taxonomy)
    clade_rank = tax.get_most_specific_rank()
    
    if not clade_rank:
        raise ValueError(f"Cannot determine clade rank from taxonomy: {clade_taxonomy}")
    
    logger.info(f"  ✓ Clade rank: {tax.get_rank_name(clade_rank)}")
    
    # Create database directory. Megatree subclade/backbone DBs nest under
    #   their parent taxon's dir; everything else stays flat (db_dir_for).
    sm_for_path = subclade_meta or {}
    nest_under = sm_for_path.get("parent_taxon") if (
        sm_for_path.get("is_subclade") or sm_for_path.get("is_backbone")) else None
    db_dir = db_dir_for(output_dir, clade_name, nest_under)

    if force and db_dir.exists():
        logger.info(f"  ↻ Removing existing database (--force)")
        shutil.rmtree(db_dir)
    
    db_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"  → Database directory: {db_dir}")
    
    # Create database_config.json
    config = {
        "created": datetime.now().isoformat(),
        "last_updated": datetime.now().isoformat(),
        "version": "1.0",
        "database_type": "hierarchical_taxonomy",
        "clade_name": clade_name,
        "clade_taxonomy": clade_taxonomy,
        "clade_rank": clade_rank,
        "clade_rank_name": tax.get_rank_name(clade_rank),
        # Relative to db_dir, same convention as source_genome_dir/sketch_file/
        # members_file below, so it survives a different bind-mount mapping
        # than at creation time (see relative_symlink). Resolved back to
        # absolute on load by assembly_router.py's _load_single_database.
        "orthophyl_source": os.path.relpath(orthophyl_dir.resolve(), db_dir.resolve()),
        "n_genomes": validation['n_genomes'],
        "has_hmms": validation['has_hmms'],
        "has_trees": validation['has_trees'],
        "validation_details": validation['validation_details'],
        # Add available tree methods and data types for ReLeaf compatibility
        "available_tree_methods": sorted(list(validation['tree_methods'])),
        "available_data_types": sorted(list(validation['tree_data_types'])),
        "available_trees": [
            {
                "filename": t['name'],
                "method": t['method'],
                "data_type": t['data_type']
            }
            for t in validation['available_trees']
        ],
        # Taxon source metadata (for auto-gather feature)
        "source_taxon_name": None,
        "source_taxid": None,
        "source_rank": None,
        # Provenance: was clade_taxonomy NCBI-derived, or user-supplied (local-genome-
        #   ingest mode)? And was CheckM2 QC applied to the input genomes?
        "taxonomy_source": taxonomy_source,
        "qc_applied": qc_applied,
        "assembly_accessions": [],
        "n_assemblies_at_creation": validation['n_genomes'],
        "quality_filters": {
            "min_completeness": None,
            "max_contamination": None,
            "min_n50": None
        }
    }

    # Subclade metadata (MASH-partitioned oversized taxa). Defaults keep a normal
    #   DB byte-identical to before; megatree subclades carry these fields so the
    #   parent taxon / subclade identity are recorded alongside the built tree.
    sm = subclade_meta or {}
    config["is_subclade"] = bool(sm.get("is_subclade", False))
    config["parent_taxon"] = sm.get("parent_taxon")
    config["subclade_id"] = sm.get("subclade_id")
    # A DB is "built" (real tree) unless we're registering it lazily.
    config["built"] = (not register_only)
    src_genome_dir = sm.get("source_genome_dir")
    if src_genome_dir:
        # Stored relative to db_dir so it survives a different bind-mount
        # mapping between build time and load time (see relative_symlink).
        config["source_genome_dir"] = os.path.relpath(
            Path(src_genome_dir).resolve(), db_dir.resolve())
    else:
        config["source_genome_dir"] = None
    # Marks the megatree's sparse overview DB, distinct from its dense subclades.
    config["is_backbone"] = bool(sm.get("is_backbone", False))
    # sketch_file / members_file are set to the in-DB copies below (if provided).
    config["sketch_file"] = None
    config["members_file"] = None
    if register_only:
        config["has_trees"] = False

    # Copy the subclade's MASH sketch + member list into the DB dir so the router
    #   can compare a query against them without reaching back to the run dir.
    src_sketch = sm.get("sketch_file")
    if src_sketch and os.path.exists(src_sketch):
        dest_sketch = db_dir / "subclade_sketch.msh"
        shutil.copy(src_sketch, dest_sketch)
        # Relative (bare filename): resolved against db_dir on load, so this
        # survives a different bind-mount mapping than at creation time.
        config["sketch_file"] = dest_sketch.name
        logger.info("  ✓ Copied subclade sketch")
    src_members = sm.get("members_file")
    if src_members and os.path.exists(src_members):
        dest_members = db_dir / "subclade_members.txt"
        shutil.copy(src_members, dest_members)
        config["members_file"] = dest_members.name
        logger.info("  ✓ Copied subclade member list")

    # Optional assembly metadata (e.g. all_asm_acc_metadata), used by
    #   python_scripts/annotate_tree_by_country.py for --geo figures. Not
    #   subclade-specific, so it's read from the top-level metadata_file arg
    #   rather than subclade_meta.
    config["metadata_file"] = None
    if metadata_file and os.path.exists(metadata_file):
        dest_metadata = db_dir / "assembly_metadata.tsv"
        shutil.copy(metadata_file, dest_metadata)
        config["metadata_file"] = dest_metadata.name
        logger.info("  ✓ Copied assembly metadata")

    config_file = db_dir / "database_config.json"
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)

    logger.info(f"  ✓ Created config")

    # Copy or link phylogeny
    if register_only:
        # No tree yet -- lazy registration writes a placeholder and skips the symlink.
        (db_dir / "phylogeny.nwk").write_text(
            f"# Placeholder tree for {clade_name} (registered lazily, built=false)\n")
        logger.info("  ⚠ Registered lazily: placeholder tree, no orthophyl_run symlink")
    else:
        if validation['tree_file']:
            dest_tree = db_dir / "phylogeny.nwk"
            shutil.copy(validation['tree_file'], dest_tree)
            logger.info(f"  ✓ Copied tree: {validation['tree_file'].name}")
        else:
            (db_dir / "phylogeny.nwk").write_text(f"# Placeholder tree for {clade_name}\n")
            logger.info("  ⚠ Created placeholder tree (no tree file found)")

        # Create symlink to OrthoPhyl run
        orthophyl_link = db_dir / "orthophyl_run"
        if orthophyl_link.exists():
            orthophyl_link.unlink()

        try:
            relative_symlink(orthophyl_link, orthophyl_dir)
            logger.info(f"  ✓ Created symlink: orthophyl_run")
        except OSError as e:
            logger.warning(f"  ⚠ Could not create symlink: {e}")
            logger.info("  → Copying essential files instead...")
            orthophyl_link.mkdir(exist_ok=True)

            # Copy HMMs if found
            if validation['hmm_dir']:
                dest_hmm = orthophyl_link / "hmms"
                shutil.copytree(validation['hmm_dir'], dest_hmm, dirs_exist_ok=True)
                logger.info(f"  ✓ Copied HMMs")

    # Create genome list
    genome_list_file = db_dir / "genome_list.txt"
    genome_names = []

    if register_only:
        # No OrthoPhyl run yet -- the genome list IS the subclade's member list.
        # config["members_file"] is now relative to db_dir; resolve it here.
        if config.get("members_file"):
            members_path = db_dir / config["members_file"]
            if members_path.exists():
                with open(members_path) as f:
                    genome_names = [line.strip() for line in f if line.strip()]
    else:
        # Try to extract genome names
        for source in ['genome_list', 'all_input_list']:
            for base_dir in [orthophyl_dir, orthophyl_dir / "store"]:
                source_file = base_dir / source
                if source_file.exists():
                    with open(source_file, 'r') as f:
                        genome_names = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                    if genome_names:
                        break
            if genome_names:
                break

        if not genome_names:
            # Fallback to protein files
            for prot_dir in [orthophyl_dir / "annots_prots", orthophyl_dir / "store" / "annots_prots"]:
                if prot_dir.exists():
                    genome_names = [f.stem for f in prot_dir.glob("*.faa")]
                    if genome_names:
                        break

    with open(genome_list_file, 'w') as f:
        f.write(f"# Genomes in {clade_name} database\n")
        f.write(f"# Created: {datetime.now().isoformat()}\n")
        for genome in sorted(genome_names):
            f.write(f"{genome}\n")
    
    logger.info(f"  ✓ Created genome list: {len(genome_names)} genomes")
    
    # Create README
    readme = db_dir / "README.txt"
    with open(readme, 'w') as f:
        f.write(f"Hierarchical Taxonomy Database: {clade_name}\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Created: {config['created']}\n")
        f.write(f"Clade: {clade_name}\n")
        f.write(f"Taxonomy: {clade_taxonomy}\n")
        f.write(f"Rank: {config['clade_rank_name']}\n")
        f.write(f"Genomes: {validation['n_genomes']}\n\n")
        
        f.write("Validation Details:\n")
        for detail in validation['validation_details']:
            f.write(f"  - {detail}\n")
        f.write("\n")
        
        f.write("Source OrthoPhyl Run:\n")
        f.write(f"  {orthophyl_dir}\n\n")
        
        f.write("Use with assembly_router.py:\n")
        f.write(f"  python assembly_router.py \\\n")
        f.write(f"      --assembly genome.fna \\\n")
        f.write(f"      --taxonomy \"d__...; ...\" \\\n")
        f.write(f"      --database-dir {db_dir}\n\n")
        
        if validation['warnings']:
            f.write("Warnings:\n")
            for warning in validation['warnings']:
                f.write(f"  - {warning}\n")
    
    logger.info(f"  ✓ Created README")
    logger.info(f"✓ Database created successfully\n")
    
    return db_dir


def get_existing_databases(output_dir: Path) -> Set[str]:
    """
    Get set of clade names for existing databases.
    
    Returns:
        Set of clade names
    """
    existing = set()

    for db_dir in iter_database_dirs(output_dir):
        config_file = db_dir / "database_config.json"
        try:
            with open(config_file, 'r') as f:
                config = json.load(f)
            existing.add(config['clade_name'])
        except:
            pass

    return existing


def remove_databases(names: List[str], output_dir: Path, confirm: bool = False) -> Dict:
    """
    Remove one or more clade databases by name from output_dir.

    Resolves each name via database_exists() (the same safe-name transform
    used everywhere else in this file). Default (confirm=False) is a DRY
    RUN: nothing is deleted, callers get back what WOULD be removed. With
    confirm=True, each resolved _db directory is actually shutil.rmtree'd
    (safe against its orthophyl_run symlink -- rmtree unlinks a symlink
    entry rather than following it into the target), and the master index
    is refreshed exactly once at the end from the post-removal directory
    state (never once per name).

    A name that does not resolve to an existing _db directory is reported
    in 'not_found' and does not stop processing of the remaining names,
    matching this file's existing --update-mode tone of warning and
    continuing rather than hard-failing on one bad entry.

    output_dir not existing at all means every name is 'not_found' -- no
    directory is created as a side effect of a remove call.

    Returns a dict: {'removed': [...], 'not_found': [...], 'dry_run': bool},
    where each 'removed' entry is {'clade_name', 'database_dir',
    'clade_taxonomy', 'n_genomes'} (taxonomy/n_genomes are None if the
    config was missing or unreadable -- the directory is still resolved
    and reported, just without that context).
    """
    output_dir = Path(output_dir)
    removed = []
    not_found = []

    if not output_dir.exists():
        return {'removed': [], 'not_found': list(names), 'dry_run': not confirm}

    for name in names:
        db_dir = database_exists(name, output_dir)
        if db_dir is None:
            logger.warning(f"  ⊙ No database found for clade '{name}'; skipping")
            not_found.append(name)
            continue

        clade_taxonomy = None
        n_genomes = None
        config_file = db_dir / "database_config.json"
        if config_file.exists():
            try:
                with open(config_file) as f:
                    config = json.load(f)
                clade_taxonomy = config.get('clade_taxonomy')
                n_genomes = config.get('n_genomes')
            except (json.JSONDecodeError, OSError):
                pass

        if confirm:
            shutil.rmtree(db_dir)
            logger.info(f"  ✓ Removed {db_dir}")
        else:
            logger.info(f"  [DRY RUN] Would remove {db_dir} "
                        f"(taxonomy={clade_taxonomy}, n_genomes={n_genomes})")

        removed.append({
            'clade_name': name,
            'database_dir': db_dir,
            'clade_taxonomy': clade_taxonomy,
            'n_genomes': n_genomes,
        })

    if confirm and removed:
        all_dbs = _scan_all_databases_for_index(output_dir)
        create_master_index(all_dbs, output_dir)

    return {'removed': removed, 'not_found': not_found, 'dry_run': not confirm}


def _scan_all_databases_for_index(output_dir: Path) -> List[Dict]:
    """
    Scan every *_db directory under output_dir and read back the 5-key dict
    shape create_master_index expects (clade_name, clade_taxonomy,
    database_dir, n_genomes, taxonomy_source). Skips a _db dir whose config
    is missing or unreadable rather than raising.

    This is the read-only rescan used to refresh the master index after any
    change to the set of databases on disk (add, update, or remove) -- kept
    separate from get_existing_databases (which returns only a Set[str] of
    names, for a different purpose).
    """
    all_dbs = []
    for db_dir in iter_database_dirs(output_dir):
        config_file = db_dir / "database_config.json"
        try:
            with open(config_file) as f:
                config = json.load(f)
            all_dbs.append({
                'clade_name': config['clade_name'],
                'clade_taxonomy': config['clade_taxonomy'],
                'database_dir': db_dir,
                'n_genomes': config['n_genomes'],
                'taxonomy_source': config.get('taxonomy_source', 'ncbi'),
            })
        except (json.JSONDecodeError, KeyError, OSError):
            continue
    return all_dbs


def create_master_index(databases: List[Dict], output_dir: Path):
    """
    Create master index of all databases.
    """
    index_file = output_dir / "database_index.json"
    
    index = {
        "created": datetime.now().isoformat(),
        "version": "1.0",
        "n_databases": len(databases),
        "databases": []
    }
    
    for db in databases:
        tax = GTDBTaxonomy(db['clade_taxonomy'])

        entry = {
            "clade_name": db['clade_name'],
            "clade_taxonomy": db['clade_taxonomy'],
            "clade_rank": tax.get_most_specific_rank(),
            "clade_rank_name": tax.get_rank_name(tax.get_most_specific_rank()),
            # Relative to output_dir so the index survives a different
            # bind-mount mapping than at write time (see assembly_router.py's
            # _load_databases, which resolves it against self.database_dir).
            "database_dir": os.path.relpath(
                Path(db['database_dir']).resolve(), output_dir.resolve()),
            "n_genomes": db['n_genomes']
        }
        # Surface subclade/built status in the index when known (read back from
        #   each DB's config by the callers that assemble `databases`).
        if 'is_subclade' in db:
            entry['is_subclade'] = db['is_subclade']
        if 'built' in db:
            entry['built'] = db['built']
        if 'taxonomy_source' in db:
            entry['taxonomy_source'] = db['taxonomy_source']
        index["databases"].append(entry)
    
    with open(index_file, 'w') as f:
        json.dump(index, f, indent=2)
    
    logger.info(f"✓ Created master index: {index_file}")
    
    # Human-readable summary
    summary_file = output_dir / "database_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("HIERARCHICAL TAXONOMY DATABASES - SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"Created: {index['created']}\n")
        f.write(f"Total databases: {len(databases)}\n\n")
        
        f.write("Available Databases:\n")
        f.write("-" * 70 + "\n")
        
        for db in sorted(databases, key=lambda x: x['clade_name']):
            tax = GTDBTaxonomy(db['clade_taxonomy'])
            rank_name = tax.get_rank_name(tax.get_most_specific_rank())
            
            f.write(f"\n{db['clade_name']} ({rank_name})\n")
            f.write(f"  Taxonomy: {db['clade_taxonomy']}\n")
            f.write(f"  Genomes: {db['n_genomes']}\n")
            f.write(f"  Database: {db['database_dir']}\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("USAGE\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("Route an assembly:\n\n")
        f.write("  python assembly_router.py \\\n")
        f.write("      --assembly genome.fna \\\n")
        f.write("      --taxonomy \"d__Bacteria;p__...; ...\" \\\n")
        f.write("      --database-dir [one of the databases above]\n")
    
    logger.info(f"✓ Created summary: {summary_file}")


def parse_input_table(input_file: Path) -> List[Dict]:
    """
    Parse input TSV: clade_name<TAB>orthophyl_dir<TAB>clade_taxonomy

    A relative orthophyl_dir is resolved against the TSV's own directory
    (the wrapper's _create_database_entry writes orthophyl_runs.tsv this
    way, relative to --database-dir, so the TSV survives a different
    bind-mount mapping than at write time). An absolute path is used as-is,
    for backward compatibility with pre-existing TSVs.
    """
    runs = []
    input_file = Path(input_file)
    tsv_dir = input_file.resolve().parent

    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue

            if line_num == 1 and 'clade' in line.lower():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 3:
                logger.warning(f"Line {line_num}: expected 3 fields, got {len(fields)}. Skipping.")
                continue

            clade_name = fields[0].strip()
            orthophyl_dir = Path(fields[1].strip())
            if not orthophyl_dir.is_absolute():
                orthophyl_dir = tsv_dir / orthophyl_dir
            clade_taxonomy = fields[2].strip()
            
            if not clade_name or not clade_taxonomy:
                logger.warning(f"Line {line_num}: empty field(s), skipping")
                continue
            
            runs.append({
                'clade_name': clade_name,
                'orthophyl_dir': orthophyl_dir,
                'clade_taxonomy': clade_taxonomy
            })
    
    logger.info(f"Parsed {len(runs)} entries from {input_file}")
    return runs


def _infer_strip_prefix(output_dir: Path) -> Optional[str]:
    """Infer an old container mount prefix (e.g. "/databases") from
    database_index.json, for --migrate.

    database_index.json has stored database_dir relative to output_dir since
    commit 852ddd8; a DB directory created before that commit still has an
    absolute, container-local path baked into its index entry (and several
    config path fields). If an entry's database_dir is absolute and does not
    exist as given, but output_dir/<basename of that path> does exist, the
    entry's parent directory is that old prefix.
    """
    index_file = output_dir / "database_index.json"
    if not index_file.exists():
        return None
    try:
        with open(index_file) as f:
            index = json.load(f)
    except (json.JSONDecodeError, OSError):
        return None
    for entry in index.get("databases", []):
        raw = entry.get("database_dir")
        if not raw:
            continue
        p = Path(raw)
        if p.is_absolute() and not p.exists():
            candidate = output_dir / p.name
            if candidate.exists():
                return str(p.parent)
    return None


def _migrate_remap(value: Optional[str], strip_prefix: Optional[str],
                    output_dir: Path) -> Optional[str]:
    """Rewrite a value living under strip_prefix to live under output_dir
    instead; pass through unchanged otherwise (including None/no prefix)."""
    if not value or not strip_prefix:
        return value
    if value == strip_prefix:
        return str(output_dir)
    prefix_with_sep = strip_prefix.rstrip(os.sep) + os.sep
    if value.startswith(prefix_with_sep):
        return str(output_dir / value[len(prefix_with_sep):])
    return value


def _migrate_resolve_field(value: Optional[str], db_dir: Path,
                            strip_prefix: Optional[str],
                            output_dir: Path) -> Optional[Path]:
    """Resolve a config path field to a Path anchored at db_dir's CURRENT
    (pre-move) location -- handles both the current relative-to-db_dir
    convention and a legacy absolute/mount-prefixed value.

    A value that was absolute BEFORE remapping is anchored at output_dir by
    _migrate_remap (which may itself be relative to cwd, e.g. in tests) --
    it must not then be re-joined with db_dir, or the result double-nests
    under db_dir and silently fails every exists() check below.
    """
    if not value:
        return None
    was_absolute = Path(value).is_absolute()
    remapped = _migrate_remap(value, strip_prefix, output_dir)
    p = Path(remapped)
    if not was_absolute:
        p = db_dir / p
    return p


def _migrate_relativize_runs_tsv(output_dir: Path, strip_prefix: Optional[str]) -> None:
    """Rewrite orthophyl_runs.tsv's orthophyl_dir column (field 2) relative
    to output_dir, remapping any legacy mount-prefixed absolute path first."""
    tsv_path = output_dir / "orthophyl_runs.tsv"
    if not tsv_path.exists():
        return
    lines = tsv_path.read_text().splitlines()
    new_lines = []
    for line in lines:
        if line.startswith('#') or not line.strip():
            new_lines.append(line)
            continue
        fields = line.split('\t')
        if len(fields) >= 2:
            remapped = _migrate_remap(fields[1], strip_prefix, output_dir)
            p = Path(remapped)
            if p.is_absolute():
                try:
                    p = Path(os.path.relpath(p, output_dir.resolve()))
                except ValueError:
                    pass
            fields[1] = str(p)
        new_lines.append('\t'.join(fields))
    tsv_path.write_text('\n'.join(new_lines) + ('\n' if new_lines else ''))


def migrate_databases(output_dir: Path, strip_prefix: Optional[str] = None,
                       dry_run: bool = True) -> Dict:
    """One-off repair for a pre-nesting / pre-portability database directory.

    Three problems fixed together (the second needs paths the first would
    otherwise move out from under, so both are planned from each DB's
    config before anything is touched on disk):
      1. A legacy absolute container mount prefix (e.g. "/databases/...",
         baked into configs/index written before commit 852ddd8) is
         remapped back to a real path under output_dir -- inferred from
         database_index.json unless strip_prefix is given explicitly.
      2. Megatree subclade/backbone DBs (is_subclade/is_backbone +
         parent_taxon) are moved from the old flat layout into the nested
         <output_dir>/<parent_taxon>/<name>_db layout.
      3. Every per-DB path field (sketch_file, members_file, metadata_file,
         source_genome_dir, orthophyl_source) is rewritten relative to the
         DB's (possibly new) db_dir, and orthophyl_run is re-created with
         relative_symlink -- a move changes the correct relative offset
         even when the symlink's target hasn't itself moved.

    Dangling orthophyl_run targets (the source .pipeline_runs dir was
    deleted) are reported, never deleted or skipped -- the DB's own config/
    sketch/tree still stand alone and are still migrated; remove the DB
    afterwards with --remove --yes if it's unusable without its run dir.

    dry_run=True (the default) computes and reports every move/rewrite
    without changing anything on disk. Idempotent: a second run against
    already-migrated DBs reports no moves and no changes.

    Returns {'strip_prefix', 'dry_run', 'moved': [...], 'rewritten': [...],
    'dangling_links': [...]}.
    """
    output_dir = Path(output_dir)
    if strip_prefix is None:
        strip_prefix = _infer_strip_prefix(output_dir)
    logger.info(f"  Mount prefix to strip: {strip_prefix or '(none inferred)'}")

    report = {'strip_prefix': strip_prefix, 'dry_run': dry_run,
              'moved': [], 'rewritten': [], 'dangling_links': []}

    PATH_FIELDS = ('sketch_file', 'members_file', 'metadata_file',
                   'source_genome_dir', 'orthophyl_source')

    # Plan everything from each DB's CURRENT (pre-move) location first.
    plans = []
    for db_dir in list(iter_database_dirs(output_dir)):
        config_file = db_dir / "database_config.json"
        try:
            with open(config_file) as f:
                config = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            logger.warning(f"  ⊙ Skipping {db_dir}: unreadable config ({e})")
            continue

        resolved = {
            field: _migrate_resolve_field(config.get(field), db_dir,
                                           strip_prefix, output_dir)
            for field in PATH_FIELDS
        }

        link = db_dir / "orthophyl_run"
        link_target = None
        if link.is_symlink():
            raw_target = os.readlink(link)
            if os.path.isabs(raw_target):
                link_target = Path(_migrate_remap(raw_target, strip_prefix, output_dir))
            else:
                link_target = (link.parent / raw_target).resolve()

        nest_under = config.get("parent_taxon") if (
            config.get("is_subclade") or config.get("is_backbone")) else None
        new_db_dir = db_dir_for(output_dir, config.get('clade_name', db_dir.name),
                                 nest_under)

        plans.append({
            'db_dir': db_dir, 'new_db_dir': new_db_dir, 'config': config,
            'resolved': resolved, 'link_target': link_target,
        })

    for plan in plans:
        db_dir, new_db_dir = plan['db_dir'], plan['new_db_dir']
        clade_name = plan['config'].get('clade_name', db_dir.name)

        if db_dir != new_db_dir:
            report['moved'].append({'clade_name': clade_name,
                                     'from': str(db_dir), 'to': str(new_db_dir)})
            prefix = "[DRY RUN] Would move" if dry_run else "Moving"
            logger.info(f"  {prefix} {db_dir} -> {new_db_dir}")
            if not dry_run:
                new_db_dir.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(db_dir), str(new_db_dir))

        changes = {}
        for field, abs_path in plan['resolved'].items():
            if abs_path is None:
                continue
            # A field pointing inside db_dir (e.g. sketch_file, copied in at
            # DB creation) moves along with it, so its CURRENT location by
            # the time this runs (in a live pass, db_dir may already have
            # been moved to new_db_dir above) is under new_db_dir, not
            # abs_path. Compare un-resolved (abs_path's db_dir-relative
            # fields were built as a literal db_dir / p join, never
            # resolve()'d) so this stays a pure string computation and
            # doesn't need db_dir to still exist on disk.
            try:
                rel_to_old = abs_path.relative_to(db_dir)
                final_abs = new_db_dir / rel_to_old
            except ValueError:
                final_abs = abs_path
            if final_abs.exists():
                new_val = os.path.relpath(final_abs.resolve(), new_db_dir.resolve())
            else:
                new_val = plan['config'].get(field)
            if new_val != plan['config'].get(field):
                changes[field] = new_val

        link_target = plan['link_target']
        if link_target is not None and not link_target.exists():
            report['dangling_links'].append({
                'clade_name': clade_name, 'db_dir': str(new_db_dir),
                'target': str(link_target),
            })

        if changes or (link_target is not None and link_target.exists()):
            report['rewritten'].append({'clade_name': clade_name,
                                         'db_dir': str(new_db_dir),
                                         'changes': changes})

        if dry_run:
            continue

        if changes:
            plan['config'].update(changes)
            with open(new_db_dir / "database_config.json", 'w') as f:
                json.dump(plan['config'], f, indent=2)

        if link_target is not None and link_target.exists():
            new_link = new_db_dir / "orthophyl_run"
            # Re-create unconditionally: a move changes the correct relative
            # offset even when the target itself hasn't changed.
            if new_link.exists() or new_link.is_symlink():
                new_link.unlink()
            relative_symlink(new_link, link_target)
        elif link_target is not None:
            logger.warning(f"  ⚠ Dangling orthophyl_run target for "
                            f"{clade_name}: {link_target} (left as-is)")

    if not dry_run:
        _migrate_relativize_runs_tsv(output_dir, strip_prefix)
        all_dbs = _scan_all_databases_for_index(output_dir)
        if all_dbs:
            create_master_index(all_dbs, output_dir)

    return report


def main():
    parser = argparse.ArgumentParser(
        description="Create/update hierarchical taxonomy databases from OrthoPhyl runs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Initial creation
  python OP_database_tool.py \\
      --input orthophyl_runs.tsv \\
      --output-dir taxonomy_databases/
  
  # Add new entries (skip existing)
  python OP_database_tool.py \\
      --input orthophyl_runs.tsv \\
      --output-dir taxonomy_databases/ \\
      --update
  
  # Force rebuild all
  python OP_database_tool.py \\
      --input orthophyl_runs.tsv \\
      --output-dir taxonomy_databases/ \\
      --force

Input Format (TSV):
  clade_name<TAB>orthophyl_dir<TAB>clade_taxonomy
  
  Gaiellales	/data/gaiellales_run	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales
  Escherichia	/data/ecoli_run	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;g__Escherichia
        """
    )
    
    parser.add_argument(
        '--input',
        help='TSV: clade_name, orthophyl_dir, clade_taxonomy '
             '(not required with --single-clade)'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for databases'
    )
    parser.add_argument(
        '--update',
        action='store_true',
        help='Update mode: skip existing databases, add new ones only'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force rebuild all databases (overwrite existing)'
    )

    # Register a single clade without a TSV (used for megatree subclade DBs).
    parser.add_argument(
        '--single-clade',
        nargs=3,
        metavar=('NAME', 'TAXONOMY', 'ORTHOPHYL_DIR'),
        help='Register ONE clade directly (name, taxonomy, orthophyl_dir) instead '
             'of parsing --input. With --register-only, ORTHOPHYL_DIR may be a '
             'placeholder (unused).'
    )
    # Subclade metadata (MASH-partitioned oversized taxa).
    parser.add_argument('--is-subclade', action='store_true',
                        help='Mark this DB as a subclade of a larger taxon.')
    parser.add_argument('--is-backbone', action='store_true',
                        help="Mark this DB as a megatree's sparse backbone "
                             '(overview) database, as opposed to a dense subclade.')
    parser.add_argument('--parent-taxon', help='Parent taxon name for a subclade.')
    parser.add_argument('--subclade-id', type=int, help='1-based subclade index.')
    parser.add_argument('--sketch-file', help='MASH .msh sketch for the subclade.')
    parser.add_argument('--members-file', help='Text file of subclade member basenames.')
    parser.add_argument('--source-genome-dir',
                        help='Directory holding the subclade raw member FASTAs '
                             '(for lazy build re-staging).')
    parser.add_argument('--metadata-file',
                        help='Optional assembly metadata TSV (e.g. '
                             'all_asm_acc_metadata from gather_filter_asms.sh) to '
                             'copy into the DB for later use (e.g. --geo figures). '
                             'Not required; the router does not read it.')
    parser.add_argument('--n-genomes', type=int, default=0,
                        help='Member count to record for the subclade (also used '
                             'for a lazily-registered subclade).')
    parser.add_argument('--taxonomy-source', choices=['ncbi', 'user_supplied'],
                        default='ncbi',
                        help='Provenance of clade_taxonomy (default: ncbi). '
                             'user_supplied marks a taxonomy that was not '
                             'NCBI-assigned (e.g. local-genome-ingest mode).')
    parser.add_argument('--qc-not-applied', action='store_true',
                        help='Record that CheckM2 QC was NOT run on the input '
                             'genomes (default: QC is assumed applied).')
    parser.add_argument('--register-only', action='store_true',
                        help='Register a subclade DB with built=false (no tree yet) '
                             'for lazy build-on-demand. Requires --is-subclade.')

    # Remove one or more existing clade databases.
    parser.add_argument('--remove', nargs='+', metavar='NAME',
                        help='Remove one or more clade databases by name. Default is '
                             'a DRY RUN (prints what would be removed, deletes '
                             'nothing); pass --yes to actually delete. Cannot be '
                             'combined with --input or --single-clade.')
    parser.add_argument('--yes', action='store_true',
                        help='Confirm --remove actually deletes the resolved '
                             'database directories and refreshes the master index. '
                             'Has no effect without --remove.')

    # One-off repair for a pre-nesting / pre-portability database directory.
    parser.add_argument('--migrate', action='store_true',
                        help='Repair an existing database directory in place: '
                             'remap a legacy absolute container mount prefix '
                             'back to a real path, move megatree subclade/'
                             'backbone DBs into the nested <parent>/<name>_db '
                             'layout, and rewrite path fields relative to '
                             'each DB. Cannot be combined with --input, '
                             '--single-clade, or --remove.')
    parser.add_argument('--strip-prefix', metavar='PREFIX',
                        help='Old mount prefix to strip during --migrate '
                             '(e.g. "/databases"). If omitted, inferred from '
                             'database_index.json.')
    parser.add_argument('--dry-run', action='store_true',
                        help='With --migrate, print every move/rewrite '
                             'without changing anything on disk.')

    args = parser.parse_args()

    if args.migrate and (args.input or args.single_clade or args.remove):
        parser.error("--migrate cannot be combined with --input, --single-clade, or --remove")
    if args.strip_prefix and not args.migrate:
        parser.error("--strip-prefix only applies to --migrate")
    if args.dry_run and not args.migrate:
        parser.error("--dry-run only applies to --migrate")
    if args.migrate:
        report = migrate_databases(Path(args.output_dir),
                                     strip_prefix=args.strip_prefix,
                                     dry_run=args.dry_run)
        logger.info("")
        logger.info("=" * 70)
        logger.info(f"MIGRATION {'PREVIEW (dry run)' if report['dry_run'] else 'COMPLETE'}")
        logger.info("=" * 70)
        logger.info(f"Moved: {len(report['moved'])}")
        for m in report['moved']:
            logger.info(f"  - {m['clade_name']}: {m['from']} -> {m['to']}")
        logger.info(f"Rewritten: {len(report['rewritten'])}")
        for r in report['rewritten']:
            if r['changes']:
                logger.info(f"  - {r['clade_name']} ({r['db_dir']}): {r['changes']}")
        if report['dangling_links']:
            logger.warning(f"Dangling orthophyl_run links (left in place): "
                            f"{len(report['dangling_links'])}")
            for d in report['dangling_links']:
                logger.warning(f"  - {d['clade_name']} ({d['db_dir']}) -> {d['target']}")
        return 0

    if args.remove and (args.input or args.single_clade):
        parser.error("--remove cannot be combined with --input or --single-clade")
    if args.yes and not args.remove:
        parser.error("--yes only applies to --remove")
    if args.remove:
        result = remove_databases(args.remove, Path(args.output_dir), confirm=args.yes)
        for r in result['removed']:
            action = "Removed" if args.yes else "[DRY RUN] Would remove"
            logger.info(f"{action}: {r['clade_name']} ({r['database_dir']})")
        for name in result['not_found']:
            logger.warning(f"Not found, skipped: {name}")
        return 0

    if args.update and args.force:
        parser.error("Cannot use --update and --force together")
    if not args.input and not args.single_clade:
        parser.error("Provide either --input <tsv> or --single-clade NAME TAX DIR")
    if args.register_only and not args.is_subclade:
        parser.error("--register-only requires --is-subclade")

    # Assemble subclade_meta once (passed to create_database_for_run below).
    subclade_meta = None
    if args.is_subclade or args.is_backbone:
        subclade_meta = {
            'is_subclade': args.is_subclade,
            'is_backbone': args.is_backbone,
            'parent_taxon': args.parent_taxon,
            'subclade_id': args.subclade_id,
            'sketch_file': args.sketch_file,
            'members_file': args.members_file,
            'source_genome_dir': args.source_genome_dir,
            'n_genomes': args.n_genomes,
        }

    # --single-clade short-circuits the TSV workflow: build/register exactly one DB.
    if args.single_clade:
        name, taxonomy, op_dir = args.single_clade
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        try:
            db_dir = create_database_for_run(
                Path(op_dir), taxonomy, name, output_dir,
                force=args.force, subclade_meta=subclade_meta,
                taxonomy_source=args.taxonomy_source,
                qc_applied=not args.qc_not_applied,
                register_only=args.register_only,
                metadata_file=args.metadata_file,
            )
        except FileExistsError:
            logger.info("Database already exists; nothing to do (use --force to rebuild)")
            return 0
        # Refresh the master index across all DBs present.
        # NOTE: is_subclade/built are intentionally NOT surfaced by the scan --
        #   the router re-reads each DB's own database_config.json directly
        #   (assembly_router.py:_load_single_database) rather than trusting
        #   this index for anything but the directory path, so it doesn't
        #   need to populate them.
        all_dbs = _scan_all_databases_for_index(output_dir)
        if all_dbs:
            create_master_index(all_dbs, output_dir)
        logger.info(f"✓ {'Registered' if args.register_only else 'Created'} {db_dir}")
        return 0
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("Hierarchical Taxonomy Database Builder")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {output_dir}")
    if args.update:
        logger.info("Mode: UPDATE (skip existing)")
    elif args.force:
        logger.info("Mode: FORCE (rebuild all)")
    else:
        logger.info("Mode: CREATE (error on existing)")
    logger.info("")
    
    # Get existing databases if in update mode
    existing_clades = set()
    if args.update:
        existing_clades = get_existing_databases(output_dir)
        if existing_clades:
            logger.info(f"Found {len(existing_clades)} existing databases:")
            for clade in sorted(existing_clades):
                logger.info(f"  - {clade}")
            logger.info("")
    
    # Parse input
    try:
        runs = parse_input_table(Path(args.input))
    except Exception as e:
        logger.error(f"Failed to parse input: {e}")
        return 1
    
    if not runs:
        logger.error("No valid entries in input")
        return 1
    
    # Filter if in update mode
    if args.update:
        original_count = len(runs)
        runs = [r for r in runs if r['clade_name'] not in existing_clades]
        skipped = original_count - len(runs)
        if skipped > 0:
            logger.info(f"Skipping {skipped} existing databases")
        if not runs:
            logger.info("No new databases to create")
            # Still update index
            all_dbs = _scan_all_databases_for_index(output_dir)
            if all_dbs:
                create_master_index(all_dbs, output_dir)
            return 0
        logger.info(f"Processing {len(runs)} new databases\n")
    
    # Create databases
    databases = []
    failed = []
    skipped = []
    
    for i, run in enumerate(runs, 1):
        logger.info(f"[{i}/{len(runs)}] {run['clade_name']}")
        
        try:
            db_dir = create_database_for_run(
                run['orthophyl_dir'],
                run['clade_taxonomy'],
                run['clade_name'],
                output_dir,
                force=args.force,
                subclade_meta=subclade_meta,
                taxonomy_source=args.taxonomy_source,
                qc_applied=not args.qc_not_applied,
                metadata_file=args.metadata_file,
            )
            
            # Read back config
            with open(db_dir / "database_config.json", 'r') as f:
                config = json.load(f)
            
            databases.append({
                'clade_name': run['clade_name'],
                'clade_taxonomy': run['clade_taxonomy'],
                'database_dir': db_dir,
                'n_genomes': config['n_genomes'],
                'taxonomy_source': config.get('taxonomy_source', 'ncbi'),
            })
            
        except FileExistsError:
            skipped.append(run['clade_name'])
            logger.info("")
        except Exception as e:
            logger.error(f"✗ Failed: {e}")
            failed.append(run['clade_name'])
            logger.info("")
    
    # In update mode, also load existing databases for index (skip any
    # clade already present from the runs we just processed above).
    if args.update:
        existing_names = {d['clade_name'] for d in databases}
        for scanned in _scan_all_databases_for_index(output_dir):
            if scanned['clade_name'] not in existing_names:
                databases.append(scanned)

    # Create master index
    if databases:
        logger.info("Creating master index...")
        create_master_index(databases, output_dir)
    
    # Final summary
    logger.info("")
    logger.info("=" * 70)
    logger.info("DATABASE BUILD COMPLETE")
    logger.info("=" * 70)
    new_databases = len(databases) - len(existing_clades) if args.update else len(databases)
    if skipped:
        logger.info(f"Skipped (already exist): {len(skipped)}")
    if failed:
        logger.info(f"Failed: {len(failed)}")
        for name in failed:
            logger.info(f"  - {name}")
    
    logger.info(f"\nTotal databases: {len(databases)}")
    logger.info(f"Location: {output_dir}")
    logger.info(f"Index: {output_dir / 'database_index.json'}")
    
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
