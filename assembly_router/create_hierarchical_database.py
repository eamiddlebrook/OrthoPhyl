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
    python create_hierarchical_database.py \
        --input orthophyl_runs.tsv \
        --output-dir taxonomy_databases/
    
    # Update with new entries (skips existing)
    python create_hierarchical_database.py \
        --input orthophyl_runs.tsv \
        --output-dir taxonomy_databases/ \
        --update
    
    # Force rebuild all
    python create_hierarchical_database.py \
        --input orthophyl_runs.tsv \
        --output-dir taxonomy_databases/ \
        --force

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
    safe_name = clade_name.replace(' ', '_').replace('/', '_')
    db_dir = output_dir / f"{safe_name}_db"
    
    if db_dir.exists() and (db_dir / "database_config.json").exists():
        return db_dir
    return None


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
) -> Path:
    """
    Create a database directory for a single OrthoPhyl run.

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
    
    # Create database directory
    safe_name = clade_name.replace(' ', '_').replace('/', '_')
    db_dir = output_dir / f"{safe_name}_db"
    
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
        "orthophyl_source": str(orthophyl_dir.resolve()),
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
    config["source_genome_dir"] = sm.get("source_genome_dir")
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
        config["sketch_file"] = str(dest_sketch)
        logger.info("  ✓ Copied subclade sketch")
    src_members = sm.get("members_file")
    if src_members and os.path.exists(src_members):
        dest_members = db_dir / "subclade_members.txt"
        shutil.copy(src_members, dest_members)
        config["members_file"] = str(dest_members)
        logger.info("  ✓ Copied subclade member list")

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
            orthophyl_link.symlink_to(orthophyl_dir.resolve())
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
        if config.get("members_file") and os.path.exists(config["members_file"]):
            with open(config["members_file"]) as f:
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
    
    for db_dir in output_dir.glob("*_db"):
        config_file = db_dir / "database_config.json"
        if config_file.exists():
            try:
                with open(config_file, 'r') as f:
                    config = json.load(f)
                existing.add(config['clade_name'])
            except:
                pass
    
    return existing


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
            "database_dir": str(db['database_dir']),
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
    """
    runs = []
    
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


def main():
    parser = argparse.ArgumentParser(
        description="Create/update hierarchical taxonomy databases from OrthoPhyl runs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Initial creation
  python create_hierarchical_database.py \\
      --input orthophyl_runs.tsv \\
      --output-dir taxonomy_databases/
  
  # Add new entries (skip existing)
  python create_hierarchical_database.py \\
      --input orthophyl_runs.tsv \\
      --output-dir taxonomy_databases/ \\
      --update
  
  # Force rebuild all
  python create_hierarchical_database.py \\
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

    args = parser.parse_args()

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
            )
        except FileExistsError:
            logger.info("Database already exists; nothing to do (use --force to rebuild)")
            return 0
        # Refresh the master index across all DBs present.
        all_dbs = []
        for d in output_dir.glob("*_db"):
            cfg = d / "database_config.json"
            if cfg.exists():
                with open(cfg) as f:
                    c = json.load(f)
                all_dbs.append({
                    'clade_name': c['clade_name'],
                    'clade_taxonomy': c['clade_taxonomy'],
                    'database_dir': d,
                    'n_genomes': c['n_genomes'],
                    'taxonomy_source': c.get('taxonomy_source', 'ncbi'),
                    # NOTE: is_subclade/built are intentionally NOT surfaced here --
                    #   the router re-reads each DB's own database_config.json
                    #   directly (assembly_router.py:_load_single_database) rather
                    #   than trusting this index for anything but the directory
                    #   path, so this loop does not need to populate them.
                })
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
            all_dbs = []
            for db_dir in output_dir.glob("*_db"):
                config_file = db_dir / "database_config.json"
                if config_file.exists():
                    with open(config_file) as f:
                        config = json.load(f)
                    all_dbs.append({
                        'clade_name': config['clade_name'],
                        'clade_taxonomy': config['clade_taxonomy'],
                        'database_dir': db_dir,
                        'n_genomes': config['n_genomes'],
                        'taxonomy_source': config.get('taxonomy_source', 'ncbi'),
                    })
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
    
    # In update mode, also load existing databases for index
    if args.update:
        for db_dir in output_dir.glob("*_db"):
            config_file = db_dir / "database_config.json"
            if config_file.exists():
                with open(config_file) as f:
                    config = json.load(f)
                # Only add if not already in our new list
                if not any(d['clade_name'] == config['clade_name'] for d in databases):
                    databases.append({
                        'clade_name': config['clade_name'],
                        'clade_taxonomy': config['clade_taxonomy'],
                        'database_dir': db_dir,
                        'n_genomes': config['n_genomes'],
                        'taxonomy_source': config.get('taxonomy_source', 'ncbi'),
                    })
    
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
