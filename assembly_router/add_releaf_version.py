#!/usr/bin/env python3
"""
Add ReLeaf Output as New Database Version

This script takes a ReLeaf output directory and creates a new version
of an existing database that can be used for subsequent ReLeaf runs.

Usage:
    python add_releaf_version.py \
        --database-dir databases/Rhizobiaceae_db/ \
        --releaf-output /path/to/store/ReLeaf_dir/ \
        --version-name v2_releaf_2024-01-15

Author: HGTool/OrthoPhyl
Date: 2025
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional
import logging
from datetime import datetime
import shutil

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def get_current_version(db_dir: Path) -> str:
    """Get current version name from database."""
    current_link = db_dir / "current"
    if current_link.exists() and current_link.is_symlink():
        return current_link.readlink().name
    
    # Fallback: look for v1
    if (db_dir / "v1_orthophyl_initial").exists():
        return "v1_orthophyl_initial"
    
    raise ValueError(f"Cannot determine current version in {db_dir}")


def get_next_version_number(db_dir: Path) -> int:
    """Get next version number."""
    versions = [d for d in db_dir.glob("v*") if d.is_dir()]
    if not versions:
        return 1
    
    version_numbers = []
    for v in versions:
        try:
            num = int(v.name.split('_')[0][1:])  # Extract number from v1_, v2_, etc.
            version_numbers.append(num)
        except:
            pass
    
    return max(version_numbers) + 1 if version_numbers else 1


def validate_releaf_output(releaf_dir: Path) -> Dict:
    """Validate ReLeaf output directory."""
    releaf_dir = Path(releaf_dir)
    
    if not releaf_dir.exists():
        raise FileNotFoundError(f"ReLeaf directory not found: {releaf_dir}")
    
    validation = {
        'valid': False,
        'has_alignments': False,
        'has_trees': False,
        'warnings': []
    }
    
    # Check for new alignments
    prot_aln = releaf_dir / "new_prot_alignments.trm.nm"
    cds_aln = releaf_dir / "new_CDS_alignments.trm.nm"
    
    if prot_aln.exists() or cds_aln.exists():
        validation['has_alignments'] = True
        logger.info(f"  ✓ Found updated alignments")
    else:
        validation['warnings'].append("No updated alignments found")
    
    # Check for new trees
    tree_dir = releaf_dir / "new_trees"
    if tree_dir.exists():
        tree_files = list(tree_dir.glob("*.tree*"))
        if tree_files:
            validation['has_trees'] = True
            logger.info(f"  ✓ Found {len(tree_files)} new tree(s)")
        else:
            validation['warnings'].append("No new trees found")
    else:
        validation['warnings'].append("new_trees directory not found")
    
    validation['valid'] = validation['has_alignments'] and validation['has_trees']
    
    return validation


def count_genomes_in_releaf(releaf_dir: Path) -> tuple:
    """Count genomes in ReLeaf output. Returns (new_genomes, total_genomes)."""
    genome_list = releaf_dir / "genome_list"
    
    if not genome_list.exists():
        # Try to count from annotations
        annots_dir = releaf_dir / "annots_prots"
        if annots_dir.exists():
            new_genomes = len(list(annots_dir.glob("*.faa")))
            logger.info(f"  Counted {new_genomes} new genomes from annotations")
            return new_genomes, None
        return 0, None
    
    # Count from genome_list
    with open(genome_list, 'r') as f:
        all_genomes = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    # TODO: Need to compare with parent version to get new_genomes count
    return 0, len(all_genomes)


def create_composite_orthophyl_dir(
    base_version_dir: Path,
    releaf_output_dir: Path,
    composite_dir: Path
):
    """
    Create a composite OrthoPhyl structure for ReLeaf.
    
    Combines unchanged files from base version with updated files from ReLeaf.
    """
    logger.info(f"  Creating composite structure...")
    
    composite_dir.mkdir(parents=True, exist_ok=True)
    
    base_orthophyl = base_version_dir / "orthophyl_run"
    if not base_orthophyl.exists():
        raise FileNotFoundError(f"Base orthophyl_run not found: {base_orthophyl}")
    
    # 1. HMMs - from base version (unchanged)
    hmm_source = base_orthophyl / "OG_alignmentsToHMM" / "hmms_final"
    if hmm_source.exists():
        hmm_dest = composite_dir / "OG_alignmentsToHMM" / "hmms_final"
        hmm_dest.parent.mkdir(parents=True, exist_ok=True)
        if hmm_dest.exists():
            hmm_dest.unlink()
        hmm_dest.symlink_to(hmm_source.resolve())
        logger.info(f"    ✓ Linked HMMs from base version")
    
    # 2. Alignments - from ReLeaf output (UPDATED)
    phylo_dir = composite_dir / "phylo_current"
    phylo_dir.mkdir(parents=True, exist_ok=True)
    
    prot_aln = releaf_output_dir / "new_prot_alignments.trm.nm"
    if prot_aln.exists():
        dest = phylo_dir / "AlignmentsProts.trm"
        if dest.exists():
            dest.unlink()
        dest.symlink_to(prot_aln.resolve())
        logger.info(f"    ✓ Linked protein alignments from ReLeaf")
    
    cds_aln = releaf_output_dir / "new_CDS_alignments.trm.nm"
    if cds_aln.exists():
        dest = phylo_dir / "AlignmentsCDS.trm"
        if dest.exists():
            dest.unlink()
        dest.symlink_to(cds_aln.resolve())
        logger.info(f"    ✓ Linked CDS alignments from ReLeaf")
    
    # 3. Trimmed columns - from base version (unchanged)
    trim_source = base_orthophyl / "phylo_current" / "trimmed_columns"
    if trim_source.exists():
        trim_dest = phylo_dir / "trimmed_columns"
        if trim_dest.exists():
            trim_dest.unlink()
        trim_dest.symlink_to(trim_source.resolve())
        logger.info(f"    ✓ Linked trimmed columns from base version")
    
    # 4. Trees - from ReLeaf output (UPDATED)
    trees_dir = composite_dir / "FINAL_SPECIES_TREES"
    trees_dir.mkdir(parents=True, exist_ok=True)
    
    releaf_trees = releaf_output_dir / "new_trees"
    if releaf_trees.exists():
        for tree_file in releaf_trees.glob("*.tree*"):
            # Rename: remove .addasm suffix
            new_name = tree_file.name.replace(".addasm", "")
            dest = trees_dir / new_name
            if dest.exists():
                dest.unlink()
            dest.symlink_to(tree_file.resolve())
        logger.info(f"    ✓ Linked trees from ReLeaf")
    
    # 5. IQTree info - from base version
    iqtree_source = base_orthophyl / "phylo_current" / "SpeciesTree" / "iqtree"
    if iqtree_source.exists():
        iqtree_dest = phylo_dir / "SpeciesTree" / "iqtree"
        iqtree_dest.parent.mkdir(parents=True, exist_ok=True)
        if iqtree_dest.exists():
            iqtree_dest.unlink()
        iqtree_dest.symlink_to(iqtree_source.resolve())
        logger.info(f"    ✓ Linked IQTree info from base version")
    
    # 6. SCO sets - from base version (unchanged)
    for sco_file in base_orthophyl.glob("phylo_current/SCO_*"):
        dest = phylo_dir / sco_file.name
        if dest.exists():
            dest.unlink()
        dest.symlink_to(sco_file.resolve())
    logger.info(f"    ✓ Linked SCO sets from base version")


def create_releaf_version(
    db_dir: Path,
    releaf_output: Path,
    version_name: Optional[str] = None
):
    """Create a new database version from ReLeaf output."""
    
    logger.info(f"\nCreating new ReLeaf version for database: {db_dir.name}")
    
    # Validate inputs
    if not db_dir.exists():
        raise FileNotFoundError(f"Database not found: {db_dir}")
    
    # Validate ReLeaf output
    validation = validate_releaf_output(releaf_output)
    if not validation['valid']:
        raise ValueError(f"Invalid ReLeaf output: {releaf_output}\nWarnings: {validation['warnings']}")
    
    # Get current version
    current_version = get_current_version(db_dir)
    logger.info(f"  Current version: {current_version}")
    
    # Load current version info
    current_version_dir = db_dir / current_version
    with open(current_version_dir / "version_info.json", 'r') as f:
        parent_info = json.load(f)
    
    # Determine new version name
    if not version_name:
        version_num = get_next_version_number(db_dir)
        version_name = f"v{version_num}_releaf_{datetime.now().strftime('%Y-%m-%d')}"
    
    new_version_dir = db_dir / version_name
    if new_version_dir.exists():
        raise FileExistsError(f"Version already exists: {new_version_dir}")
    
    new_version_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"  New version: {version_name}")
    
    # Count genomes
    new_genomes, total_genomes = count_genomes_in_releaf(releaf_output)
    if total_genomes is None:
        total_genomes = parent_info['total_genomes'] + new_genomes
    
    # Create version info
    version_info = {
        "version": version_name,
        "version_number": get_next_version_number(db_dir),
        "version_type": "releaf",
        "created": datetime.now().isoformat(),
        "parent_version": current_version,
        "base_genomes": parent_info['total_genomes'],
        "added_genomes": new_genomes,
        "total_genomes": total_genomes,
        "source_type": "ReLeaf",
        "source_path": str(releaf_output.resolve())
    }
    
    with open(new_version_dir / "version_info.json", 'w') as f:
        json.dump(version_info, f, indent=2)
    
    logger.info(f"  ✓ Created version info")
    logger.info(f"    Base genomes: {parent_info['total_genomes']}")
    logger.info(f"    New genomes: {new_genomes}")
    logger.info(f"    Total genomes: {total_genomes}")
    
    # Create symlink to parent version
    parent_link = new_version_dir / "base_version"
    parent_link.symlink_to(f"../{current_version}")
    
    # Create symlink to ReLeaf output
    releaf_link = new_version_dir / "releaf_output"
    releaf_link.symlink_to(releaf_output.resolve())
    
    # Create composite orthophyl_run directory
    composite_dir = new_version_dir / "orthophyl_run_composite"
    create_composite_orthophyl_dir(current_version_dir, releaf_output, composite_dir)
    
    # Create orthophyl_run symlink pointing to composite
    orthophyl_link = new_version_dir / "orthophyl_run"
    orthophyl_link.symlink_to("orthophyl_run_composite")
    
    # Copy database config from parent and update
    with open(current_version_dir / "database_config.json", 'r') as f:
        config = json.load(f)
    
    config['n_genomes'] = total_genomes
    config['last_updated'] = datetime.now().isoformat()
    config['version'] = version_name
    
    with open(new_version_dir / "database_config.json", 'w') as f:
        json.dump(config, f, indent=2)
    
    # Link to best tree from ReLeaf
    tree_files = list((releaf_output / "new_trees").glob("*.tree*"))
    if tree_files:
        # Prefer iqtree.CDS
        best_tree = None
        for t in tree_files:
            if 'iqtree' in t.name.lower() and 'cds' in t.name.lower():
                best_tree = t
                break
        if not best_tree:
            best_tree = tree_files[0]
        
        phylogeny_link = new_version_dir / "phylogeny.nwk"
        phylogeny_link.symlink_to(best_tree.resolve())
        logger.info(f"  ✓ Linked phylogeny: {best_tree.name}")
    
    # Update genome list
    genome_list_file = new_version_dir / "genome_list.txt"
    if (releaf_output / "genome_list").exists():
        shutil.copy(releaf_output / "genome_list", genome_list_file)
    else:
        # Merge parent genome list with new genomes
        with open(current_version_dir / "genome_list.txt", 'r') as f:
            old_genomes = set(line.strip() for line in f if line.strip() and not line.startswith('#'))
        
        # Get new genomes from annotations
        annots_dir = releaf_output / "annots_prots"
        if annots_dir.exists():
            new_genome_ids = [f.stem for f in annots_dir.glob("*.faa")]
            all_genomes = sorted(old_genomes | set(new_genome_ids))
        else:
            all_genomes = sorted(old_genomes)
        
        with open(genome_list_file, 'w') as f:
            f.write(f"# Genomes in {db_dir.name} {version_name}\n")
            f.write(f"# Created: {datetime.now().isoformat()}\n")
            f.write(f"# Parent: {current_version} ({parent_info['total_genomes']} genomes)\n")
            f.write(f"# Added: {new_genomes} genomes\n")
            for genome in all_genomes:
                f.write(f"{genome}\n")
    
    logger.info(f"  ✓ Created genome list")
    
    # Update current symlink
    current_link = db_dir / "current"
    if current_link.exists():
        current_link.unlink()
    current_link.symlink_to(version_name)
    logger.info(f"  ✓ Updated 'current' to {version_name}")
    
    logger.info(f"\n✓ Successfully created version {version_name}")
    logger.info(f"  Location: {new_version_dir}")
    logger.info(f"  Composite structure: {composite_dir}")
    
    return new_version_dir


def list_versions(db_dir: Path):
    """List all versions in a database."""
    if not db_dir.exists():
        raise FileNotFoundError(f"Database not found: {db_dir}")
    
    versions = sorted([d for d in db_dir.glob("v*") if d.is_dir()])
    
    if not versions:
        logger.info(f"No versions found in {db_dir}")
        return
    
    current = get_current_version(db_dir)
    
    print(f"\nVersions in {db_dir.name}:")
    print("=" * 70)
    
    for v_dir in versions:
        version_info_file = v_dir / "version_info.json"
        if not version_info_file.exists():
            continue
        
        with open(version_info_file, 'r') as f:
            info = json.load(f)
        
        is_current = " (current)" if v_dir.name == current else ""
        print(f"\n{info['version']}{is_current}")
        print(f"  Type: {info['version_type']}")
        print(f"  Created: {info['created']}")
        print(f"  Genomes: {info['total_genomes']} ({info['added_genomes']} added)")
        if info['parent_version']:
            print(f"  Parent: {info['parent_version']}")


def main():
    parser = argparse.ArgumentParser(
        description="Add ReLeaf output as new database version",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--database-dir',
        required=True,
        help='Database directory (e.g., databases/Rhizobiaceae_db/)'
    )
    parser.add_argument(
        '--releaf-output',
        help='ReLeaf output directory (e.g., /path/to/store/ReLeaf_dir/)'
    )
    parser.add_argument(
        '--version-name',
        help='Custom version name (default: auto-generated)'
    )
    parser.add_argument(
        '--list-versions',
        action='store_true',
        help='List all versions in the database'
    )
    
    args = parser.parse_args()
    
    db_dir = Path(args.database_dir)
    
    try:
        if args.list_versions:
            list_versions(db_dir)
        elif args.releaf_output:
            create_releaf_version(
                db_dir=db_dir,
                releaf_output=Path(args.releaf_output),
                version_name=args.version_name
            )
        else:
            parser.error("Either --releaf-output or --list-versions required")
        
        return 0
        
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())