#!/usr/bin/env python3
"""
Create Hierarchical Taxonomy Database from OrthoPhyl Runs

This script creates a queryable database from one or more OrthoPhyl output directories.
Each OrthoPhyl run represents a single taxonomic clade (e.g., order, family, genus).

The database allows assembly_router_hierarchical.py to determine which OrthoPhyl
database an assembly belongs to based on hierarchical taxonomy matching.

Input Format (TSV):
    clade_name<TAB>orthophyl_dir<TAB>clade_taxonomy
    
    Example:
    Gaiellales	/path/to/gaiellales_run	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales
    Escherichia	/path/to/ecoli_run	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia

Usage:
    python create_hierarchical_database.py \
        --input orthophyl_runs.tsv \
        --output-dir taxonomy_databases/

Author: Generated for HGTool
Date: 2024
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
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
        'warnings': []
    }
    
    # Check for HMM profiles (required for ReLeaf)
    hmm_locations = [
        orthophyl_dir / "OG_alignmentsToHMM" / "hmms_final",
        orthophyl_dir / "annots_prots.fixed" / "OrthoFinder" / "Results_ortho" / "MultipleSequenceAlignments"
    ]
    
    for hmm_dir in hmm_locations:
        if hmm_dir.exists():
            hmm_files = list(hmm_dir.glob("*.hmm")) + list(hmm_dir.glob("*.fa"))
            if hmm_files:
                validation['has_hmms'] = True
                validation['hmm_dir'] = hmm_dir
                logger.info(f"  Found {len(hmm_files)} HMM/alignment files in {hmm_dir.relative_to(orthophyl_dir)}")
                break
    
    if not validation['has_hmms']:
        validation['warnings'].append("No HMM profiles found - ReLeaf will not work")
    
    # Check for alignments
    alignment_locations = [
        orthophyl_dir / "phylo_current" / "AlignmentsProts.trm",
        orthophyl_dir / "phylo_current" / "AlignmentsCDS.trm"
    ]
    
    for aln_dir in alignment_locations:
        if aln_dir.exists():
            validation['has_alignments'] = True
            break
    
    # Check for species trees
    tree_dir = orthophyl_dir / "FINAL_SPECIES_TREES"
    if tree_dir.exists():
        tree_files = list(tree_dir.glob("*.tree"))
        if tree_files:
            validation['has_trees'] = True
            # Prefer IQ-TREE output
            for tree_file in tree_files:
                if 'iqtree' in tree_file.name.lower():
                    validation['tree_file'] = tree_file
                    break
            if not validation['tree_file']:
                validation['tree_file'] = tree_files[0]
            logger.info(f"  Found tree: {validation['tree_file'].name}")
    
    if not validation['has_trees']:
        validation['warnings'].append("No species tree found")
    
    # Count genomes
    genome_list = orthophyl_dir / "genome_list"
    all_input_list = orthophyl_dir / "all_input_list"
    
    if genome_list.exists():
        with open(genome_list, 'r') as f:
            validation['n_genomes'] = len([l for l in f if l.strip()])
    elif all_input_list.exists():
        with open(all_input_list, 'r') as f:
            validation['n_genomes'] = len([l for l in f if l.strip()])
    else:
        # Try to count from protein files
        prot_dir = orthophyl_dir / "annots_prots"
        if prot_dir.exists():
            validation['n_genomes'] = len(list(prot_dir.glob("*.faa")))
    
    logger.info(f"  Estimated genomes: {validation['n_genomes']}")
    
    # Overall validation
    validation['valid'] = validation['has_hmms'] or validation['has_alignments']
    
    return validation


def create_database_for_run(
    orthophyl_dir: Path,
    clade_taxonomy: str,
    clade_name: str,
    output_dir: Path
) -> Path:
    """
    Create a database directory for a single OrthoPhyl run.
    
    Returns:
        Path to created database directory
    """
    logger.info(f"Creating database for: {clade_name}")
    logger.info(f"  Taxonomy: {clade_taxonomy}")
    logger.info(f"  Source: {orthophyl_dir}")
    
    # Validate OrthoPhyl run
    validation = validate_orthophyl_run(orthophyl_dir)
    
    if not validation['valid']:
        raise ValueError(f"Invalid OrthoPhyl run: {orthophyl_dir}")
    
    if validation['warnings']:
        for warning in validation['warnings']:
            logger.warning(f"  Warning: {warning}")
    
    # Parse taxonomy to determine clade rank
    tax = GTDBTaxonomy(clade_taxonomy)
    clade_rank = tax.get_most_specific_rank()
    
    if not clade_rank:
        raise ValueError(f"Cannot determine clade rank from taxonomy: {clade_taxonomy}")
    
    logger.info(f"  Clade rank: {tax.get_rank_name(clade_rank)}")
    
    # Create database directory
    # Use safe filename from clade name
    safe_name = clade_name.replace(' ', '_').replace('/', '_')
    db_dir = output_dir / f"{safe_name}_db"
    db_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"  Database directory: {db_dir}")
    
    # Create database_config.json
    config = {
        "created": datetime.now().isoformat(),
        "version": "1.0",
        "database_type": "hierarchical_taxonomy",
        "clade_name": clade_name,
        "clade_taxonomy": clade_taxonomy,
        "clade_rank": clade_rank,
        "clade_rank_name": tax.get_rank_name(clade_rank),
        "orthophyl_source": str(orthophyl_dir.resolve()),
        "n_genomes": validation['n_genomes'],
        "has_hmms": validation['has_hmms'],
        "has_trees": validation['has_trees']
    }
    
    config_file = db_dir / "database_config.json"
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    logger.info(f"  Created config: {config_file.name}")
    
    # Copy or link phylogeny
    if validation['tree_file']:
        dest_tree = db_dir / "phylogeny.nwk"
        shutil.copy(validation['tree_file'], dest_tree)
        logger.info(f"  Copied tree: {validation['tree_file'].name} → {dest_tree.name}")
    else:
        # Create placeholder
        (db_dir / "phylogeny.nwk").write_text(f"# Placeholder tree for {clade_name}\n")
        logger.info("  Created placeholder tree")
    
    # Create symlink to OrthoPhyl run
    orthophyl_link = db_dir / "orthophyl_run"
    if orthophyl_link.exists():
        orthophyl_link.unlink()
    
    try:
        orthophyl_link.symlink_to(orthophyl_dir.resolve())
        logger.info(f"  Created symlink: orthophyl_run → {orthophyl_dir}")
    except OSError as e:
        logger.warning(f"  Could not create symlink: {e}")
        logger.info("  Creating directory and copying key files...")
        orthophyl_link.mkdir(exist_ok=True)
        
        # Copy essential files/directories
        if validation['hmm_dir']:
            dest_hmm = orthophyl_link / "hmms"
            shutil.copytree(validation['hmm_dir'], dest_hmm, dirs_exist_ok=True)
    
    # Create genome list
    genome_list_file = db_dir / "genome_list.txt"
    genome_names = []
    
    # Try to extract genome names from various sources
    for source in ['genome_list', 'all_input_list']:
        source_file = orthophyl_dir / source
        if source_file.exists():
            with open(source_file, 'r') as f:
                genome_names = [line.strip() for line in f if line.strip()]
            break
    
    if not genome_names:
        # Try from protein directory
        prot_dir = orthophyl_dir / "annots_prots"
        if prot_dir.exists():
            genome_names = [f.stem for f in prot_dir.glob("*.faa")]
    
    with open(genome_list_file, 'w') as f:
        f.write(f"# Genomes in {clade_name} database\n")
        f.write(f"# Created: {datetime.now().isoformat()}\n")
        for genome in sorted(genome_names):
            f.write(f"{genome}\n")
    
    logger.info(f"  Created genome list: {len(genome_names)} genomes")
    
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
        
        f.write("Source OrthoPhyl Run:\n")
        f.write(f"  {orthophyl_dir}\n\n")
        
        f.write("Use with assembly_router_hierarchical.py:\n")
        f.write(f"  python assembly_router_hierarchical.py \\\n")
        f.write(f"      --assembly genome.fna \\\n")
        f.write(f"      --taxonomy \"d__...; ...\" \\\n")
        f.write(f"      --database-dir {db_dir}\n\n")
        
        if validation['warnings']:
            f.write("Warnings:\n")
            for warning in validation['warnings']:
                f.write(f"  - {warning}\n")
    
    logger.info(f"  Created README")
    logger.info(f"✓ Database created successfully\n")
    
    return db_dir


def create_master_index(databases: List[Dict], output_dir: Path):
    """
    Create master index of all databases for easy lookup.
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
        
        index["databases"].append({
            "clade_name": db['clade_name'],
            "clade_taxonomy": db['clade_taxonomy'],
            "clade_rank": tax.get_most_specific_rank(),
            "clade_rank_name": tax.get_rank_name(tax.get_most_specific_rank()),
            "database_dir": str(db['database_dir']),
            "n_genomes": db['n_genomes']
        })
    
    with open(index_file, 'w') as f:
        json.dump(index, f, indent=2)
    
    logger.info(f"Created master index: {index_file}")
    
    # Also create human-readable summary
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
        
        f.write("Route an assembly to appropriate database:\n\n")
        f.write("  python assembly_router_hierarchical.py \\\n")
        f.write("      --assembly genome.fna \\\n")
        f.write("      --taxonomy \"d__Bacteria;p__...; ...\" \\\n")
        f.write("      --database-dir [one of the database directories above]\n\n")
    
    logger.info(f"Created summary: {summary_file}")


def parse_input_table(input_file: Path) -> List[Dict]:
    """
    Parse input TSV file.
    
    Format: clade_name<TAB>orthophyl_dir<TAB>clade_taxonomy
    
    Returns:
        List of dictionaries with clade_name, orthophyl_dir, clade_taxonomy
    """
    runs = []
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue
            
            # Skip header if present
            if line_num == 1 and 'clade' in line.lower():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                logger.warning(f"Line {line_num}: expected 3 fields (clade_name, orthophyl_dir, clade_taxonomy), got {len(fields)}. Skipping.")
                continue
            
            clade_name = fields[0].strip()
            orthophyl_dir = Path(fields[1].strip())
            clade_taxonomy = fields[2].strip()
            
            if not clade_name:
                logger.warning(f"Line {line_num}: empty clade_name, skipping")
                continue
            
            if not clade_taxonomy:
                logger.warning(f"Line {line_num}: empty clade_taxonomy, skipping")
                continue
            
            runs.append({
                'clade_name': clade_name,
                'orthophyl_dir': orthophyl_dir,
                'clade_taxonomy': clade_taxonomy
            })
    
    logger.info(f"Parsed {len(runs)} OrthoPhyl runs from {input_file}")
    return runs


def main():
    parser = argparse.ArgumentParser(
        description="Create hierarchical taxonomy databases from OrthoPhyl runs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create databases from input table
  python create_hierarchical_database.py \\
      --input orthophyl_runs.tsv \\
      --output-dir taxonomy_databases/

Input Table Format (TSV):
  clade_name<TAB>orthophyl_dir<TAB>clade_taxonomy
  
  Example:
  Gaiellales	/data/gaiellales_run	d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales
  Escherichia	/data/ecoli_run	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia
  Actinomycetota	/data/actino_run	d__Bacteria;p__Actinomycetota

Output:
  taxonomy_databases/
  ├── Gaiellales_db/
  │   ├── database_config.json
  │   ├── phylogeny.nwk
  │   ├── orthophyl_run -> /data/gaiellales_run
  │   ├── genome_list.txt
  │   └── README.txt
  ├── Escherichia_db/
  │   └── ...
  ├── database_index.json
  └── database_summary.txt
        """
    )
    
    parser.add_argument(
        '--input',
        required=True,
        help='TSV file: clade_name, orthophyl_dir, clade_taxonomy'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for databases'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing databases'
    )
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("Creating Hierarchical Taxonomy Databases")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {output_dir}")
    logger.info("")
    
    # Parse input
    try:
        runs = parse_input_table(Path(args.input))
    except Exception as e:
        logger.error(f"Failed to parse input: {e}")
        return 1
    
    if not runs:
        logger.error("No valid OrthoPhyl runs found in input")
        return 1
    
    # Create databases
    databases = []
    failed = []
    
    for i, run in enumerate(runs, 1):
        logger.info(f"[{i}/{len(runs)}] Processing: {run['clade_name']}")
        
        try:
            db_dir = create_database_for_run(
                run['orthophyl_dir'],
                run['clade_taxonomy'],
                run['clade_name'],
                output_dir
            )
            
            # Read back config to get accurate counts
            with open(db_dir / "database_config.json", 'r') as f:
                config = json.load(f)
            
            databases.append({
                'clade_name': run['clade_name'],
                'clade_taxonomy': run['clade_taxonomy'],
                'database_dir': db_dir,
                'n_genomes': config['n_genomes']
            })
            
        except Exception as e:
            logger.error(f"Failed to create database for {run['clade_name']}: {e}")
            failed.append(run['clade_name'])
    
    # Create master index
    if databases:
        logger.info("")
        logger.info("Creating master index...")
        create_master_index(databases, output_dir)
    
    # Final summary
    logger.info("")
    logger.info("=" * 70)
    logger.info("DATABASE CREATION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Successfully created: {len(databases)} databases")
    
    if failed:
        logger.info(f"Failed: {len(failed)} databases")
        for name in failed:
            logger.info(f"  - {name}")
    
    logger.info(f"\nDatabases location: {output_dir}")
    logger.info(f"Master index: {output_dir / 'database_index.json'}")
    logger.info(f"Summary: {output_dir / 'database_summary.txt'}")
    
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
