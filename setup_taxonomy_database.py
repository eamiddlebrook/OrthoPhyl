#!/usr/bin/env python3
"""
Setup Taxonomy Database - Create reference database from OrthoPhyl run

This script creates a taxonomy database from an existing OrthoPhyl run output.
The database is used by assembly_router.py to determine whether new assemblies
should be added via ReLeaf or require a new OrthoPhyl run.

Usage:
    python setup_taxonomy_database.py \
        --orthophyl-dir /path/to/orthophyl_output/ \
        --genome-metadata genomes_metadata.tsv \
        --output-dir taxonomy_db/

Genome Metadata Format (TSV):
    genome_id<TAB>taxonomy<TAB>level
    
    Example:
    GCF_000005845.2	Escherichia coli K-12	strain
    GCF_000008865.2	Escherichia coli	species
    GCF_000482265.1	Escherichia coli	species

Author: Generated for HGTool
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Set
from collections import defaultdict
from datetime import datetime
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)


def parse_genome_metadata(metadata_file: Path) -> Dict[str, Dict]:
    """
    Parse genome metadata file.
    
    Returns:
        Dict mapping genome_id to {taxonomy, level}
    """
    metadata = {}
    
    with open(metadata_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue
            
            # Skip header
            if line_num == 1 and line.lower().startswith('genome'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 2:
                logger.warning(f"Line {line_num}: insufficient fields")
                continue
            
            genome_id = fields[0]
            taxonomy = fields[1]
            level = fields[2] if len(fields) > 2 else 'strain'
            
            metadata[genome_id] = {
                'taxonomy': taxonomy,
                'level': level
            }
    
    logger.info(f"Loaded metadata for {len(metadata)} genomes")
    return metadata


def build_taxonomy_map(metadata: Dict[str, Dict]) -> Dict[str, Dict]:
    """
    Build taxonomy map from genome metadata.
    
    Groups genomes by taxonomy at different levels.
    
    Returns:
        Dict mapping taxonomy → {level, genomes}
    """
    taxonomy_map = defaultdict(lambda: {'level': None, 'genomes': []})
    
    for genome_id, info in metadata.items():
        taxonomy = info['taxonomy']
        level = info['level']
        
        # Add to exact taxonomy
        taxonomy_map[taxonomy]['genomes'].append(genome_id)
        taxonomy_map[taxonomy]['level'] = level
        
        # Also add to higher levels
        parts = taxonomy.split()
        
        # Add to species level (first two words)
        if len(parts) >= 2 and level == 'strain':
            species = ' '.join(parts[:2])
            if species != taxonomy:  # Don't duplicate
                taxonomy_map[species]['genomes'].append(genome_id)
                taxonomy_map[species]['level'] = 'species'
        
        # Add to genus level (first word)
        if len(parts) >= 1 and level in ['species', 'strain']:
            genus = parts[0]
            if genus != taxonomy:  # Don't duplicate
                taxonomy_map[genus]['genomes'].append(genome_id)
                taxonomy_map[genus]['level'] = 'genus'
    
    # Remove duplicates and sort
    for taxonomy in taxonomy_map:
        taxonomy_map[taxonomy]['genomes'] = sorted(set(taxonomy_map[taxonomy]['genomes']))
    
    return dict(taxonomy_map)


def extract_tree_from_orthophyl(orthophyl_dir: Path) -> Path:
    """
    Find the main species tree from OrthoPhyl output.
    
    Returns:
        Path to tree file
    """
    final_trees = orthophyl_dir / "FINAL_SPECIES_TREES"
    
    if not final_trees.exists():
        raise FileNotFoundError(f"Cannot find FINAL_SPECIES_TREES in {orthophyl_dir}")
    
    # Look for IQ-TREE output (preferred)
    tree_patterns = [
        "iqtree.SCO_strict.CDS.tree",
        "iqtree.SCO_strict.PROT.tree",
        "fastTree.SCO_strict.CDS.tree",
        "*.tree"
    ]
    
    for pattern in tree_patterns:
        trees = list(final_trees.glob(pattern))
        if trees:
            logger.info(f"Using tree: {trees[0]}")
            return trees[0]
    
    raise FileNotFoundError(f"No suitable tree found in {final_trees}")


def create_taxonomy_database(
    orthophyl_dir: Path,
    metadata_file: Path,
    output_dir: Path
):
    """
    Create complete taxonomy database.
    
    Args:
        orthophyl_dir: OrthoPhyl output directory
        metadata_file: Genome metadata TSV
        output_dir: Output database directory
    """
    logger.info("Creating taxonomy database...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse metadata
    metadata = parse_genome_metadata(metadata_file)
    
    # Build taxonomy map
    taxonomy_map = build_taxonomy_map(metadata)
    logger.info(f"Built taxonomy map with {len(taxonomy_map)} entries")
    
    # Write taxonomy_map.tsv
    taxonomy_map_file = output_dir / "taxonomy_map.tsv"
    with open(taxonomy_map_file, 'w') as f:
        f.write("# Taxonomy database for assembly routing\n")
        f.write("# Format: taxonomy<TAB>level<TAB>genome_ids (comma-separated)\n")
        f.write("taxonomy\tlevel\tgenomes\n")
        
        for taxonomy in sorted(taxonomy_map.keys()):
            info = taxonomy_map[taxonomy]
            genomes_str = ','.join(info['genomes'])
            f.write(f"{taxonomy}\t{info['level']}\t{genomes_str}\n")
    
    logger.info(f"Wrote taxonomy map to {taxonomy_map_file}")
    
    # Copy/link phylogeny
    try:
        tree_file = extract_tree_from_orthophyl(orthophyl_dir)
        dest_tree = output_dir / "phylogeny.nwk"
        
        # Copy tree content
        with open(tree_file, 'r') as src:
            tree_content = src.read()
        with open(dest_tree, 'w') as dst:
            dst.write(tree_content)
        
        logger.info(f"Copied phylogeny to {dest_tree}")
    except Exception as e:
        logger.warning(f"Could not extract tree: {e}")
        logger.info("Creating placeholder tree")
        (output_dir / "phylogeny.nwk").write_text("# Placeholder tree\n")
    
    # Link to OrthoPhyl run
    orthophyl_link = output_dir / "orthophyl_run"
    if orthophyl_link.exists():
        orthophyl_link.unlink()
    
    # Create symlink
    orthophyl_link.symlink_to(orthophyl_dir.resolve())
    logger.info(f"Linked OrthoPhyl run: {orthophyl_link} → {orthophyl_dir}")
    
    # Create metadata.json
    metadata_json = {
        "created": datetime.now().isoformat(),
        "version": "1.0",
        "description": "Taxonomy database for assembly routing",
        "orthophyl_dir": str(orthophyl_dir.resolve()),
        "n_genomes": len(metadata),
        "n_taxa": len(taxonomy_map),
        "taxonomic_levels": sorted(set(t['level'] for t in taxonomy_map.values())),
        "database_type": "OrthoPhyl_reference"
    }
    
    with open(output_dir / "metadata.json", 'w') as f:
        json.dump(metadata_json, f, indent=2)
    
    logger.info(f"Wrote metadata to {output_dir / 'metadata.json'}")
    
    # Generate summary
    summary_file = output_dir / "database_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("TAXONOMY DATABASE SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"Created: {metadata_json['created']}\n")
        f.write(f"OrthoPhyl run: {orthophyl_dir}\n")
        f.write(f"Total genomes: {len(metadata)}\n")
        f.write(f"Total taxonomic entries: {len(taxonomy_map)}\n\n")
        
        # Count by level
        level_counts = defaultdict(int)
        for info in taxonomy_map.values():
            level_counts[info['level']] += 1
        
        f.write("Taxonomy levels:\n")
        for level, count in sorted(level_counts.items()):
            f.write(f"  {level}: {count}\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("SAMPLE ENTRIES\n")
        f.write("=" * 70 + "\n\n")
        
        # Show first 10 entries
        for taxonomy in sorted(taxonomy_map.keys())[:10]:
            info = taxonomy_map[taxonomy]
            f.write(f"{taxonomy} ({info['level']}): {len(info['genomes'])} genomes\n")
        
        if len(taxonomy_map) > 10:
            f.write(f"\n... and {len(taxonomy_map) - 10} more\n")
    
    logger.info(f"Wrote summary to {summary_file}")
    
    print("\n" + "=" * 70)
    print("DATABASE CREATION COMPLETE")
    print("=" * 70)
    print(f"\nDatabase location: {output_dir}")
    print(f"Total genomes: {len(metadata)}")
    print(f"Total taxa: {len(taxonomy_map)}")
    print("\nUse with assembly_router.py:")
    print(f"  python assembly_router.py \\")
    print(f"      --assembly new_genome.fna \\")
    print(f"      --taxonomy \"Species name\" \\")
    print(f"      --database-dir {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Create taxonomy database from OrthoPhyl run",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--orthophyl-dir',
        required=True,
        help='OrthoPhyl output directory (contains FINAL_SPECIES_TREES/)'
    )
    parser.add_argument(
        '--genome-metadata',
        required=True,
        help='TSV file: genome_id<TAB>taxonomy<TAB>level'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output database directory'
    )
    
    args = parser.parse_args()
    
    try:
        create_taxonomy_database(
            Path(args.orthophyl_dir),
            Path(args.genome_metadata),
            Path(args.output_dir)
        )
        return 0
    except Exception as e:
        logger.error(f"Failed to create database: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())