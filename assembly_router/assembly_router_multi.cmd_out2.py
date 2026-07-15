#!/usr/bin/env python3
"""
Multi-Database Assembly Router - Query multiple databases and find best match

This script queries all available databases in a directory and routes assemblies
to the best matching database (ReLeaf) or suggests creating a new OrthoPhyl run
for novel taxa with commands using gather_filter_asms.sh.

Usage:
    python assembly_router_multi.py \
        --assembly genome.fna \
        --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \
        --database-dir /path/to/databases/ \
        --gather-filter-script utils/gather_filter_asms.sh \
        --output-dir results/

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
import re

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class GTDBTaxonomy:
    """Parse and manipulate GTDB taxonomy strings."""
    
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
        return self.levels.get(rank)
    
    def get_most_specific_rank(self) -> Optional[str]:
        for rank in reversed(self.RANKS):
            if self.levels.get(rank):
                return rank
        return None
    
    def get_rank_name(self, rank: str) -> str:
        return self.RANK_NAMES.get(rank, rank)
    
    def is_within_clade(self, clade: 'GTDBTaxonomy', clade_rank: str) -> bool:
        """Check if this taxonomy belongs to specified clade."""
        for rank in self.RANKS:
            if rank == clade_rank:
                return self.get_rank(rank) == clade.get_rank(rank)
            
            rank_idx = self.RANKS.index(rank)
            clade_idx = self.RANKS.index(clade_rank)
            
            if rank_idx < clade_idx:
                if self.get_rank(rank) != clade.get_rank(rank):
                    return False
        return True
    
    def get_taxonomy_string_at_rank(self, rank: str) -> str:
        """Get taxonomy string up to specified rank."""
        parts = []
        for r in self.RANKS:
            value = self.get_rank(r)
            parts.append(f"{r}__{value if value else ''}")
            if r == rank:
                break
        return ';'.join(parts)
    
    def __str__(self):
        return self.raw


class MultiDatabaseRouter:
    """Routes assemblies by querying multiple databases."""
    
    def __init__(
        self,
        database_dir: Path,
        output_dir: Path,
        gather_filter_script: Optional[Path] = None,
        threads: int = 8
    ):
        self.database_dir = Path(database_dir)
        self.output_dir = Path(output_dir)
        self.gather_filter_script = Path(gather_filter_script) if gather_filter_script else None
        self.threads = threads
        
        # Load all databases
        self.databases = self._load_databases()
        
        if not self.databases:
            raise ValueError(f"No databases found in {database_dir}")
        
        logger.info(f"Loaded {len(self.databases)} databases")
        for db in self.databases:
            logger.info(f"  - {db['clade_name']} ({db['rank_name']}, {db['n_genomes']} genomes)")
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup file logging
        log_file = self.output_dir / f"routing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logger.addHandler(file_handler)
    
    def _load_databases(self) -> List[Dict]:
        """Load all database configs from directory."""
        databases = []
        
        # Check for database_index.json
        index_file = self.database_dir / "database_index.json"
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)
            for db_info in index['databases']:
                db_dir = Path(db_info['database_dir'])
                if db_dir.exists():
                    databases.append(self._load_single_database(db_dir))
        else:
            # Scan for *_db directories
            for db_dir in self.database_dir.glob("*_db"):
                config_file = db_dir / "database_config.json"
                if config_file.exists():
                    databases.append(self._load_single_database(db_dir))
        
        return databases
    
    def _load_single_database(self, db_dir: Path) -> Dict:
        """Load a single database config."""
        config_file = db_dir / "database_config.json"
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        tax = GTDBTaxonomy(config['clade_taxonomy'])
        
        return {
            'db_dir': db_dir,
            'clade_name': config['clade_name'],
            'clade_taxonomy': config['clade_taxonomy'],
            'clade_rank': config['clade_rank'],
            'rank_name': config['clade_rank_name'],
            'n_genomes': config['n_genomes'],
            'taxonomy_obj': tax
        }
    
    def find_matching_databases(self, query_taxonomy: str) -> List[Dict]:
        """Find all databases that match the query taxonomy."""
        query_tax = GTDBTaxonomy(query_taxonomy)
        matches = []
        
        for db in self.databases:
            clade_tax = db['taxonomy_obj']
            clade_rank = db['clade_rank']
            
            if query_tax.is_within_clade(clade_tax, clade_rank):
                matches.append({
                    'database': db,
                    'specificity': GTDBTaxonomy.RANKS.index(clade_rank)  # Higher = more specific
                })
        
        # Sort by specificity (most specific first)
        matches.sort(key=lambda x: x['specificity'], reverse=True)
        
        return matches
    
    def route_assembly(
        self,
        assembly_path: Path,
        taxonomy: str,
        assembly_id: Optional[str] = None
    ) -> Dict:
        """Route assembly to best matching database or suggest new OrthoPhyl run."""
        if assembly_id is None:
            assembly_id = Path(assembly_path).stem
        
        logger.info(f"\nRouting assembly: {assembly_id}")
        logger.info(f"Query taxonomy: {taxonomy}")
        
        # Find matching databases
        matches = self.find_matching_databases(taxonomy)
        
        if matches:
            # Use most specific match
            best_match = matches[0]['database']
            
            logger.info(f"✓ MATCH FOUND: {best_match['clade_name']}")
            logger.info(f"  Rank: {best_match['rank_name']}")
            logger.info(f"  Database: {best_match['db_dir'].name}")
            logger.info(f"  Genomes: {best_match['n_genomes']}")
            
            if len(matches) > 1:
                logger.info(f"  (Also matches {len(matches)-1} other database(s) at broader levels)")
            
            return self._route_to_releaf(assembly_path, assembly_id, taxonomy, best_match)
        else:
            logger.info("✗ NO MATCH FOUND in any database")
            logger.info("  → New OrthoPhyl run required")
            
            return self._route_to_orthophyl(assembly_path, assembly_id, taxonomy)
    
    def _route_to_releaf(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str,
        database: Dict
    ) -> Dict:
        """Generate ReLeaf routing decision."""
        # Prepare ReLeaf input
        releaf_input_dir = self.output_dir / "releaf_input"
        releaf_input_dir.mkdir(parents=True, exist_ok=True)
        
        import shutil
        dest = releaf_input_dir / f"{assembly_id}.fna"
        if not dest.exists():
            shutil.copy(assembly_path, dest)
        
        # Build ReLeaf command
        releaf_cmd = (
            f"./ReLeaf.sh \\\n"
            f"    --store {database['db_dir'] / 'orthophyl_run'} \\\n"
            f"    --input_genomes {releaf_input_dir} \\\n"
            f"    -t {self.threads} \\\n"
            f"    --tree_method iqtree \\\n"
            f"    --TREE_DATA CDS"
        )
        
        decision = {
            'pipeline': 'ReLeaf',
            'reason': f"Taxonomy matches {database['clade_name']} at {database['rank_name']} level",
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': taxonomy,
            'matched_database': database['clade_name'],
            'matched_rank': database['rank_name'],
            'database_dir': str(database['db_dir']),
            'database_genomes': database['n_genomes'],
            'command': releaf_cmd
        }
        
        self._save_decision(decision)
        return decision
    
    def _route_to_orthophyl(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str
    ) -> Dict:
        """Generate OrthoPhyl routing decision with gather_filter_asms.sh."""
        query_tax = GTDBTaxonomy(taxonomy)
        specific_rank = query_tax.get_most_specific_rank()
        
        # Determine what taxonomic level to download
        download_rank = specific_rank
        download_value = query_tax.get_rank(specific_rank) if specific_rank else None
        
        # If at species level and species is undefined, go up to genus
        if specific_rank == 's' and not download_value:
            download_rank = 'g'
            download_value = query_tax.get_rank('g')
        
        # If still no value, go up to family
        if not download_value and download_rank in ['s', 'g']:
            download_rank = 'f'
            download_value = query_tax.get_rank('f')
        
        rank_name = query_tax.get_rank_name(download_rank) if download_rank else 'unknown'
        
        # Prepare OrthoPhyl input directory
        orthophyl_input_dir = self.output_dir / "orthophyl_new_genomes"
        orthophyl_input_dir.mkdir(parents=True, exist_ok=True)
        
        import shutil
        dest = orthophyl_input_dir / Path(assembly_path).name
        if not dest.exists():
            shutil.copy(assembly_path, dest)
        
        # Build taxonomy string for gather_filter_asms.sh
        # The script can use taxon name directly
        
        # Build commands
        commands = []
        
        # Step 1: Download and filter related genomes using gather_filter_asms.sh
        if self.gather_filter_script and self.gather_filter_script.exists():
            # gather_filter_asms.sh usage: script taxon output_dir threads
            # The taxon can be a taxonomy string or taxon name
            gather_cmd = (
                f"# Step 1: Download and filter related genomes from NCBI\n"
                f"# This will download genomes, run CheckM QC, and filter by quality\n"
                f"{self.gather_filter_script} \\\n"
                f"    {download_value} \\\n"
                f"    {orthophyl_input_dir} \\\n"
                f"    {self.threads}\n"
                f"\n"
                f"# Note: The script will create {orthophyl_input_dir}/genomes_to_keep/\n"
                f"# with high-quality, non-redundant genomes.\n"
                f"# Quality filters (adjust in script if needed):\n"
                f"#   - Completeness >= 95%\n"
                f"#   - Contamination <= 1.0%\n"
                f"#   - Duplication <= 2%\n"
                f"#   - Removes RefSeq/GenBank redundancy\n"
                f"#   - Filters by assembly stats (N50, GC content, length)"
            )
            commands.append(gather_cmd)
            
            # Update orthophyl input directory to the filtered genomes
            orthophyl_input_actual = f"{orthophyl_input_dir}/genomes_to_keep"
        else:
            gather_cmd = (
                f"# Step 1: Download related genomes manually\n"
                f"# Search NCBI for: {download_value} ({rank_name})\n"
                f"# Recommended: Use utils/gather_filter_asms.sh for automatic download & QC\n"
                f"# Place genome files in: {orthophyl_input_dir}/"
            )
            commands.append(gather_cmd)
            orthophyl_input_actual = str(orthophyl_input_dir)
        
        # Step 2: Run OrthoPhyl
        orthophyl_cmd = (
            f"# Step 2: Run OrthoPhyl on expanded genome set\n"
            f"# Note: Your query genome ({Path(assembly_path).name}) should be in {orthophyl_input_actual}\n"
            f"./OrthoPhyl.sh \\\n"
            f"    -g {orthophyl_input_actual} \\\n"
            f"    -o {self.output_dir / 'orthophyl_output'} \\\n"
            f"    -t {self.threads} \\\n"
            f"    --tree_method iqtree \\\n"
            f"    --TREE_DATA CDS \\\n"
            f"    --use_partitions true"
        )
        commands.append(orthophyl_cmd)
        
        # Step 3: Create database
        tax_string_for_download = query_tax.get_taxonomy_string_at_rank(download_rank)
        db_cmd = (
            f"# Step 3: Create new database entry\n"
            f"# Add this line to your orthophyl_runs.tsv:\n"
            f"{download_value}\t{self.output_dir / 'orthophyl_output'}\t{tax_string_for_download}\n"
            f"\n"
            f"# Then rebuild the database index:\n"
            f"python create_hierarchical_database_v2.py \\\n"
            f"    --input orthophyl_runs.tsv \\\n"
            f"    --output-dir {self.database_dir} \\\n"
            f"    --update"
        )
        commands.append(db_cmd)
        
        full_command = "\n\n".join(commands)
        
        decision = {
            'pipeline': 'OrthoPhyl',
            'reason': 'Novel taxonomy not represented in any database',
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': taxonomy,
            'download_taxonomy': tax_string_for_download,
            'download_rank': rank_name,
            'download_value': download_value,
            'suggestion': f"Create new database for {rank_name} '{download_value}'",
            'command': full_command
        }
        
        self._save_decision(decision)
        return decision
    
    def _save_decision(self, decision: Dict):
        """Save routing decision to files."""
        decision_file = self.output_dir / f"routing_decision_{decision['assembly_id']}.json"
        
        with open(decision_file, 'w') as f:
            json.dump(decision, f, indent=2)
        
        logger.info(f"\n→ Decision saved: {decision_file}")
        
        # Human-readable summary
        summary_file = self.output_dir / f"routing_summary_{decision['assembly_id']}.txt"
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("ASSEMBLY ROUTING DECISION (Multi-Database Query)\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Assembly ID: {decision['assembly_id']}\n")
            f.write(f"Assembly: {decision['assembly']}\n")
            f.write(f"Query Taxonomy: {decision['query_taxonomy']}\n")
            f.write(f"\nPipeline: {decision['pipeline']}\n")
            f.write(f"Reason: {decision['reason']}\n\n")
            
            if decision['pipeline'] == 'ReLeaf':
                f.write("Match Details:\n")
                f.write(f"  Database: {decision['matched_database']}\n")
                f.write(f"  Matched at: {decision['matched_rank']} level\n")
                f.write(f"  Database directory: {decision['database_dir']}\n")
                f.write(f"  Contains: {decision['database_genomes']} genomes\n")
            else:
                f.write("Action Required:\n")
                f.write(f"  {decision['suggestion']}\n")
                f.write(f"  Download taxonomy: {decision['download_taxonomy']}\n")
                f.write(f"  Target rank: {decision['download_rank']}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("COMMAND(S) TO RUN\n")
            f.write("=" * 70 + "\n\n")
            f.write(decision['command'] + "\n")
        
        logger.info(f"→ Summary saved: {summary_file}")
    
    def batch_route(self, input_table: Path) -> List[Dict]:
        """Route multiple assemblies from table."""
        logger.info(f"\nBatch routing from {input_table}")
        
        decisions = []
        
        with open(input_table, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('#') and not first_line.lower().startswith('assembly'):
                f.seek(0)
            
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    logger.warning(f"Line {line_num}: insufficient fields")
                    continue
                
                assembly_path = Path(fields[0]).expanduser()
                taxonomy = fields[1].strip('"')  # Remove quotes if present
                assembly_id = fields[2] if len(fields) > 2 else None
                
                try:
                    decision = self.route_assembly(assembly_path, taxonomy, assembly_id)
                    decisions.append(decision)
                except Exception as e:
                    logger.error(f"Failed to route {assembly_path}: {e}")
        
        self._generate_batch_summary(decisions)
        return decisions
    
    def _generate_batch_summary(self, decisions: List[Dict]):
        """Generate batch summary."""
        summary_file = self.output_dir / "batch_routing_summary.txt"
        
        n_releaf = sum(1 for d in decisions if d['pipeline'] == 'ReLeaf')
        n_orthophyl = sum(1 for d in decisions if d['pipeline'] == 'OrthoPhyl')
        
        # Group by database
        by_database = {}
        for d in decisions:
            if d['pipeline'] == 'ReLeaf':
                db = d['matched_database']
                if db not in by_database:
                    by_database[db] = []
                by_database[db].append(d['assembly_id'])
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("BATCH ROUTING SUMMARY (Multi-Database)\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Available databases: {len(self.databases)}\n")
            for db in self.databases:
                f.write(f"  - {db['clade_name']} ({db['rank_name']}, {db['n_genomes']} genomes)\n")
            f.write("\n")
            
            f.write(f"Total assemblies: {len(decisions)}\n")
            f.write(f"  → ReLeaf (matched): {n_releaf}\n")
            f.write(f"  → OrthoPhyl (novel): {n_orthophyl}\n\n")
            
            if n_releaf > 0:
                f.write("ReLeaf Routing by Database:\n")
                f.write("-" * 70 + "\n")
                for db_name, assemblies in sorted(by_database.items()):
                    f.write(f"\n{db_name} ({len(assemblies)} assemblies):\n")
                    for asm in assemblies:
                        f.write(f"  - {asm}\n")
                f.write("\n")
            
            if n_orthophyl > 0:
                f.write("OrthoPhyl Assemblies (novel taxa):\n")
                f.write("-" * 70 + "\n")
                for d in decisions:
                    if d['pipeline'] == 'OrthoPhyl':
                        f.write(f"  {d['assembly_id']:30s} - {d['download_rank']}: {d['download_value']}\n")
                f.write("\n")
        
        logger.info(f"\n→ Batch summary: {summary_file}")
        print(f"\nRouting complete: {n_releaf} → ReLeaf, {n_orthophyl} → OrthoPhyl")


def main():
    parser = argparse.ArgumentParser(
        description="Multi-database assembly router with automatic best-match selection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Route single assembly (queries all databases)
  python assembly_router_multi.py \\
      --assembly genome.fna \\
      --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \\
      --database-dir /path/to/databases/ \\
      --output-dir results/
  
  # With gather_filter_asms.sh for downloading genomes
  python assembly_router_multi.py \\
      --assembly genome.fna \\
      --taxonomy "d__Bacteria;..." \\
      --database-dir /path/to/databases/ \\
      --gather-filter-script utils/gather_filter_asms.sh \\
      --output-dir results/
  
  # Batch mode
  python assembly_router_multi.py \\
      --batch assemblies.tsv \\
      --database-dir /path/to/databases/ \\
      --gather-filter-script utils/gather_filter_asms.sh \\
      --output-dir results/
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--assembly', help='Single assembly FASTA')
    input_group.add_argument('--batch', help='TSV: assembly_path, taxonomy, [id]')
    
    parser.add_argument('--taxonomy', help='GTDB taxonomy (for --assembly)')
    parser.add_argument('--assembly-id', help='Assembly identifier')
    parser.add_argument(
        '--database-dir',
        required=True,
        help='Directory containing multiple *_db databases'
    )
    parser.add_argument(
        '--output-dir',
        default='routing_output',
        help='Output directory'
    )
    parser.add_argument(
        '--gather-filter-script',
        help='Path to gather_filter_asms.sh for downloading genomes'
    )
    parser.add_argument('-t', '--threads', type=int, default=8, help='Threads')
    
    args = parser.parse_args()
    
    if args.assembly and not args.taxonomy:
        parser.error("--taxonomy required with --assembly")
    
    try:
        router = MultiDatabaseRouter(
            database_dir=Path(args.database_dir),
            output_dir=Path(args.output_dir),
            gather_filter_script=Path(args.gather_filter_script) if args.gather_filter_script else None,
            threads=args.threads
        )
    except Exception as e:
        logger.error(f"Failed to initialize: {e}")
        return 1
    
    try:
        if args.batch:
            decisions = router.batch_route(Path(args.batch))
        else:
            decision = router.route_assembly(
                Path(args.assembly).expanduser(),
                args.taxonomy.strip('"'),  # Remove quotes if present
                args.assembly_id
            )
            decisions = [decision]
        
        print("\n" + "=" * 70)
        print("ROUTING COMPLETE")
        print("=" * 70)
        
        for decision in decisions:
            print(f"\n{decision['assembly_id']}:")
            print(f"  Pipeline: {decision['pipeline']}")
            if decision['pipeline'] == 'ReLeaf':
                print(f"  Database: {decision['matched_database']} ({decision['matched_rank']})")
            else:
                print(f"  Action: {decision['suggestion']}")
            print(f"\n  Commands saved to: routing_summary_{decision['assembly_id']}.txt")
        
        return 0
        
    except Exception as e:
        logger.error(f"Routing failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())