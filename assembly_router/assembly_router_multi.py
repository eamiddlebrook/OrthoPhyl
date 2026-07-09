#!/usr/bin/env python3
"""
Multi-Database Assembly Router - Search all databases for best match

This script searches through all available taxonomy databases to find the best
match for a new assembly based on hierarchical GTDB taxonomy.

Usage:
    python assembly_router_multi.py \
        --assembly genome.fna \
        --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \
        --database-dir /path/to/all_databases/ \
        --output-dir results/

The script will:
1. Find all *_db/ directories in --database-dir
2. Query each database to see if taxonomy matches
3. Report best match(es) and recommend action

Author: Generated for HGTool
Date: 2024
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
        return self.levels.get(rank)
    
    def get_most_specific_rank(self) -> Optional[str]:
        for rank in reversed(self.RANKS):
            if self.levels.get(rank):
                return rank
        return None
    
    def get_rank_name(self, rank: str) -> str:
        return self.RANK_NAMES.get(rank, rank)
    
    def matches_at_rank(self, other: 'GTDBTaxonomy', rank: str) -> bool:
        self_val = self.get_rank(rank)
        other_val = other.get_rank(rank)
        
        if not self_val or not other_val:
            return False
        
        return self_val == other_val
    
    def is_within_clade(self, clade: 'GTDBTaxonomy', clade_rank: str) -> bool:
        """Check if this taxonomy belongs to specified clade."""
        # Must match at the clade rank and all higher ranks
        for rank in self.RANKS:
            if rank == clade_rank:
                return self.matches_at_rank(clade, rank)
            
            rank_idx = self.RANKS.index(rank)
            clade_idx = self.RANKS.index(clade_rank)
            
            if rank_idx < clade_idx:
                if not self.matches_at_rank(clade, rank):
                    return False
        
        return True
    
    def __str__(self):
        return self.raw


class DatabaseInfo:
    """Information about a single taxonomy database."""
    
    def __init__(self, db_dir: Path):
        self.db_dir = Path(db_dir)
        self.config_file = self.db_dir / "database_config.json"
        
        if not self.config_file.exists():
            raise FileNotFoundError(f"No config found in {db_dir}")
        
        with open(self.config_file, 'r') as f:
            self.config = json.load(f)
        
        self.clade_name = self.config['clade_name']
        self.clade_taxonomy = GTDBTaxonomy(self.config['clade_taxonomy'])
        self.clade_rank = self.config['clade_rank']
        self.n_genomes = self.config.get('n_genomes', 0)
    
    def query_taxonomy(self, query_taxonomy_str: str) -> Dict:
        """Check if query belongs to this database's clade."""
        query_tax = GTDBTaxonomy(query_taxonomy_str)
        
        is_match = query_tax.is_within_clade(self.clade_taxonomy, self.clade_rank)
        
        if is_match:
            query_specific_rank = query_tax.get_most_specific_rank()
            
            return {
                'matched': True,
                'database': self.db_dir,
                'clade_name': self.clade_name,
                'clade_rank': self.clade_rank,
                'clade_rank_name': self.clade_taxonomy.get_rank_name(self.clade_rank),
                'clade_taxonomy': str(self.clade_taxonomy),
                'query_taxonomy': query_taxonomy_str,
                'query_specific_rank': query_specific_rank,
                'query_rank_name': query_tax.get_rank_name(query_specific_rank) if query_specific_rank else 'unknown',
                'n_genomes': self.n_genomes,
                'match_score': self._calculate_match_score(query_tax)
            }
        else:
            return {
                'matched': False,
                'database': self.db_dir,
                'clade_name': self.clade_name,
                'clade_taxonomy': str(self.clade_taxonomy),
                'query_taxonomy': query_taxonomy_str
            }
    
    def _calculate_match_score(self, query_tax: GTDBTaxonomy) -> int:
        """Calculate how well query matches this database (higher = better)."""
        score = 0
        
        # Award points for each matching rank
        for i, rank in enumerate(GTDBTaxonomy.RANKS):
            if query_tax.matches_at_rank(self.clade_taxonomy, rank):
                score += (i + 1) * 10  # Higher ranks worth more
        
        return score


class MultiDatabaseRouter:
    """Routes assemblies by searching multiple databases."""
    
    def __init__(self, database_dir: Path, output_dir: Path, threads: int = 8):
        self.database_dir = Path(database_dir)
        self.output_dir = Path(output_dir)
        self.threads = threads
        
        # Find all databases
        self.databases = self._discover_databases()
        
        if not self.databases:
            raise ValueError(f"No databases found in {database_dir}")
        
        logger.info(f"Found {len(self.databases)} databases in {database_dir}")
        for db in self.databases:
            logger.info(f"  - {db.clade_name} ({db.clade_taxonomy.get_rank_name(db.clade_rank)})")
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging to file
        log_file = self.output_dir / f"routing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logger.addHandler(file_handler)
    
    def _discover_databases(self) -> List[DatabaseInfo]:
        """Find all valid databases in directory."""
        databases = []
        
        # Look for *_db directories
        for db_dir in self.database_dir.glob("*_db"):
            if db_dir.is_dir():
                try:
                    db_info = DatabaseInfo(db_dir)
                    databases.append(db_info)
                except Exception as e:
                    logger.warning(f"Skipping {db_dir.name}: {e}")
        
        return databases
    
    def route_assembly(
        self,
        assembly_path: Path,
        taxonomy: str,
        assembly_id: Optional[str] = None
    ) -> Dict:
        """Route assembly by querying all databases."""
        if assembly_id is None:
            assembly_id = Path(assembly_path).stem
        
        logger.info(f"Routing assembly: {assembly_id}")
        logger.info(f"Taxonomy: {taxonomy}")
        logger.info(f"Searching {len(self.databases)} databases...")
        
        # Query all databases
        matches = []
        for db in self.databases:
            result = db.query_taxonomy(taxonomy)
            if result['matched']:
                matches.append(result)
        
        # Determine action based on matches
        if len(matches) == 0:
            return self._route_no_match(assembly_path, assembly_id, taxonomy)
        elif len(matches) == 1:
            return self._route_single_match(assembly_path, assembly_id, matches[0])
        else:
            return self._route_multiple_matches(assembly_path, assembly_id, matches)
    
    def _route_no_match(self, assembly_path: Path, assembly_id: str, taxonomy: str) -> Dict:
        """No database matches."""
        logger.info("✗ No matching databases found")
        logger.info("  → Requires new OrthoPhyl run for this clade")
        
        query_tax = GTDBTaxonomy(taxonomy)
        specific_rank = query_tax.get_most_specific_rank()
        
        decision = {
            'pipeline': 'OrthoPhyl',
            'reason': 'No existing database covers this taxonomy',
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': taxonomy,
            'query_rank': specific_rank,
            'query_rank_name': query_tax.get_rank_name(specific_rank) if specific_rank else 'unknown',
            'matches_found': 0,
            'available_clades': [db.clade_name for db in self.databases],
            'suggestion': self._suggest_new_database(taxonomy),
            'command': self._build_orthophyl_command(assembly_path, taxonomy)
        }
        
        self._save_decision(decision)
        return decision
    
    def _route_single_match(self, assembly_path: Path, assembly_id: str, match: Dict) -> Dict:
        """Single database matches."""
        logger.info(f"✓ Matched database: {match['clade_name']}")
        logger.info(f"  Clade: {match['clade_taxonomy']}")
        logger.info(f"  Rank: {match['clade_rank_name']}")
        logger.info(f"  Genomes: {match['n_genomes']}")
        logger.info(f"  → Use ReLeaf to add to existing phylogeny")
        
        decision = {
            'pipeline': 'ReLeaf',
            'reason': f"Taxonomy belongs to {match['clade_name']} database",
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': match['query_taxonomy'],
            'matched_database': str(match['database']),
            'clade_name': match['clade_name'],
            'clade_taxonomy': match['clade_taxonomy'],
            'clade_rank': match['clade_rank'],
            'clade_rank_name': match['clade_rank_name'],
            'database_genomes': match['n_genomes'],
            'matches_found': 1,
            'match_score': match['match_score'],
            'command': self._build_releaf_command(assembly_path, assembly_id, match['database'])
        }
        
        self._save_decision(decision)
        return decision
    
    def _route_multiple_matches(self, assembly_path: Path, assembly_id: str, matches: List[Dict]) -> Dict:
        """Multiple databases match - choose best."""
        logger.info(f"⚠ Multiple matches found ({len(matches)} databases)")
        
        # Sort by match score (higher = more specific match)
        matches_sorted = sorted(matches, key=lambda m: m['match_score'], reverse=True)
        best_match = matches_sorted[0]
        
        logger.info(f"  Best match: {best_match['clade_name']} (score: {best_match['match_score']})")
        for match in matches_sorted[1:]:
            logger.info(f"  Also matches: {match['clade_name']} (score: {match['match_score']})")
        
        logger.info(f"  → Recommending: {best_match['clade_name']}")
        
        decision = {
            'pipeline': 'ReLeaf',
            'reason': f"Best match among {len(matches)} databases: {best_match['clade_name']}",
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': best_match['query_taxonomy'],
            'matched_database': str(best_match['database']),
            'clade_name': best_match['clade_name'],
            'clade_taxonomy': best_match['clade_taxonomy'],
            'clade_rank': best_match['clade_rank'],
            'clade_rank_name': best_match['clade_rank_name'],
            'database_genomes': best_match['n_genomes'],
            'matches_found': len(matches),
            'match_score': best_match['match_score'],
            'all_matches': [
                {
                    'clade_name': m['clade_name'],
                    'clade_rank': m['clade_rank_name'],
                    'score': m['match_score'],
                    'database': str(m['database'])
                }
                for m in matches_sorted
            ],
            'command': self._build_releaf_command(assembly_path, assembly_id, best_match['database'])
        }
        
        self._save_decision(decision)
        return decision
    
    def _suggest_new_database(self, taxonomy: str) -> str:
        """Suggest what database to create."""
        query_tax = GTDBTaxonomy(taxonomy)
        specific_rank = query_tax.get_most_specific_rank()
        
        if not specific_rank:
            return "Run OrthoPhyl with related genomes"
        
        rank_name = query_tax.get_rank_name(specific_rank)
        rank_value = query_tax.get_rank(specific_rank)
        
        return f"Create new database at {rank_name} level: {rank_value}"
    
    def _build_releaf_command(self, assembly_path: Path, assembly_id: str, database_dir: Path) -> str:
        """Build ReLeaf command."""
        releaf_input_dir = self.output_dir / "releaf_input"
        releaf_input_dir.mkdir(parents=True, exist_ok=True)
        
        dest = releaf_input_dir / f"{assembly_id}.fna"
        if not dest.exists():
            shutil.copy(assembly_path, dest)
        
        orthophyl_run = database_dir / "orthophyl_run"
        
        cmd = [
            "./ReLeaf.sh",
            "--store", str(orthophyl_run),
            "--input_genomes", str(releaf_input_dir),
            "-t", str(self.threads),
            "--tree_method", "iqtree",
            "--TREE_DATA", "CDS"
        ]
        
        return ' '.join(cmd)
    
    def _build_orthophyl_command(self, assembly_path: Path, taxonomy: str) -> str:
        """Build OrthoPhyl command."""
        orthophyl_input_dir = self.output_dir / "orthophyl_new_clade"
        orthophyl_input_dir.mkdir(parents=True, exist_ok=True)
        
        dest = orthophyl_input_dir / Path(assembly_path).name
        if not dest.exists():
            shutil.copy(assembly_path, dest)
        
        query_tax = GTDBTaxonomy(taxonomy)
        specific_rank = query_tax.get_most_specific_rank()
        
        download_note = "# Add related genomes to: " + str(orthophyl_input_dir)
        if specific_rank:
            rank_val = query_tax.get_rank(specific_rank)
            rank_name = query_tax.get_rank_name(specific_rank)
            download_note += f"\\n# Suggested: Download genomes from {rank_name} '{rank_val}'"
        
        cmd = [
            "# " + download_note,
            "\\n./OrthoPhyl.sh",
            "-g", str(orthophyl_input_dir),
            "-r", Path(assembly_path).stem,
            "-o", str(self.output_dir / "orthophyl_output"),
            "-t", str(self.threads),
            "--tree_method", "iqtree",
            "--TREE_DATA", "CDS"
        ]
        
        return ' '.join(cmd)
    
    def _save_decision(self, decision: Dict):
        """Save routing decision."""
        decision_file = self.output_dir / f"routing_decision_{decision['assembly_id']}.json"
        
        with open(decision_file, 'w') as f:
            json.dump(decision, f, indent=2)
        
        logger.info(f"Decision saved to {decision_file}")
        
        # Human-readable summary
        summary_file = self.output_dir / f"routing_summary_{decision['assembly_id']}.txt"
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("ASSEMBLY ROUTING DECISION (Multi-Database Search)\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Assembly ID: {decision['assembly_id']}\n")
            f.write(f"Assembly: {decision['assembly']}\n")
            f.write(f"Query Taxonomy: {decision['query_taxonomy']}\n")
            f.write(f"Pipeline: {decision['pipeline']}\n")
            f.write(f"Reason: {decision['reason']}\n\n")
            
            f.write(f"Databases searched: {len(self.databases)}\n")
            f.write(f"Matches found: {decision['matches_found']}\n\n")
            
            if decision['pipeline'] == 'ReLeaf':
                f.write("Matched Database:\n")
                f.write(f"  Name: {decision['clade_name']}\n")
                f.write(f"  Taxonomy: {decision['clade_taxonomy']}\n")
                f.write(f"  Rank: {decision['clade_rank_name']}\n")
                f.write(f"  Genomes: {decision['database_genomes']}\n")
                f.write(f"  Location: {decision['matched_database']}\n")
                
                if decision['matches_found'] > 1:
                    f.write("\nAlternative matches:\n")
                    for alt in decision.get('all_matches', [])[1:]:
                        f.write(f"  - {alt['clade_name']} ({alt['clade_rank']}, score: {alt['score']})\n")
            else:
                f.write("No matching databases found.\n")
                f.write(f"Suggestion: {decision['suggestion']}\n")
                f.write(f"\nAvailable databases:\n")
                for clade in decision.get('available_clades', []):
                    f.write(f"  - {clade}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("COMMAND TO RUN\n")
            f.write("=" * 70 + "\n\n")
            f.write(decision['command'] + "\n")
        
        logger.info(f"Summary saved to {summary_file}")
    
    def batch_route(self, input_table: Path) -> List[Dict]:
        """Route multiple assemblies."""
        logger.info(f"Batch routing from {input_table}")
        
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
                
                assembly_path = Path(fields[0])
                taxonomy = fields[1]
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
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("BATCH ROUTING SUMMARY (Multi-Database)\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Databases searched: {len(self.databases)}\n")
            for db in self.databases:
                f.write(f"  - {db.clade_name}\n")
            f.write("\n")
            
            f.write(f"Total assemblies: {len(decisions)}\n")
            f.write(f"  → ReLeaf (matched): {n_releaf}\n")
            f.write(f"  → OrthoPhyl (no match): {n_orthophyl}\n\n")
            
            if n_releaf > 0:
                f.write("ReLeaf Assemblies:\n")
                f.write("-" * 70 + "\n")
                for d in decisions:
                    if d['pipeline'] == 'ReLeaf':
                        f.write(f"  {d['assembly_id']:30s} → {d['clade_name']}\n")
                f.write("\n")
            
            if n_orthophyl > 0:
                f.write("OrthoPhyl Assemblies (no database match):\n")
                f.write("-" * 70 + "\n")
                for d in decisions:
                    if d['pipeline'] == 'OrthoPhyl':
                        f.write(f"  {d['assembly_id']:30s}\n")
                f.write("\n")
        
        logger.info(f"Batch summary: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Route assemblies by searching all available databases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Route single assembly (searches all databases)
  python assembly_router_multi.py \\
      --assembly genome.fna \\
      --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \\
      --database-dir /path/to/all_databases/ \\
      --output-dir results/

  # Batch route
  python assembly_router_multi.py \\
      --batch assemblies.tsv \\
      --database-dir /path/to/all_databases/ \\
      --output-dir results/
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--assembly', help='Single assembly FASTA')
    input_group.add_argument('--batch', help='TSV: assembly_path, taxonomy, [id]')
    
    parser.add_argument('--taxonomy', help='GTDB taxonomy (for --assembly)')
    parser.add_argument('--assembly-id', help='Assembly identifier')
    parser.add_argument('--database-dir', required=True, help='Directory containing *_db/ subdirectories')
    parser.add_argument('--output-dir', default='routing_output', help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Threads')
    
    args = parser.parse_args()
    
    if args.assembly and not args.taxonomy:
        parser.error("--taxonomy required with --assembly")
    
    try:
        router = MultiDatabaseRouter(
            database_dir=Path(args.database_dir),
            output_dir=Path(args.output_dir),
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
                Path(args.assembly),
                args.taxonomy,
                args.assembly_id
            )
            decisions = [decision]
        
        print("\n" + "=" * 70)
        print("ROUTING COMPLETE")
        print("=" * 70)
        
        for decision in decisions:
            print(f"\n{decision['assembly_id']}:")
            print(f"  Pipeline: {decision['pipeline']}")
            print(f"  Reason: {decision['reason']}")
            if decision['pipeline'] == 'ReLeaf':
                print(f"  Database: {decision['clade_name']}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Routing failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
