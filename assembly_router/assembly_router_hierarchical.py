#!/usr/bin/env python3
"""
Assembly Router - Hierarchical Taxonomy Matching for OrthoPhyl/ReLeaf

Routes assemblies based on hierarchical taxonomy matching (e.g., GTDB format).
Each OrthoPhyl run is associated with a single taxonomic level/clade.

New assemblies are checked to see if they belong to that clade:
- If YES → ReLeaf (add to existing phylogeny)
- If NO  → OrthoPhyl (create new analysis or suggest alternative)

GTDB Taxonomy Format:
    d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__

Usage:
    python assembly_router_hierarchical.py \\
        --assembly assembly.fna \\
        --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \\
        --database-dir reference_db/ \\
        --output-dir results/

Database Structure:
    database_dir/
    ├── database_config.json      (single clade definition per database)
    ├── phylogeny.nwk             (reference tree)
    ├── orthophyl_run/            (OrthoPhyl outputs with HMMs)
    └── genome_list.txt           (genomes in this database)

Author: Generated for HGTool
Date: 2024
"""

import os
import sys
import json
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import logging
from datetime import datetime
import re

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class GTDBTaxonomy:
    """
    Parse and manipulate GTDB-style hierarchical taxonomy strings.
    
    Format: d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
    """
    
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
        """
        Parse GTDB taxonomy string.
        
        Args:
            taxonomy_string: GTDB format taxonomy
        """
        self.raw = taxonomy_string.strip()
        self.levels = self._parse_taxonomy(self.raw)
    
    def _parse_taxonomy(self, taxonomy_string: str) -> Dict[str, str]:
        """Parse taxonomy string into dictionary."""
        levels = {}
        
        # Split by semicolon
        parts = taxonomy_string.split(';')
        
        for part in parts:
            part = part.strip()
            if not part:
                continue
            
            # Extract rank prefix (d__, p__, c__, etc.)
            match = re.match(r'^([dpcofgs])__(.*)$', part)
            if match:
                rank = match.group(1)
                name = match.group(2).strip()
                levels[rank] = name if name else None
            else:
                logger.warning(f"Could not parse taxonomy part: {part}")
        
        return levels
    
    def get_rank(self, rank: str) -> Optional[str]:
        """Get taxonomy at specific rank (e.g., 'o' for order)."""
        return self.levels.get(rank)
    
    def get_rank_name(self, rank: str) -> str:
        """Get full name of rank (e.g., 'o' → 'order')."""
        return self.RANK_NAMES.get(rank, rank)
    
    def is_defined_at_rank(self, rank: str) -> bool:
        """Check if taxonomy is defined (not empty) at this rank."""
        value = self.levels.get(rank)
        return value is not None and value != ''
    
    def matches_at_rank(self, other: 'GTDBTaxonomy', rank: str) -> bool:
        """
        Check if two taxonomies match at specified rank.
        
        Args:
            other: Another GTDBTaxonomy object
            rank: Rank to compare (d, p, c, o, f, g, s)
        
        Returns:
            True if both are defined and match at this rank
        """
        self_val = self.get_rank(rank)
        other_val = other.get_rank(rank)
        
        # Both must be defined
        if not self_val or not other_val:
            return False
        
        return self_val == other_val
    
    def is_within_clade(self, clade: 'GTDBTaxonomy', clade_rank: str) -> bool:
        """
        Check if this taxonomy belongs to specified clade.
        
        Args:
            clade: The clade taxonomy
            clade_rank: The rank at which the clade is defined
        
        Returns:
            True if this taxonomy is within the clade
        """
        # Must match at the clade rank and all higher ranks
        for rank in self.RANKS:
            if rank == clade_rank:
                # Check if we match at the defining rank
                return self.matches_at_rank(clade, rank)
            
            # Check all ranks above the clade rank
            rank_idx = self.RANKS.index(rank)
            clade_idx = self.RANKS.index(clade_rank)
            
            if rank_idx < clade_idx:
                if not self.matches_at_rank(clade, rank):
                    return False
        
        return True
    
    def get_taxonomy_at_rank(self, rank: str) -> str:
        """
        Get taxonomy string up to and including specified rank.
        
        Args:
            rank: Rank to truncate at (d, p, c, o, f, g, s)
        
        Returns:
            Truncated taxonomy string
        """
        parts = []
        for r in self.RANKS:
            value = self.get_rank(r)
            if value:
                parts.append(f"{r}__{value}")
            else:
                parts.append(f"{r}__")
            
            if r == rank:
                break
        
        return ';'.join(parts)
    
    def get_most_specific_rank(self) -> Optional[str]:
        """Get the most specific (lowest) defined rank."""
        for rank in reversed(self.RANKS):
            if self.is_defined_at_rank(rank):
                return rank
        return None
    
    def __str__(self):
        return self.raw
    
    def __repr__(self):
        return f"GTDBTaxonomy('{self.raw}')"


class TaxonomyDatabase:
    """
    Manages a hierarchical taxonomy database for a single OrthoPhyl run.
    """
    
    def __init__(self, database_dir: Path):
        """
        Initialize taxonomy database.
        
        Args:
            database_dir: Path to database directory
        """
        self.database_dir = Path(database_dir)
        self.config_file = self.database_dir / "database_config.json"
        self.phylogeny_file = self.database_dir / "phylogeny.nwk"
        self.genome_list_file = self.database_dir / "genome_list.txt"
        self.orthophyl_dir = self.database_dir / "orthophyl_run"
        
        # Validate and load
        self._validate_database()
        self.config = self._load_config()
        
        # Parse clade taxonomy
        self.clade_taxonomy = GTDBTaxonomy(self.config['clade_taxonomy'])
        self.clade_rank = self.config['clade_rank']
        
        logger.info(f"Loaded database for clade: {self.config['clade_name']}")
        logger.info(f"  Rank: {self.clade_taxonomy.get_rank_name(self.clade_rank)}")
        logger.info(f"  Taxonomy: {self.clade_taxonomy}")
        logger.info(f"  Genomes: {self.config['n_genomes']}")
    
    def _validate_database(self):
        """Validate required database files."""
        if not self.database_dir.exists():
            raise FileNotFoundError(f"Database directory not found: {self.database_dir}")
        
        required_files = [
            self.config_file,
            self.phylogeny_file,
            self.orthophyl_dir
        ]
        
        missing = [f for f in required_files if not Path(f).exists()]
        if missing:
            raise FileNotFoundError(
                f"Missing required database files:\\n" +
                "\\n".join(f"  - {f}" for f in missing)
            )
        
        # Check for HMM profiles
        hmm_dir = self.orthophyl_dir / "OG_alignmentsToHMM" / "hmms_final"
        if not hmm_dir.exists():
            hmm_dir = self.orthophyl_dir / "annots_prots.fixed" / "OrthoFinder" / \
                      "Results_ortho" / "MultipleSequenceAlignments"
        
        if not hmm_dir.exists():
            raise FileNotFoundError(
                f"Could not find HMM profiles in {self.orthophyl_dir}"
            )
    
    def _load_config(self) -> Dict:
        """Load database configuration."""
        with open(self.config_file, 'r') as f:
            return json.load(f)
    
    def query_taxonomy(self, query_taxonomy_str: str) -> Dict:
        """
        Check if query taxonomy belongs to this database's clade.
        
        Args:
            query_taxonomy_str: GTDB taxonomy string to query
        
        Returns:
            Dictionary with match information
        """
        query_tax = GTDBTaxonomy(query_taxonomy_str)
        
        # Check if query is within this clade
        is_match = query_tax.is_within_clade(self.clade_taxonomy, self.clade_rank)
        
        if is_match:
            # Determine match specificity
            query_specific_rank = query_tax.get_most_specific_rank()
            
            return {
                'matched': True,
                'clade_name': self.config['clade_name'],
                'clade_rank': self.clade_rank,
                'clade_rank_name': self.clade_taxonomy.get_rank_name(self.clade_rank),
                'clade_taxonomy': str(self.clade_taxonomy),
                'query_taxonomy': query_taxonomy_str,
                'query_specific_rank': query_specific_rank,
                'query_rank_name': query_tax.get_rank_name(query_specific_rank) if query_specific_rank else 'unknown',
                'confidence': 'high' if query_specific_rank else 'low',
                'n_genomes': self.config['n_genomes']
            }
        else:
            # Find where they diverge
            divergence_rank = self._find_divergence_rank(query_tax)
            
            return {
                'matched': False,
                'clade_name': self.config['clade_name'],
                'clade_rank': self.clade_rank,
                'clade_taxonomy': str(self.clade_taxonomy),
                'query_taxonomy': query_taxonomy_str,
                'divergence_rank': divergence_rank,
                'divergence_rank_name': self.clade_taxonomy.get_rank_name(divergence_rank) if divergence_rank else 'domain',
                'reason': self._explain_mismatch(query_tax, divergence_rank)
            }
    
    def _find_divergence_rank(self, query_tax: GTDBTaxonomy) -> Optional[str]:
        """Find the rank where query diverges from clade."""
        for rank in GTDBTaxonomy.RANKS:
            if not query_tax.matches_at_rank(self.clade_taxonomy, rank):
                return rank
        return None
    
    def _explain_mismatch(self, query_tax: GTDBTaxonomy, divergence_rank: Optional[str]) -> str:
        """Generate human-readable explanation of mismatch."""
        if not divergence_rank:
            return "Query taxonomy is empty or invalid"
        
        rank_name = self.clade_taxonomy.get_rank_name(divergence_rank)
        
        query_val = query_tax.get_rank(divergence_rank) or "undefined"
        clade_val = self.clade_taxonomy.get_rank(divergence_rank) or "undefined"
        
        return f"Diverges at {rank_name}: query has '{query_val}' but clade is '{clade_val}'"


class AssemblyRouter:
    """
    Routes assemblies based on hierarchical taxonomy matching.
    """
    
    def __init__(
        self,
        database_dir: Path,
        output_dir: Path,
        threads: int = 8,
        strict_matching: bool = False
    ):
        """
        Initialize router.
        
        Args:
            database_dir: Path to taxonomy database
            output_dir: Output directory
            threads: Number of threads
            strict_matching: If True, require exact rank matching
        """
        self.database = TaxonomyDatabase(database_dir)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.strict_matching = strict_matching
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging to file
        log_file = self.output_dir / f"routing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logger.addHandler(file_handler)
    
    def route_assembly(
        self,
        assembly_path: Path,
        taxonomy: str,
        assembly_id: Optional[str] = None
    ) -> Dict:
        """
        Route assembly based on hierarchical taxonomy.
        
        Args:
            assembly_path: Path to assembly
            taxonomy: GTDB taxonomy string
            assembly_id: Optional identifier
        
        Returns:
            Routing decision dictionary
        """
        if assembly_id is None:
            assembly_id = Path(assembly_path).stem
        
        logger.info(f"Routing assembly: {assembly_id}")
        logger.info(f"Taxonomy: {taxonomy}")
        
        # Query database
        match_result = self.database.query_taxonomy(taxonomy)
        
        # Route based on match
        if match_result['matched']:
            return self._route_to_releaf(assembly_path, assembly_id, match_result)
        else:
            return self._route_to_orthophyl(assembly_path, assembly_id, taxonomy, match_result)
    
    def _route_to_releaf(
        self,
        assembly_path: Path,
        assembly_id: str,
        match_result: Dict
    ) -> Dict:
        """Route to ReLeaf."""
        logger.info(f"✓ Taxonomy MATCHES database clade")
        logger.info(f"  Clade: {match_result['clade_name']} ({match_result['clade_rank_name']})")
        logger.info(f"  Query defined at: {match_result['query_rank_name']} level")
        logger.info(f"  Confidence: {match_result['confidence']}")
        
        decision = {
            'pipeline': 'ReLeaf',
            'reason': f"Taxonomy belongs to {match_result['clade_name']} clade",
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': match_result['query_taxonomy'],
            'clade_name': match_result['clade_name'],
            'clade_rank': match_result['clade_rank'],
            'clade_taxonomy': match_result['clade_taxonomy'],
            'match_confidence': match_result['confidence'],
            'query_specific_rank': match_result['query_specific_rank'],
            'database_genomes': match_result['n_genomes'],
            'command': self._build_releaf_command(assembly_path, assembly_id)
        }
        
        self._save_decision(decision)
        return decision
    
    def _route_to_orthophyl(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str,
        match_result: Dict
    ) -> Dict:
        """Route to OrthoPhyl."""
        logger.info(f"✗ Taxonomy DOES NOT MATCH database clade")
        logger.info(f"  Database clade: {match_result['clade_name']} ({match_result['clade_taxonomy']})")
        logger.info(f"  Mismatch reason: {match_result['reason']}")
        logger.info(f"  Divergence at: {match_result.get('divergence_rank_name', 'unknown')}")
        logger.info(f"  → Requires new OrthoPhyl run or different database")
        
        decision = {
            'pipeline': 'OrthoPhyl',
            'reason': 'Taxonomy outside of database clade',
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': taxonomy,
            'clade_name': match_result['clade_name'],
            'clade_taxonomy': match_result['clade_taxonomy'],
            'mismatch_reason': match_result['reason'],
            'divergence_rank': match_result.get('divergence_rank'),
            'divergence_rank_name': match_result.get('divergence_rank_name'),
            'suggestion': self._suggest_action(taxonomy, match_result),
            'command': self._build_orthophyl_command(assembly_path, taxonomy)
        }
        
        self._save_decision(decision)
        return decision
    
    def _suggest_action(self, query_taxonomy: str, match_result: Dict) -> str:
        """Suggest what user should do for mismatched taxonomy."""
        suggestions = []
        
        suggestions.append(
            f"Create new OrthoPhyl database for the {match_result.get('divergence_rank_name', 'taxonomic')} "
            f"level containing this assembly"
        )
        
        query_tax = GTDBTaxonomy(query_taxonomy)
        specific_rank = query_tax.get_most_specific_rank()
        
        if specific_rank:
            rank_name = query_tax.get_rank_name(specific_rank)
            rank_value = query_tax.get_rank(specific_rank)
            suggestions.append(
                f"Download related genomes from {rank_name} '{rank_value}' and run OrthoPhyl"
            )
        
        suggestions.append(
            "Check if another database exists for this clade"
        )
        
        return " | ".join(suggestions)
    
    def _build_releaf_command(self, assembly_path: Path, assembly_id: str) -> str:
        """Build ReLeaf command."""
        releaf_input_dir = self.output_dir / "releaf_input"
        releaf_input_dir.mkdir(parents=True, exist_ok=True)
        
        dest = releaf_input_dir / f"{assembly_id}.fna"
        if not dest.exists():
            import shutil
            shutil.copy(assembly_path, dest)
        
        cmd = [
            "./ReLeaf.sh",
            "--store", str(self.database.orthophyl_dir),
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
            import shutil
            shutil.copy(assembly_path, dest)
        
        # Extract taxonomic group for download
        query_tax = GTDBTaxonomy(taxonomy)
        specific_rank = query_tax.get_most_specific_rank()
        
        download_note = "# Add genomes manually to: " + str(orthophyl_input_dir)
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
            f.write("=" * 70 + "\\n")
            f.write("ASSEMBLY ROUTING DECISION (Hierarchical Taxonomy)\\n")
            f.write("=" * 70 + "\\n\\n")
            
            f.write(f"Assembly ID: {decision['assembly_id']}\\n")
            f.write(f"Assembly: {decision['assembly']}\\n")
            f.write(f"Query Taxonomy: {decision['query_taxonomy']}\\n")
            f.write(f"Pipeline: {decision['pipeline']}\\n")
            f.write(f"Reason: {decision['reason']}\\n\\n")
            
            f.write("Database Information:\\n")
            f.write(f"  Clade: {decision['clade_name']}\\n")
            f.write(f"  Clade Taxonomy: {decision['clade_taxonomy']}\\n")
            
            if decision['pipeline'] == 'ReLeaf':
                f.write(f"\\nMatch Details:\\n")
                f.write(f"  Confidence: {decision['match_confidence']}\\n")
                f.write(f"  Query defined at: {decision['query_specific_rank']} level\\n")
                f.write(f"  Database contains: {decision['database_genomes']} genomes\\n")
            else:
                f.write(f"\\nMismatch Details:\\n")
                f.write(f"  Reason: {decision['mismatch_reason']}\\n")
                f.write(f"  Divergence at: {decision.get('divergence_rank_name', 'unknown')}\\n")
                f.write(f"\\nSuggested Action:\\n")
                f.write(f"  {decision['suggestion']}\\n")
            
            f.write("\\n" + "=" * 70 + "\\n")
            f.write("COMMAND TO RUN\\n")
            f.write("=" * 70 + "\\n\\n")
            f.write(decision['command'] + "\\n")
        
        logger.info(f"Summary saved to {summary_file}")
    
    def batch_route(self, input_table: Path) -> List[Dict]:
        """Route multiple assemblies from table."""
        logger.info(f"Batch routing from {input_table}")
        
        decisions = []
        
        with open(input_table, 'r') as f:
            # Skip header if present
            first_line = f.readline().strip()
            if not first_line.startswith('#') and not first_line.lower().startswith('assembly'):
                f.seek(0)
            
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\\t')
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
            f.write("=" * 70 + "\\n")
            f.write("BATCH ROUTING SUMMARY (Hierarchical Taxonomy)\\n")
            f.write("=" * 70 + "\\n\\n")
            
            f.write(f"Database Clade: {self.database.config['clade_name']}\\n")
            f.write(f"Clade Taxonomy: {self.database.clade_taxonomy}\\n\\n")
            
            f.write(f"Total assemblies: {len(decisions)}\\n")
            f.write(f"  → ReLeaf (within clade): {n_releaf}\\n")
            f.write(f"  → OrthoPhyl (outside clade): {n_orthophyl}\\n\\n")
            
            if n_releaf > 0:
                f.write("ReLeaf Assemblies (within clade):\\n")
                f.write("-" * 70 + "\\n")
                for d in decisions:
                    if d['pipeline'] == 'ReLeaf':
                        f.write(f"  {d['assembly_id']:30s} [{d['match_confidence']}]\\n")
                f.write("\\n")
            
            if n_orthophyl > 0:
                f.write("OrthoPhyl Assemblies (outside clade):\\n")
                f.write("-" * 70 + "\\n")
                for d in decisions:
                    if d['pipeline'] == 'OrthoPhyl':
                        f.write(f"  {d['assembly_id']:30s} - {d['mismatch_reason']}\\n")
                f.write("\\n")
        
        logger.info(f"Batch summary saved to {summary_file}")
        print(f"\\nRouting complete: {n_releaf} → ReLeaf, {n_orthophyl} → OrthoPhyl")


def main():
    parser = argparse.ArgumentParser(
        description="Route assemblies using hierarchical taxonomy matching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Route single assembly
  python assembly_router_hierarchical.py \\
      --assembly genome.fna \\
      --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \\
      --database-dir gaiellales_db/ \\
      --output-dir results/

  # Batch route
  python assembly_router_hierarchical.py \\
      --batch assemblies.tsv \\
      --database-dir gaiellales_db/ \\
      --output-dir results/
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--assembly', help='Single assembly FASTA')
    input_group.add_argument('--batch', help='TSV: assembly_path, taxonomy, [id]')
    
    parser.add_argument('--taxonomy', help='GTDB taxonomy string (for --assembly)')
    parser.add_argument('--assembly-id', help='Assembly identifier')
    parser.add_argument('--database-dir', required=True, help='Database directory')
    parser.add_argument('--output-dir', default='routing_output', help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Threads')
    parser.add_argument('--strict', action='store_true', help='Strict rank matching')
    
    args = parser.parse_args()
    
    if args.assembly and not args.taxonomy:
        parser.error("--taxonomy required with --assembly")
    
    try:
        router = AssemblyRouter(
            database_dir=Path(args.database_dir),
            output_dir=Path(args.output_dir),
            threads=args.threads,
            strict_matching=args.strict
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
        
        print("\\n" + "=" * 70)
        print("ROUTING COMPLETE")
        print("=" * 70)
        
        for decision in decisions:
            print(f"\\n{decision['assembly_id']}:")
            print(f"  Pipeline: {decision['pipeline']}")
            print(f"  Reason: {decision['reason']}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Routing failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
"