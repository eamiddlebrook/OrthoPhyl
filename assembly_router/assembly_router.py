#!/usr/bin/env python3
"""
Assembly Router - Route assemblies to ReLeaf or OrthoPhyl based on taxonomy

This script determines whether a new assembly should be:
1. Added to existing phylogeny using ReLeaf (if taxonomy is represented)
2. Trigger new OrthoPhyl run with expanded genome set (if taxonomy is novel)

Usage:
    python assembly_router.py \
        --assembly assembly.fna \
        --taxonomy "Escherichia coli" \
        --database-dir /path/to/reference_db/ \
        --output-dir results/ \
        [options]

Database Directory Structure (expected):
    database_dir/
    ├── taxonomy_map.tsv          (taxon → representative genomes)
    ├── phylogeny.nwk             (reference tree)
    ├── orthophyl_run/            (OrthoPhyl outputs)
    │   ├── hmms_final/
    │   ├── AlignmentsProts.trm/
    │   └── ...
    └── metadata.json             (database info)

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

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class TaxonomyDatabase:
    """
    Handles loading and querying the reference taxonomy database.
    """
    
    def __init__(self, database_dir: Path):
        """
        Initialize taxonomy database.
        
        Args:
            database_dir: Path to reference database directory
        """
        self.database_dir = Path(database_dir)
        self.taxonomy_map_file = self.database_dir / "taxonomy_map.tsv"
        self.phylogeny_file = self.database_dir / "phylogeny.nwk"
        self.metadata_file = self.database_dir / "metadata.json"
        self.orthophyl_dir = self.database_dir / "orthophyl_run"
        
        # Validate database
        self._validate_database()
        
        # Load taxonomy map
        self.taxonomy_map = self._load_taxonomy_map()
        
        # Load metadata
        self.metadata = self._load_metadata()
        
        logger.info(f"Loaded taxonomy database from {database_dir}")
        logger.info(f"Database contains {len(self.taxonomy_map)} taxonomic groups")
    
    def _validate_database(self):
        """Validate that required database files exist."""
        if not self.database_dir.exists():
            raise FileNotFoundError(f"Database directory not found: {self.database_dir}")
        
        required_files = [
            self.taxonomy_map_file,
            self.phylogeny_file,
            self.orthophyl_dir
        ]
        
        missing = [f for f in required_files if not Path(f).exists()]
        if missing:
            raise FileNotFoundError(
                f"Missing required database files:\n" + 
                "\n".join(f"  - {f}" for f in missing)
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
    
    def _load_taxonomy_map(self) -> Dict[str, Dict]:
        """
        Load taxonomy mapping file.
        
        Returns:
            Dictionary mapping taxonomy strings to genome info
        """
        taxonomy_map = {}
        
        with open(self.taxonomy_map_file, 'r') as f:
            header = f.readline().strip().split('\t')
            
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                taxonomy = fields[0]
                level = fields[1]  # genus, species, strain
                genomes = fields[2].split(',')
                
                taxonomy_map[taxonomy] = {
                    'level': level,
                    'genomes': genomes,
                    'count': len(genomes)
                }
        
        return taxonomy_map
    
    def _load_metadata(self) -> Dict:
        """Load database metadata."""
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        else:
            logger.warning("No metadata.json found, using defaults")
            return {
                'created': 'unknown',
                'version': '1.0',
                'n_genomes': len(self.get_all_genomes()),
                'taxonomic_levels': ['genus', 'species', 'strain']
            }
    
    def get_all_genomes(self) -> Set[str]:
        """Get set of all genomes in database."""
        all_genomes = set()
        for tax_info in self.taxonomy_map.values():
            all_genomes.update(tax_info['genomes'])
        return all_genomes
    
    def query_taxonomy(self, taxonomy: str, match_level: str = 'auto') -> Optional[Dict]:
        """
        Query database for taxonomy.
        
        Args:
            taxonomy: Taxonomy string (e.g., "Escherichia coli")
            match_level: Match at 'genus', 'species', 'strain', or 'auto'
        
        Returns:
            Dictionary with match info, or None if no match
        """
        # Direct match
        if taxonomy in self.taxonomy_map:
            return {
                'matched': True,
                'match_type': 'exact',
                'taxonomy': taxonomy,
                'info': self.taxonomy_map[taxonomy]
            }
        
        # Try hierarchical matching if auto
        if match_level == 'auto':
            parts = taxonomy.split()
            
            # Try species-level (first two words)
            if len(parts) >= 2:
                species = ' '.join(parts[:2])
                if species in self.taxonomy_map:
                    return {
                        'matched': True,
                        'match_type': 'species',
                        'taxonomy': species,
                        'query': taxonomy,
                        'info': self.taxonomy_map[species]
                    }
            
            # Try genus-level (first word)
            if len(parts) >= 1:
                genus = parts[0]
                if genus in self.taxonomy_map:
                    return {
                        'matched': True,
                        'match_type': 'genus',
                        'taxonomy': genus,
                        'query': taxonomy,
                        'info': self.taxonomy_map[genus]
                    }
        
        # No match
        return {
            'matched': False,
            'query': taxonomy,
            'suggestions': self._find_similar_taxa(taxonomy)
        }
    
    def _find_similar_taxa(self, taxonomy: str, max_results: int = 5) -> List[str]:
        """Find taxonomically similar entries."""
        parts = taxonomy.lower().split()
        
        candidates = []
        for tax in self.taxonomy_map.keys():
            tax_parts = tax.lower().split()
            
            # Calculate simple word overlap
            overlap = len(set(parts) & set(tax_parts))
            if overlap > 0:
                candidates.append((tax, overlap))
        
        # Sort by overlap and return top matches
        candidates.sort(key=lambda x: x[1], reverse=True)
        return [tax for tax, score in candidates[:max_results]]


class AssemblyRouter:
    """
    Routes assemblies to ReLeaf or OrthoPhyl based on taxonomy coverage.
    """
    
    def __init__(
        self,
        database_dir: Path,
        output_dir: Path,
        ncbi_datasets_path: Optional[str] = None,
        threads: int = 8,
        min_genomes_for_orthophyl: int = 5
    ):
        """
        Initialize router.
        
        Args:
            database_dir: Path to reference taxonomy database
            output_dir: Output directory for results
            ncbi_datasets_path: Path to ncbi-datasets CLI tool
            threads: Number of threads for analysis
            min_genomes_for_orthophyl: Min additional genomes to trigger OrthoPhyl
        """
        self.database = TaxonomyDatabase(database_dir)
        self.output_dir = Path(output_dir)
        self.ncbi_datasets = ncbi_datasets_path or "datasets"
        self.threads = threads
        self.min_genomes_for_orthophyl = min_genomes_for_orthophyl
        
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
        Route assembly to appropriate pipeline.
        
        Args:
            assembly_path: Path to assembly FASTA
            taxonomy: Taxonomic classification
            assembly_id: Optional assembly identifier
        
        Returns:
            Dictionary with routing decision and parameters
        """
        if assembly_id is None:
            assembly_id = Path(assembly_path).stem
        
        logger.info(f"Routing assembly: {assembly_id}")
        logger.info(f"Taxonomy: {taxonomy}")
        
        # Query database
        match_result = self.database.query_taxonomy(taxonomy)
        
        # Decision logic
        if match_result['matched']:
            # Taxonomy is represented - use ReLeaf
            return self._route_to_releaf(assembly_path, assembly_id, taxonomy, match_result)
        else:
            # Novel taxonomy - need OrthoPhyl with expanded set
            return self._route_to_orthophyl(assembly_path, assembly_id, taxonomy, match_result)
    
    def _route_to_releaf(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str,
        match_result: Dict
    ) -> Dict:
        """Route to ReLeaf pipeline."""
        logger.info(f"✓ Taxonomy matched in database ({match_result['match_type']} level)")
        logger.info(f"  Matched: {match_result['taxonomy']}")
        logger.info(f"  Representative genomes: {match_result['info']['count']}")
        
        decision = {
            'pipeline': 'ReLeaf',
            'reason': f"Taxonomy represented in database at {match_result['match_type']} level",
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'taxonomy': taxonomy,
            'matched_taxonomy': match_result['taxonomy'],
            'match_type': match_result['match_type'],
            'database_genomes': match_result['info']['genomes'],
            'command': self._build_releaf_command(assembly_path, assembly_id)
        }
        
        # Save decision
        self._save_decision(decision)
        
        return decision
    
    def _route_to_orthophyl(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str,
        match_result: Dict
    ) -> Dict:
        """Route to OrthoPhyl with expanded genome set."""
        logger.info(f"✗ Taxonomy NOT found in database")
        logger.info(f"  Triggering OrthoPhyl with expanded genome set")
        
        if match_result.get('suggestions'):
            logger.info(f"  Similar taxa in database: {', '.join(match_result['suggestions'][:3])}")
        
        # Download related genomes
        related_genomes = self._find_and_download_related_genomes(taxonomy)
        
        decision = {
            'pipeline': 'OrthoPhyl',
            'reason': 'Novel taxonomy not represented in database',
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'taxonomy': taxonomy,
            'suggestions': match_result.get('suggestions', []),
            'related_genomes_downloaded': len(related_genomes),
            'related_genomes': related_genomes,
            'command': self._build_orthophyl_command(assembly_path, related_genomes)
        }
        
        # Save decision
        self._save_decision(decision)
        
        return decision
    
    def _find_and_download_related_genomes(
        self,
        taxonomy: str,
        max_genomes: int = 50
    ) -> List[str]:
        """
        Find and download related genomes from NCBI.
        
        Args:
            taxonomy: Target taxonomy
            max_genomes: Maximum genomes to download
        
        Returns:
            List of downloaded genome paths
        """
        logger.info(f"Searching NCBI for related genomes: {taxonomy}")
        
        # Parse taxonomy
        parts = taxonomy.split()
        if len(parts) >= 2:
            search_term = ' '.join(parts[:2])  # Use genus + species
            level = "species"
        else:
            search_term = parts[0]  # Use genus only
            level = "genus"
        
        logger.info(f"Searching at {level} level: {search_term}")
        
        # Create download directory
        download_dir = self.output_dir / "downloaded_genomes" / taxonomy.replace(' ', '_')
        download_dir.mkdir(parents=True, exist_ok=True)
        
        # Try to use ncbi-datasets CLI
        try:
            downloaded = self._download_via_ncbi_datasets(
                search_term, 
                download_dir, 
                max_genomes
            )
            
            if downloaded:
                logger.info(f"Downloaded {len(downloaded)} genomes from NCBI")
                return downloaded
            else:
                logger.warning("No genomes downloaded from NCBI")
                
        except Exception as e:
            logger.error(f"Failed to download from NCBI: {e}")
        
        # Fallback: create placeholder
        logger.warning("NCBI download failed, will need manual genome addition")
        placeholder = download_dir / "README.txt"
        placeholder.write_text(
            f"Please manually add genomes for '{taxonomy}' to this directory.\n"
            f"Search term used: {search_term}\n"
            f"Download from: https://www.ncbi.nlm.nih.gov/datasets/genome/\n"
        )
        return []
    
    def _download_via_ncbi_datasets(
        self,
        search_term: str,
        output_dir: Path,
        max_genomes: int
    ) -> List[str]:
        """
        Download genomes using ncbi-datasets CLI.
        
        Args:
            search_term: Taxonomy search term
            output_dir: Where to save genomes
            max_genomes: Maximum number to download
        
        Returns:
            List of downloaded genome paths
        """
        # Check if ncbi-datasets is available
        try:
            subprocess.run(
                [self.ncbi_datasets, "--version"],
                check=True,
                capture_output=True
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"ncbi-datasets not found at {self.ncbi_datasets}")
            logger.info("Install with: conda install -c conda-forge ncbi-datasets-cli")
            return []
        
        # Download genomes
        cmd = [
            self.ncbi_datasets,
            "download",
            "genome",
            "taxon",
            search_term,
            "--assembly-level", "complete,chromosome,scaffold",
            "--exclude-atypical",
            "--limit", str(max_genomes),
            "--filename", str(output_dir / "ncbi_dataset.zip")
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            
            # Unzip
            import zipfile
            with zipfile.ZipFile(output_dir / "ncbi_dataset.zip", 'r') as zip_ref:
                zip_ref.extractall(output_dir)
            
            # Find downloaded FASTAs
            genome_files = list((output_dir / "ncbi_dataset" / "data").rglob("*.fna"))
            
            return [str(f) for f in genome_files]
            
        except subprocess.CalledProcessError as e:
            logger.error(f"NCBI download failed: {e.stderr}")
            return []
    
    def _build_releaf_command(self, assembly_path: Path, assembly_id: str) -> str:
        """Build ReLeaf command string."""
        # Copy assembly to ReLeaf input directory
        releaf_input_dir = self.output_dir / "releaf_input"
        releaf_input_dir.mkdir(parents=True, exist_ok=True)
        
        # Symlink or copy assembly
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
    
    def _build_orthophyl_command(
        self,
        assembly_path: Path,
        related_genomes: List[str]
    ) -> str:
        """Build OrthoPhyl command string."""
        # Create genome directory with new assembly + related genomes
        orthophyl_input_dir = self.output_dir / "orthophyl_input"
        orthophyl_input_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy/symlink all genomes
        genome_paths = [assembly_path] + related_genomes
        for genome in genome_paths:
            genome_path = Path(genome)
            dest = orthophyl_input_dir / genome_path.name
            if not dest.exists():
                import shutil
                shutil.copy(genome_path, dest)
        
        cmd = [
            "./OrthoPhyl.sh",
            "-g", str(orthophyl_input_dir),
            "-r", Path(assembly_path).stem,
            "-o", str(self.output_dir / "orthophyl_output"),
            "-t", str(self.threads),
            "--tree_method", "iqtree",
            "--TREE_DATA", "CDS",
            "--use_partitions", "true"
        ]
        
        return ' '.join(cmd)
    
    def _save_decision(self, decision: Dict):
        """Save routing decision to JSON file."""
        decision_file = self.output_dir / f"routing_decision_{decision['assembly_id']}.json"
        
        with open(decision_file, 'w') as f:
            json.dump(decision, f, indent=2)
        
        logger.info(f"Decision saved to {decision_file}")
        
        # Also create human-readable summary
        summary_file = self.output_dir / f"routing_summary_{decision['assembly_id']}.txt"
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("ASSEMBLY ROUTING DECISION\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Assembly ID: {decision['assembly_id']}\n")
            f.write(f"Assembly: {decision['assembly']}\n")
            f.write(f"Taxonomy: {decision['taxonomy']}\n")
            f.write(f"Pipeline: {decision['pipeline']}\n")
            f.write(f"Reason: {decision['reason']}\n\n")
            
            if decision['pipeline'] == 'ReLeaf':
                f.write(f"Matched taxonomy: {decision['matched_taxonomy']}\n")
                f.write(f"Match type: {decision['match_type']}\n")
                f.write(f"Database genomes: {len(decision['database_genomes'])}\n")
            else:
                f.write(f"Related genomes downloaded: {decision['related_genomes_downloaded']}\n")
                if decision.get('suggestions'):
                    f.write(f"Similar taxa in database: {', '.join(decision['suggestions'][:3])}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("COMMAND TO RUN\n")
            f.write("=" * 70 + "\n\n")
            f.write(decision['command'] + "\n")
        
        logger.info(f"Summary saved to {summary_file}")
    
    def batch_route(self, input_table: Path) -> List[Dict]:
        """
        Route multiple assemblies from a table.
        
        Args:
            input_table: TSV file with columns: assembly_path, taxonomy, [assembly_id]
        
        Returns:
            List of routing decisions
        """
        logger.info(f"Batch routing from {input_table}")
        
        decisions = []
        
        with open(input_table, 'r') as f:
            # Skip header if present
            first_line = f.readline().strip()
            if not first_line.startswith('#') and not first_line.lower().startswith('assembly'):
                # No header, process first line
                f.seek(0)
            
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    logger.warning(f"Line {line_num}: insufficient fields, skipping")
                    continue
                
                assembly_path = Path(fields[0])
                taxonomy = fields[1]
                assembly_id = fields[2] if len(fields) > 2 else None
                
                try:
                    decision = self.route_assembly(assembly_path, taxonomy, assembly_id)
                    decisions.append(decision)
                except Exception as e:
                    logger.error(f"Failed to route {assembly_path}: {e}")
        
        # Generate batch summary
        self._generate_batch_summary(decisions)
        
        return decisions
    
    def _generate_batch_summary(self, decisions: List[Dict]):
        """Generate summary of batch routing."""
        summary_file = self.output_dir / "batch_routing_summary.txt"
        
        n_releaf = sum(1 for d in decisions if d['pipeline'] == 'ReLeaf')
        n_orthophyl = sum(1 for d in decisions if d['pipeline'] == 'OrthoPhyl')
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("BATCH ROUTING SUMMARY\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Total assemblies processed: {len(decisions)}\n")
            f.write(f"Routed to ReLeaf: {n_releaf}\n")
            f.write(f"Routed to OrthoPhyl: {n_orthophyl}\n\n")
            
            if n_releaf > 0:
                f.write("ReLeaf Assemblies (taxonomy represented):\n")
                f.write("-" * 70 + "\n")
                for d in decisions:
                    if d['pipeline'] == 'ReLeaf':
                        f.write(f"  {d['assembly_id']:30s} {d['taxonomy']:40s} → {d['matched_taxonomy']}\n")
                f.write("\n")
            
            if n_orthophyl > 0:
                f.write("OrthoPhyl Assemblies (novel taxonomy):\n")
                f.write("-" * 70 + "\n")
                for d in decisions:
                    if d['pipeline'] == 'OrthoPhyl':
                        f.write(f"  {d['assembly_id']:30s} {d['taxonomy']:40s} ({d['related_genomes_downloaded']} genomes)\n")
                f.write("\n")
        
        logger.info(f"Batch summary saved to {summary_file}")
        print(f"\nBatch summary: {n_releaf} → ReLeaf, {n_orthophyl} → OrthoPhyl")


def create_example_database(output_dir: Path):
    """
    Create an example taxonomy database structure.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create taxonomy_map.tsv
    taxonomy_map = output_dir / "taxonomy_map.tsv"
    with open(taxonomy_map, 'w') as f:
        f.write("taxonomy\tlevel\tgenomes\n")
        f.write("Escherichia\tgenus\tGCF_000005845.2,GCF_000008865.2\n")
        f.write("Escherichia coli\tspecies\tGCF_000005845.2,GCF_000008865.2,GCF_000482265.1\n")
        f.write("Salmonella\tgenus\tGCF_000006945.2,GCF_000011885.1\n")
        f.write("Salmonella enterica\tspecies\tGCF_000006945.2,GCF_000011885.1,GCF_000022165.1\n")
    
    # Create metadata.json
    metadata = {
        "created": datetime.now().isoformat(),
        "version": "1.0",
        "description": "Example taxonomy database",
        "n_genomes": 7,
        "n_taxa": 4
    }
    
    with open(output_dir / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    # Create placeholder files
    (output_dir / "phylogeny.nwk").write_text("(Escherichia:0.1,Salmonella:0.1);")
    (output_dir / "orthophyl_run").mkdir(exist_ok=True)
    (output_dir / "orthophyl_run" / "OG_alignmentsToHMM" / "hmms_final").mkdir(parents=True, exist_ok=True)
    
    print(f"Example database created at {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Route assemblies to ReLeaf or OrthoPhyl based on taxonomy",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument('--assembly', help='Single assembly FASTA file')
    input_group.add_argument('--batch', help='TSV file with multiple assemblies')
    input_group.add_argument('--create-example-db', metavar='DIR', help='Create example database')
    
    parser.add_argument('--taxonomy', help='Taxonomic classification (for single assembly)')
    parser.add_argument('--assembly-id', help='Assembly identifier')
    parser.add_argument('--database-dir', help='Path to reference database directory')
    parser.add_argument('--output-dir', default='routing_output', help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads')
    parser.add_argument('--execute', action='store_true', help='Execute the pipeline')
    parser.add_argument('--dry-run', action='store_true', help='Show decision only')
    
    args = parser.parse_args()
    
    # Handle example database creation
    if args.create_example_db:
        create_example_database(Path(args.create_example_db))
        return 0
    
    # Validate inputs
    if not args.batch and not args.assembly:
        parser.error("Either --assembly or --batch must be specified")
    
    if args.assembly and not args.taxonomy:
        parser.error("--taxonomy is required with --assembly")
    
    if not args.database_dir:
        parser.error("--database-dir is required")
    
    # Initialize router
    try:
        router = AssemblyRouter(
            database_dir=Path(args.database_dir),
            output_dir=Path(args.output_dir),
            threads=args.threads
        )
    except Exception as e:
        logger.error(f"Failed to initialize router: {e}")
        return 1
    
    # Route assemblies
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
        
        # Print summary
        print("\n" + "=" * 70)
        print("ROUTING COMPLETE")
        print("=" * 70)
        
        for decision in decisions:
            print(f"\n{decision['assembly_id']}:")
            print(f"  Pipeline: {decision['pipeline']}")
            print(f"  Reason: {decision['reason']}")
            print(f"  Command: {decision['command']}")
        
        if args.dry_run:
            print("\n[DRY RUN] Commands not executed")
        elif not args.execute:
            print("\n[INFO] Use --execute to run the pipeline")
        
        return 0
        
    except Exception as e:
        logger.error(f"Routing failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())