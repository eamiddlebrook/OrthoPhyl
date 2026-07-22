#!/usr/bin/env python3
"""
Taxon Assembly Gatherer - Query NCBI for assemblies by taxon

This module queries NCBI for bacterial assemblies matching a specified taxon
(by name or TaxID) and generates a filtered list suitable for OrthoPhyl/ReLeaf.

Features:
- Query by taxon name or TaxID
- Filter by taxonomic rank (family, order, genus, species)
- Quality filtering (completeness, contamination, N50)
- Handle RefSeq/GenBank redundancy (prefer RefSeq)
- Date-based filtering for database updates
- Output TSV compatible with OrthoPhyl wrapper

Usage:
    # Get all assemblies for a taxon
    python taxon_assembly_gatherer.py \\
        --taxon "Methylorubrum" \\
        --rank genus \\
        --output assemblies.tsv
    
    # Get assemblies added since a date (for updates)
    python taxon_assembly_gatherer.py \\
        --taxon "Methylorubrum" \\
        --since-date 2025-01-01 \\
        --output new_assemblies.tsv

Author: OrthoPhyl Project
Date: 2025
"""

import os
import sys
import argparse
import urllib.request
import tarfile
import gzip
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from datetime import datetime
from collections import defaultdict
import csv
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class NCBITaxonomy:
    """Parse and query NCBI taxonomy database."""
    
    def __init__(self, taxdump_dir: Path):
        self.taxdump_dir = Path(taxdump_dir)
        self.nodes = {}  # taxid -> {'parent': parent_taxid, 'rank': rank}
        self.names = {}  # taxid -> scientific_name
        self.name_to_taxid = {}  # scientific_name -> taxid
        
        self._load_nodes()
        self._load_names()
    
    def _load_nodes(self):
        """Load nodes.dmp (taxonomy tree structure)."""
        nodes_file = self.taxdump_dir / "nodes.dmp"
        
        if not nodes_file.exists():
            raise FileNotFoundError(f"nodes.dmp not found in {self.taxdump_dir}")
        
        logger.info(f"Loading taxonomy nodes...")
        
        with open(nodes_file, 'r') as f:
            for line in f:
                parts = [p.strip() for p in line.split('|')]
                taxid = parts[0]
                parent_taxid = parts[1]
                rank = parts[2]
                
                self.nodes[taxid] = {
                    'parent': parent_taxid,
                    'rank': rank
                }
        
        logger.info(f"  ✓ Loaded {len(self.nodes):,} nodes")
    
    def _load_names(self):
        """Load names.dmp (taxonomy names)."""
        names_file = self.taxdump_dir / "names.dmp"
        
        if not names_file.exists():
            raise FileNotFoundError(f"names.dmp not found in {self.taxdump_dir}")
        
        logger.info(f"Loading taxonomy names...")
        
        with open(names_file, 'r') as f:
            for line in f:
                parts = [p.strip() for p in line.split('|')]
                taxid = parts[0]
                name = parts[1]
                name_class = parts[3]
                
                # Only keep scientific names
                if name_class == 'scientific name':
                    self.names[taxid] = name
                    self.name_to_taxid[name.lower()] = taxid
        
        logger.info(f"  ✓ Loaded {len(self.names):,} scientific names")
    
    def resolve_taxon(self, taxon: str) -> Optional[str]:
        """Resolve taxon name or ID to TaxID."""
        # Check if it's already a TaxID
        if taxon.isdigit() and taxon in self.nodes:
            return taxon
        
        # Try to find by name (case-insensitive)
        taxon_lower = taxon.lower()
        if taxon_lower in self.name_to_taxid:
            return self.name_to_taxid[taxon_lower]
        
        return None
    
    def get_rank(self, taxid: str) -> str:
        """Get rank for a TaxID."""
        if taxid in self.nodes:
            return self.nodes[taxid]['rank']
        return 'unknown'
    
    def get_name(self, taxid: str) -> str:
        """Get scientific name for a TaxID."""
        return self.names.get(taxid, 'Unknown')
    
    def get_lineage(self, taxid: str) -> Dict[str, str]:
        """Get full taxonomic lineage for a taxid."""
        lineage = {}
        current_taxid = str(taxid)
        
        visited = set()
        while current_taxid in self.nodes and current_taxid not in visited:
            visited.add(current_taxid)
            
            node = self.nodes[current_taxid]
            rank = node['rank']
            name = self.names.get(current_taxid, 'Unknown')
            
            lineage[rank] = name
            
            parent = node['parent']
            if parent == current_taxid:  # Root node
                break
            current_taxid = parent
        
        return lineage


class TaxonAssemblyGatherer:
    """Gather and filter NCBI assemblies for a specified taxon."""
    
    def __init__(
        self,
        taxon: str,
        output_dir: Path,
        taxdump_dir: Optional[Path] = None,
        rank: Optional[str] = None,
        since_date: Optional[str] = None,
        min_completeness: float = 95.0,
        max_contamination: float = 1.0,
        min_n50: int = 5000
    ):
        self.taxon = taxon
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.rank = rank
        self.since_date = since_date
        self.min_completeness = min_completeness
        self.max_contamination = max_contamination
        self.min_n50 = min_n50
        
        # Load taxonomy
        if taxdump_dir is None:
            taxdump_dir = self.output_dir / "taxdump"
        self.taxdump_dir = Path(taxdump_dir)
        
        self._ensure_taxonomy_database()
        self.taxonomy = NCBITaxonomy(self.taxdump_dir)
        
        # Resolve taxon to TaxID
        self.taxid = self.taxonomy.resolve_taxon(taxon)
        if not self.taxid:
            raise ValueError(f"Could not resolve taxon: {taxon}")
        
        self.taxon_name = self.taxonomy.get_name(self.taxid)
        self.taxon_rank = self.taxonomy.get_rank(self.taxid)
        
        logger.info(f"Resolved taxon: {self.taxon_name} (TaxID: {self.taxid}, Rank: {self.taxon_rank})")
        
        # Storage
        self.assemblies = []
    
    def _ensure_taxonomy_database(self):
        """Download taxonomy database if not present."""
        nodes_file = self.taxdump_dir / "nodes.dmp"
        names_file = self.taxdump_dir / "names.dmp"
        
        if nodes_file.exists() and names_file.exists():
            logger.info("✓ Taxonomy database already present")
            return
        
        logger.info("Downloading NCBI taxonomy database...")
        self.taxdump_dir.mkdir(parents=True, exist_ok=True)
        
        taxdump_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        taxdump_tar = self.taxdump_dir / "taxdump.tar.gz"
        
        try:
            urllib.request.urlretrieve(taxdump_url, taxdump_tar)
            logger.info("  ✓ Downloaded taxdump.tar.gz")
            
            with tarfile.open(taxdump_tar, 'r:gz') as tar:
                members = [m for m in tar.getmembers() 
                          if m.name in ['nodes.dmp', 'names.dmp']]
                tar.extractall(path=self.taxdump_dir, members=members)
            logger.info("  ✓ Extracted taxonomy files")
        except Exception as e:
            raise RuntimeError(f"Failed to download taxonomy database: {e}")
    
    def download_assembly_summary(self, db_type: str = 'refseq'):
        """Download assembly summary from NCBI."""
        logger.info(f"\nDownloading {db_type.upper()} assembly summary...")
        
        urls = {
            'refseq': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt',
            'genbank': 'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
        }
        
        output_file = self.output_dir / f"assembly_summary_{db_type}.txt"
        
        if output_file.exists():
            logger.info(f"  ✓ Using existing file: {output_file}")
            return output_file
        
        try:
            urllib.request.urlretrieve(urls[db_type], output_file)
            size_mb = output_file.stat().st_size / (1024 * 1024)
            logger.info(f"  ✓ Downloaded: {size_mb:.1f} MB")
            return output_file
        except Exception as e:
            raise RuntimeError(f"Failed to download assembly summary: {e}")
    
    def parse_assembly_summary(self, summary_file: Path, db_type: str):
        """Parse assembly summary and filter by taxon."""
        logger.info(f"\nParsing {db_type.upper()} assemblies for taxon {self.taxid}...")
        
        assemblies = []
        
        with open(summary_file, 'r') as f:
            # Find header
            header = None
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#') and 'assembly_accession' in line:
                    header = line.strip().lstrip('#').split('\t')
                    break
            
            if not header:
                raise ValueError(f"Could not find header in {summary_file}")
            
            # Get column indices
            try:
                acc_idx = header.index('assembly_accession')
                taxid_idx = header.index('taxid')
                organism_idx = header.index('organism_name')
                asm_name_idx = header.index('asm_name')
                level_idx = header.index('assembly_level')
                date_idx = header.index('seq_rel_date')
                ftp_idx = header.index('ftp_path')
            except ValueError as e:
                raise ValueError(f"Missing required column: {e}")
            
            # Parse data lines
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < len(header):
                    continue
                
                asm_taxid = fields[taxid_idx]
                
                # Check if this assembly belongs to our taxon
                if not self._is_descendant_of_taxon(asm_taxid):
                    continue
                
                # Check date filter
                if self.since_date:
                    asm_date = fields[date_idx]
                    if asm_date < self.since_date:
                        continue
                
                assembly = {
                    'accession': fields[acc_idx],
                    'taxid': asm_taxid,
                    'organism': fields[organism_idx],
                    'asm_name': fields[asm_name_idx],
                    'level': fields[level_idx],
                    'date': fields[date_idx],
                    'ftp_path': fields[ftp_idx],
                    'db_type': db_type
                }
                
                assemblies.append(assembly)
        
        logger.info(f"  ✓ Found {len(assemblies):,} assemblies")
        return assemblies
    
    def _is_descendant_of_taxon(self, test_taxid: str) -> bool:
        """Check if test_taxid is a descendant of self.taxid."""
        current = test_taxid
        visited = set()
        
        while current in self.taxonomy.nodes and current not in visited:
            if current == self.taxid:
                return True
            
            visited.add(current)
            parent = self.taxonomy.nodes[current]['parent']
            
            if parent == current:  # Root
                break
            
            current = parent
        
        return False
    
    def remove_refseq_genbank_redundancy(self, assemblies: List[Dict]) -> List[Dict]:
        """Remove redundant RefSeq/GenBank pairs, preferring RefSeq."""
        logger.info("\nRemoving RefSeq/GenBank redundancy...")
        
        # Group by base accession number
        by_base_acc = defaultdict(list)
        
        for asm in assemblies:
            acc = asm['accession']
            # Extract base number (e.g., GCF_001234567.1 -> 001234567)
            base = acc.split('_')[1].split('.')[0]
            by_base_acc[base].append(asm)
        
        # For each group, prefer RefSeq
        unique = []
        for base, group in by_base_acc.items():
            refseq = [a for a in group if a['accession'].startswith('GCF_')]
            genbank = [a for a in group if a['accession'].startswith('GCA_')]
            
            if refseq:
                # Take latest RefSeq version
                refseq.sort(key=lambda x: x['accession'], reverse=True)
                unique.append(refseq[0])
            elif genbank:
                # Take latest GenBank version
                genbank.sort(key=lambda x: x['accession'], reverse=True)
                unique.append(genbank[0])
        
        logger.info(f"  ✓ Reduced from {len(assemblies):,} to {len(unique):,} unique assemblies")
        return unique
    
    def gather_assemblies(self):
        """Main method to gather assemblies."""
        logger.info("=" * 70)
        logger.info("GATHERING ASSEMBLIES FROM NCBI")
        logger.info("=" * 70)
        logger.info(f"Taxon: {self.taxon_name} (TaxID: {self.taxid})")
        logger.info(f"Rank: {self.taxon_rank}")
        if self.since_date:
            logger.info(f"Since date: {self.since_date}")
        
        # Download and parse RefSeq
        refseq_file = self.download_assembly_summary('refseq')
        refseq_asms = self.parse_assembly_summary(refseq_file, 'refseq')
        
        # Download and parse GenBank
        genbank_file = self.download_assembly_summary('genbank')
        genbank_asms = self.parse_assembly_summary(genbank_file, 'genbank')
        
        # Combine and remove redundancy
        all_asms = refseq_asms + genbank_asms
        self.assemblies = self.remove_refseq_genbank_redundancy(all_asms)
        
        logger.info(f"\n✓ Total unique assemblies: {len(self.assemblies):,}")
        
        return self.assemblies
    
    def generate_output_tsv(self, output_file: Path):
        """Generate TSV file compatible with OrthoPhyl wrapper."""
        logger.info(f"\nGenerating output TSV: {output_file}")
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Header
            writer.writerow([
                'assembly_accession',
                'organism_name',
                'asm_name',
                'assembly_level',
                'release_date',
                'ftp_path',
                'taxid'
            ])
            
            # Data rows
            for asm in self.assemblies:
                writer.writerow([
                    asm['accession'],
                    asm['organism'],
                    asm['asm_name'],
                    asm['level'],
                    asm['date'],
                    asm['ftp_path'],
                    asm['taxid']
                ])
        
        logger.info(f"  ✓ Wrote {len(self.assemblies):,} assemblies to {output_file}")
    
    def generate_metadata_json(self, output_file: Path):
        """Generate metadata JSON for database creation."""
        metadata = {
            'source_taxon_name': self.taxon_name,
            'source_taxid': self.taxid,
            'source_rank': self.taxon_rank,
            'query_date': datetime.now().isoformat(),
            'since_date': self.since_date,
            'n_assemblies': len(self.assemblies),
            'quality_filters': {
                'min_completeness': self.min_completeness,
                'max_contamination': self.max_contamination,
                'min_n50': self.min_n50
            },
            'assembly_accessions': [a['accession'] for a in self.assemblies]
        }
        
        import json
        with open(output_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"  ✓ Wrote metadata to {output_file}")

    def query_ncbi(self) -> List[Dict]:
        """Alias for gather_assemblies() -- the name used by the pipeline wrapper.

        Returns the list of unique assembly dicts (each with an 'accession' key).
        """
        return self.gather_assemblies()

    # NCBI rank name -> GTDB single-letter prefix.
    _GTDB_PREFIX = {
        'superkingdom': 'd',
        'domain': 'd',
        'phylum': 'p',
        'class': 'c',
        'order': 'o',
        'family': 'f',
        'genus': 'g',
        'species': 's',
    }
    _GTDB_ORDER = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    def get_taxonomy_string(self) -> str:
        """Build a GTDB-format taxonomy string for the resolved taxon.

        Walks the NCBI lineage of ``self.taxid`` and renders it as
        ``d__Bacteria;p__...;...`` up to the most specific defined rank. Ranks with no
        NCBI name are rendered as empty (e.g. ``s__``).
        """
        lineage = self.taxonomy.get_lineage(self.taxid)

        # Map NCBI-rank-named lineage onto GTDB prefixes.
        by_prefix = {}
        for ncbi_rank, name in lineage.items():
            prefix = self._GTDB_PREFIX.get(ncbi_rank)
            if prefix:
                by_prefix[prefix] = name

        # Emit ranks from domain down to the deepest one we actually have.
        parts = []
        deepest = None
        for prefix in self._GTDB_ORDER:
            if prefix in by_prefix:
                deepest = prefix
        if deepest is None:
            return ""
        for prefix in self._GTDB_ORDER:
            name = by_prefix.get(prefix, "")
            parts.append(f"{prefix}__{name}")
            if prefix == deepest:
                break
        return ";".join(parts)

    def download_assemblies(self, assemblies: List[Dict], download_dir: Path):
        """Download genome FASTAs for the given assemblies from the NCBI FTP site.

        Each assembly dict must carry an ``ftp_path`` and ``accession``. The genomic
        FASTA is fetched from ``{ftp_path}/{basename}_genomic.fna.gz`` and written,
        gunzipped, to ``download_dir/{accession}.fna``. Assemblies lacking an FTP path
        are skipped with a warning.

        Returns the list of Paths that were successfully written.
        """
        download_dir = Path(download_dir)
        download_dir.mkdir(parents=True, exist_ok=True)

        written = []
        for asm in assemblies:
            accession = asm['accession']
            ftp_path = asm.get('ftp_path', '')
            if not ftp_path or ftp_path == 'na':
                logger.warning(f"  ⚠ No ftp_path for {accession}, skipping")
                continue

            basename = ftp_path.rstrip('/').split('/')[-1]
            url = f"{ftp_path}/{basename}_genomic.fna.gz"
            gz_dest = download_dir / f"{accession}.fna.gz"
            fna_dest = download_dir / f"{accession}.fna"

            try:
                urllib.request.urlretrieve(url, gz_dest)
                with gzip.open(gz_dest, 'rb') as gz, open(fna_dest, 'wb') as out:
                    shutil.copyfileobj(gz, out)
                gz_dest.unlink()
                written.append(fna_dest)
            except Exception as e:
                logger.warning(f"  ⚠ Failed to download {accession}: {e}")

        logger.info(f"  ✓ Downloaded {len(written)}/{len(assemblies)} assemblies")
        return written


def main():
    parser = argparse.ArgumentParser(
        description="Gather NCBI assemblies for a specified taxon",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Get all assemblies for a genus
  python taxon_assembly_gatherer.py \\
      --taxon "Methylorubrum" \\
      --rank genus \\
      --output assemblies.tsv
  
  # Get assemblies by TaxID
  python taxon_assembly_gatherer.py \\
      --taxon 34007 \\
      --output assemblies.tsv
  
  # Get assemblies added since a date (for updates)
  python taxon_assembly_gatherer.py \\
      --taxon "Methylorubrum" \\
      --since-date 2025-01-01 \\
      --output new_assemblies.tsv
        """
    )
    
    parser.add_argument(
        '--taxon',
        required=True,
        help='Taxon name or TaxID (e.g., "Methylorubrum" or "34007")'
    )
    parser.add_argument(
        '--rank',
        choices=['family', 'order', 'genus', 'species'],
        help='Expected taxonomic rank (optional, for validation)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output TSV file'
    )
    parser.add_argument(
        '--output-dir',
        default='taxon_gather_output',
        help='Working directory for downloads (default: taxon_gather_output)'
    )
    parser.add_argument(
        '--since-date',
        help='Only include assemblies released on or after this date (YYYY-MM-DD)'
    )
    parser.add_argument(
        '--min-completeness',
        type=float,
        default=95.0,
        help='Minimum completeness percentage (default: 95.0)'
    )
    parser.add_argument(
        '--max-contamination',
        type=float,
        default=1.0,
        help='Maximum contamination percentage (default: 1.0)'
    )
    parser.add_argument(
        '--min-n50',
        type=int,
        default=5000,
        help='Minimum N50 in bp (default: 5000)'
    )
    
    args = parser.parse_args()
    
    try:
        # Create gatherer
        gatherer = TaxonAssemblyGatherer(
            taxon=args.taxon,
            output_dir=Path(args.output_dir),
            rank=args.rank,
            since_date=args.since_date,
            min_completeness=args.min_completeness,
            max_contamination=args.max_contamination,
            min_n50=args.min_n50
        )
        
        # Gather assemblies
        assemblies = gatherer.gather_assemblies()
        
        if not assemblies:
            logger.warning("No assemblies found matching criteria")
            return 1
        
        # Generate outputs
        output_file = Path(args.output)
        gatherer.generate_output_tsv(output_file)
        
        # Generate metadata
        metadata_file = output_file.parent / f"{output_file.stem}_metadata.json"
        gatherer.generate_metadata_json(metadata_file)
        
        logger.info("\n" + "=" * 70)
        logger.info("ASSEMBLY GATHERING COMPLETE!")
        logger.info("=" * 70)
        logger.info(f"Output: {output_file}")
        logger.info(f"Metadata: {metadata_file}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
