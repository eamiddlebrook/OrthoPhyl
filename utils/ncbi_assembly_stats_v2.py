#!/usr/bin/env python3
"""
NCBI Assembly Statistics v2 - With Full Taxonomic Lineages

This enhanced version:
1. Downloads NCBI taxonomy database (taxdump)
2. Parses taxonomic lineages for all assemblies
3. Generates counts by ALL taxonomic levels (phylum, class, order, family, genus, species)
4. Creates comprehensive summary tables

Usage:
    python ncbi_assembly_stats_v2.py --output-dir ncbi_stats/
    
    # With filters
    python ncbi_assembly_stats_v2.py --output-dir ncbi_stats/ --min-assemblies 10
    
    # Skip taxonomy download if already present
    python ncbi_assembly_stats_v2.py --output-dir ncbi_stats/ --skip-taxdump-download

Requirements:
    - Python 3.7+
    - Internet connection for downloading from NCBI
    - ~1GB disk space for taxonomy database and assembly summaries

Author: HGTool/OrthoPhyl Project
Date: 2025
"""

import os
import sys
import argparse
import urllib.request
import tarfile
import gzip
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set
import csv
from datetime import datetime

# NCBI FTP URLs
NCBI_URLS = {
    'refseq': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt',
    'genbank': 'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt',
    'taxdump': 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
}

# Taxonomic ranks we care about
RANKS_OF_INTEREST = ['phylum', 'class', 'order', 'family', 'genus', 'species']


class NCBITaxonomy:
    """Parse and query NCBI taxonomy database."""
    
    def __init__(self, taxdump_dir: Path):
        self.taxdump_dir = Path(taxdump_dir)
        
        # Taxonomy data structures
        self.nodes = {}  # taxid -> {'parent': parent_taxid, 'rank': rank}
        self.names = {}  # taxid -> scientific_name
        
        print("\n" + "=" * 70)
        print("LOADING NCBI TAXONOMY DATABASE")
        print("=" * 70)
        
        self._load_nodes()
        self._load_names()
        
        print(f"\n✓ Loaded {len(self.nodes):,} taxonomy nodes")
        print(f"✓ Loaded {len(self.names):,} taxonomy names")
    
    def _load_nodes(self):
        """Load nodes.dmp (taxonomy tree structure)."""
        nodes_file = self.taxdump_dir / "nodes.dmp"
        
        if not nodes_file.exists():
            raise FileNotFoundError(f"nodes.dmp not found in {self.taxdump_dir}")
        
        print(f"\nLoading nodes.dmp...")
        
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
        
        print(f"  ✓ Loaded {len(self.nodes):,} nodes")
    
    def _load_names(self):
        """Load names.dmp (taxonomy names)."""
        names_file = self.taxdump_dir / "names.dmp"
        
        if not names_file.exists():
            raise FileNotFoundError(f"names.dmp not found in {self.taxdump_dir}")
        
        print(f"\nLoading names.dmp...")
        
        with open(names_file, 'r') as f:
            for line in f:
                parts = [p.strip() for p in line.split('|')]
                taxid = parts[0]
                name = parts[1]
                name_class = parts[3]
                
                # Only keep scientific names
                if name_class == 'scientific name':
                    self.names[taxid] = name
        
        print(f"  ✓ Loaded {len(self.names):,} scientific names")
    
    def get_lineage(self, taxid: str) -> Dict[str, str]:
        """Get full taxonomic lineage for a taxid.
        
        Returns dict with rank -> name mappings.
        """
        lineage = {}
        current_taxid = str(taxid)
        
        # Walk up the tree
        visited = set()
        while current_taxid in self.nodes and current_taxid not in visited:
            visited.add(current_taxid)
            
            node = self.nodes[current_taxid]
            rank = node['rank']
            name = self.names.get(current_taxid, 'Unknown')
            
            # Store if it's a rank we care about
            if rank in RANKS_OF_INTEREST:
                lineage[rank] = name
            
            # Move to parent
            parent = node['parent']
            if parent == current_taxid:  # Root node
                break
            current_taxid = parent
        
        return lineage
    
    def get_rank_name(self, taxid: str, rank: str) -> str:
        """Get the name at a specific rank for a taxid."""
        lineage = self.get_lineage(taxid)
        return lineage.get(rank, 'Unknown')


class NCBIAssemblyStatsV2:
    """Download and analyze NCBI bacterial assembly statistics with full taxonomy."""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Storage for assembly data
        self.assemblies = {
            'refseq': [],
            'genbank': []
        }
        
        # Taxonomic summaries by rank
        self.taxonomy_counts = {
            'refseq': {rank: defaultdict(int) for rank in RANKS_OF_INTEREST},
            'genbank': {rank: defaultdict(int) for rank in RANKS_OF_INTEREST},
            'combined': {rank: defaultdict(int) for rank in RANKS_OF_INTEREST}
        }
        
        # Taxonomy database
        self.taxonomy = None
    
    def download_taxdump(self, skip_if_exists: bool = False):
        """Download and extract NCBI taxonomy database."""
        print("\n" + "=" * 70)
        print("DOWNLOADING NCBI TAXONOMY DATABASE")
        print("=" * 70)
        
        taxdump_tar = self.output_dir / "taxdump.tar.gz"
        
        # Check if already extracted
        if (self.output_dir / "nodes.dmp").exists() and (self.output_dir / "names.dmp").exists():
            if skip_if_exists:
                print("\n✓ Taxonomy database already present")
                return
            else:
                print("\n⚠ Taxonomy database exists but will re-download")
        
        # Download
        if not taxdump_tar.exists():
            print(f"\nDownloading taxdump.tar.gz...")
            print(f"  URL: {NCBI_URLS['taxdump']}")
            print(f"  Output: {taxdump_tar}")
            print(f"  Size: ~50 MB (this may take a minute)")
            
            try:
                urllib.request.urlretrieve(NCBI_URLS['taxdump'], taxdump_tar)
                size_mb = taxdump_tar.stat().st_size / (1024 * 1024)
                print(f"  ✓ Downloaded: {size_mb:.1f} MB")
            except Exception as e:
                print(f"  ✗ Download failed: {e}")
                raise
        else:
            print(f"\n✓ taxdump.tar.gz already downloaded")
        
        # Extract
        print(f"\nExtracting taxonomy files...")
        try:
            with tarfile.open(taxdump_tar, 'r:gz') as tar:
                # Extract only the files we need
                members = [m for m in tar.getmembers() 
                          if m.name in ['nodes.dmp', 'names.dmp']]
                tar.extractall(path=self.output_dir, members=members)
            print(f"  ✓ Extracted nodes.dmp and names.dmp")
        except Exception as e:
            print(f"  ✗ Extraction failed: {e}")
            raise
        
        print("\n✓ Taxonomy database ready")
    
    def load_taxonomy(self):
        """Load NCBI taxonomy database."""
        self.taxonomy = NCBITaxonomy(self.output_dir)
    
    def download_assembly_summaries(self):
        """Download assembly summary files from NCBI."""
        print("\n" + "=" * 70)
        print("DOWNLOADING ASSEMBLY SUMMARIES FROM NCBI")
        print("=" * 70)
        
        for db_type, url in NCBI_URLS.items():
            if db_type == 'taxdump':
                continue
            
            output_file = self.output_dir / f"assembly_summary_{db_type}.txt"
            
            if output_file.exists():
                print(f"\n✓ {db_type.upper()} summary already exists: {output_file}")
                print(f"  Delete file to re-download")
                continue
            
            print(f"\nDownloading {db_type.upper()} assembly summary...")
            print(f"  URL: {url}")
            print(f"  Output: {output_file}")
            
            try:
                urllib.request.urlretrieve(url, output_file)
                size_mb = output_file.stat().st_size / (1024 * 1024)
                print(f"  ✓ Downloaded: {size_mb:.1f} MB")
            except Exception as e:
                print(f"  ✗ Download failed: {e}")
                raise
        
        print("\n✓ All assembly summaries downloaded")
    
    def parse_assembly_summaries(self):
        """Parse assembly summary files and extract taxonomy information."""
        print("\n" + "=" * 70)
        print("PARSING ASSEMBLY SUMMARIES")
        print("=" * 70)
        
        for db_type in ['refseq', 'genbank']:
            summary_file = self.output_dir / f"assembly_summary_{db_type}.txt"
            
            if not summary_file.exists():
                print(f"\n⚠ Skipping {db_type.upper()}: file not found")
                continue
            
            print(f"\nParsing {db_type.upper()} assemblies...")
            
            assemblies = []
            with open(summary_file, 'r') as f:
                # Find header line (starts with # but contains column names)
                header = None
                for line in f:
                    # Skip comment lines that start with ## (double hash)
                    if line.startswith('##'):
                        continue
                    # Header line starts with single # and contains assembly_accession
                    if line.startswith('#') and 'assembly_accession' in line:
                        # Remove leading # and parse
                        header = line.strip().lstrip('#').split('\t')
                        break
                
                if not header:
                    print(f"  ✗ Could not find header line")
                    continue
                
                # Find column indices
                try:
                    asm_name_idx = header.index('asm_name')
                    organism_idx = header.index('organism_name')
                    taxid_idx = header.index('taxid')
                    species_taxid_idx = header.index('species_taxid')
                except ValueError as e:
                    print(f"  ✗ Error finding columns: {e}")
                    continue
                
                # Parse data lines
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < len(header):
                        continue
                    
                    assembly = {
                        'asm_name': fields[asm_name_idx],
                        'organism': fields[organism_idx],
                        'taxid': fields[taxid_idx],
                        'species_taxid': fields[species_taxid_idx],
                        'db_type': db_type
                    }
                    
                    assemblies.append(assembly)
            
            self.assemblies[db_type] = assemblies
            print(f"  ✓ Parsed {len(assemblies):,} {db_type.upper()} assemblies")
        
        total = sum(len(v) for v in self.assemblies.values())
        print(f"\n✓ Total assemblies parsed: {total:,}")
    
    def summarize_by_taxonomy(self):
        """Summarize assembly counts by all taxonomic levels using NCBI taxonomy."""
        print("\n" + "=" * 70)
        print("SUMMARIZING BY TAXONOMY (ALL RANKS)")
        print("=" * 70)
        
        if not self.taxonomy:
            raise RuntimeError("Taxonomy database not loaded. Call load_taxonomy() first.")
        
        for db_type in ['refseq', 'genbank']:
            if not self.assemblies[db_type]:
                continue
            
            print(f"\nProcessing {db_type.upper()} assemblies...")
            
            processed = 0
            for assembly in self.assemblies[db_type]:
                taxid = assembly['taxid']
                
                # Get full lineage
                lineage = self.taxonomy.get_lineage(taxid)
                
                # Count at each rank
                for rank in RANKS_OF_INTEREST:
                    if rank in lineage:
                        taxon_name = lineage[rank]
                        self.taxonomy_counts[db_type][rank][taxon_name] += 1
                        self.taxonomy_counts['combined'][rank][taxon_name] += 1
                
                processed += 1
                if processed % 10000 == 0:
                    print(f"  Processed {processed:,} assemblies...")
            
            print(f"  ✓ Processed {len(self.assemblies[db_type]):,} assemblies")
        
        print("\n✓ Taxonomy summary complete")
        
        # Print summary stats
        print("\nTaxonomic diversity:")
        for rank in RANKS_OF_INTEREST:
            count = len(self.taxonomy_counts['combined'][rank])
            print(f"  {rank.capitalize():<12s}: {count:>6,} unique taxa")
    
    def generate_summary_tables(self, min_assemblies: int = 1):
        """Generate summary tables by taxonomic level."""
        print("\n" + "=" * 70)
        print("GENERATING SUMMARY TABLES")
        print("=" * 70)
        
        for rank in RANKS_OF_INTEREST:
            self._generate_rank_table(rank, min_assemblies)
        
        print("\n✓ Summary tables generated")
    
    def _generate_rank_table(self, rank: str, min_assemblies: int):
        """Generate summary table for a specific taxonomic rank."""
        output_file = self.output_dir / f"assembly_counts_by_{rank}.tsv"
        
        print(f"\nGenerating {rank} table...")
        
        # Collect all taxa at this rank
        all_taxa = set()
        for db_type in ['refseq', 'genbank']:
            all_taxa.update(self.taxonomy_counts[db_type][rank].keys())
        
        # Build table rows
        rows = []
        for taxon in sorted(all_taxa):
            refseq_count = self.taxonomy_counts['refseq'][rank].get(taxon, 0)
            genbank_count = self.taxonomy_counts['genbank'][rank].get(taxon, 0)
            total = refseq_count + genbank_count
            
            # Filter by minimum assemblies
            if total < min_assemblies:
                continue
            
            rows.append({
                rank: taxon,
                'refseq': refseq_count,
                'genbank': genbank_count,
                'total': total
            })
        
        # Sort by total count (descending)
        rows.sort(key=lambda x: x['total'], reverse=True)
        
        # Write table
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[rank, 'refseq', 'genbank', 'total'], delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        
        print(f"  ✓ Wrote {len(rows):,} {rank} entries to {output_file}")
        
        # Print top 10
        print(f"\n  Top 10 {rank} by total assemblies:")
        for i, row in enumerate(rows[:10], 1):
            print(f"    {i:2d}. {row[rank]:<50s} RefSeq: {row['refseq']:>6,}  GenBank: {row['genbank']:>6,}  Total: {row['total']:>7,}")
    
    def generate_summary_report(self):
        """Generate overall summary report."""
        report_file = self.output_dir / "summary_report.txt"
        
        print(f"\nGenerating summary report...")
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("NCBI BACTERIAL ASSEMBLY STATISTICS (WITH FULL TAXONOMY)\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Overall counts
            f.write("OVERALL COUNTS:\n")
            f.write("-" * 70 + "\n")
            refseq_total = len(self.assemblies['refseq'])
            genbank_total = len(self.assemblies['genbank'])
            f.write(f"RefSeq assemblies:  {refseq_total:>10,}\n")
            f.write(f"GenBank assemblies: {genbank_total:>10,}\n")
            f.write(f"Total assemblies:   {refseq_total + genbank_total:>10,}\n\n")
            
            # Taxonomy counts by rank
            f.write("\nTAXONOMIC DIVERSITY:\n")
            f.write("-" * 70 + "\n")
            
            for rank in RANKS_OF_INTEREST:
                refseq_taxa = len(self.taxonomy_counts['refseq'][rank])
                genbank_taxa = len(self.taxonomy_counts['genbank'][rank])
                combined_taxa = len(self.taxonomy_counts['combined'][rank])
                
                f.write(f"\n{rank.upper()}:\n")
                f.write(f"  Unique {rank} in RefSeq:  {refseq_taxa:>6,}\n")
                f.write(f"  Unique {rank} in GenBank: {genbank_taxa:>6,}\n")
                f.write(f"  Unique {rank} combined:   {combined_taxa:>6,}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("GENERATED FILES:\n")
            f.write("-" * 70 + "\n")
            for rank in RANKS_OF_INTEREST:
                f.write(f"  - assembly_counts_by_{rank}.tsv\n")
            f.write("\n" + "=" * 70 + "\n")
        
        print(f"  ✓ Summary report: {report_file}")
    
    def run(self, skip_taxdump_download: bool = False, min_assemblies: int = 1):
        """Run complete analysis pipeline."""
        # Step 1: Download taxonomy database
        self.download_taxdump(skip_if_exists=skip_taxdump_download)
        
        # Step 2: Load taxonomy
        self.load_taxonomy()
        
        # Step 3: Download assembly summaries
        self.download_assembly_summaries()
        
        # Step 4: Parse assemblies
        self.parse_assembly_summaries()
        
        # Step 5: Summarize by taxonomy
        self.summarize_by_taxonomy()
        
        # Step 6: Generate tables
        self.generate_summary_tables(min_assemblies=min_assemblies)
        
        # Step 7: Generate report
        self.generate_summary_report()
        
        print("\n" + "=" * 70)
        print("ANALYSIS COMPLETE!")
        print("=" * 70)
        print(f"\nOutput directory: {self.output_dir}")
        print("\nGenerated files:")
        print("  - assembly_summary_refseq.txt")
        print("  - assembly_summary_genbank.txt")
        for rank in RANKS_OF_INTEREST:
            print(f"  - assembly_counts_by_{rank}.tsv")
        print("  - summary_report.txt")
        print("\nTaxonomy database files:")
        print("  - taxdump.tar.gz")
        print("  - nodes.dmp")
        print("  - names.dmp")


def main():
    parser = argparse.ArgumentParser(
        description="Download and analyze NCBI bacterial assembly statistics with full taxonomy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (downloads everything)
  python ncbi_assembly_stats_v2.py --output-dir ncbi_stats/
  
  # Filter to taxa with at least 10 assemblies
  python ncbi_assembly_stats_v2.py --output-dir ncbi_stats/ --min-assemblies 10
  
  # Skip taxonomy download if already present
  python ncbi_assembly_stats_v2.py --output-dir ncbi_stats/ --skip-taxdump-download

Features:
  - Downloads NCBI taxonomy database (taxdump)
  - Provides assembly counts for ALL taxonomic ranks:
    * Phylum
    * Class
    * Order
    * Family
    * Genus
    * Species
  - Generates separate TSV files for each rank
  - Shows RefSeq, GenBank, and combined counts
        """
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for downloads and results'
    )
    parser.add_argument(
        '--skip-taxdump-download',
        action='store_true',
        help='Skip downloading taxonomy database if already present'
    )
    parser.add_argument(
        '--min-assemblies',
        type=int,
        default=1,
        help='Minimum number of assemblies to include in summary tables (default: 1)'
    )
    
    args = parser.parse_args()
    
    # Create analyzer
    analyzer = NCBIAssemblyStatsV2(output_dir=args.output_dir)
    
    # Run analysis
    try:
        analyzer.run(
            skip_taxdump_download=args.skip_taxdump_download,
            min_assemblies=args.min_assemblies
        )
        return 0
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
