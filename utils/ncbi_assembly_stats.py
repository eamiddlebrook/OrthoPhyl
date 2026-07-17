#!/usr/bin/env python3
"""
NCBI Assembly Statistics - Download and Summarize Bacterial Assemblies

This script:
1. Downloads assembly summary files from NCBI for bacteria
2. Parses taxonomic information
3. Generates counts by taxonomic level (phylum, class, order, family, genus, species)
4. Creates summary tables showing RefSeq and GenBank assembly counts

Usage:
    python ncbi_assembly_stats.py --output-dir ncbi_stats/
    
    # With filters
    python ncbi_assembly_stats.py --output-dir ncbi_stats/ --min-assemblies 10
    
    # Download only (skip analysis)
    python ncbi_assembly_stats.py --output-dir ncbi_stats/ --download-only

Requirements:
    - Python 3.7+
    - Internet connection for downloading from NCBI
    - ~500MB disk space for assembly summary files

Author: HGTool/OrthoPhyl Project
Date: 2025
"""

import os
import sys
import argparse
import urllib.request
import gzip
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple
import csv
from datetime import datetime

# NCBI FTP URLs for assembly summaries
NCBI_URLS = {
    'refseq': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt',
    'genbank': 'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
}

# Taxonomic rank order
RANK_ORDER = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
RANK_PREFIXES = {
    'superkingdom': 'd__',
    'phylum': 'p__',
    'class': 'c__',
    'order': 'o__',
    'family': 'f__',
    'genus': 'g__',
    'species': 's__'
}


class NCBIAssemblyStats:
    """Download and analyze NCBI bacterial assembly statistics."""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Storage for assembly data
        self.assemblies = {
            'refseq': [],
            'genbank': []
        }
        
        # Taxonomic summaries
        self.taxonomy_counts = {
            'refseq': defaultdict(lambda: defaultdict(int)),
            'genbank': defaultdict(lambda: defaultdict(int)),
            'combined': defaultdict(lambda: defaultdict(int))
        }
    
    def download_assembly_summaries(self):
        """Download assembly summary files from NCBI."""
        print("=" * 70)
        print("DOWNLOADING ASSEMBLY SUMMARIES FROM NCBI")
        print("=" * 70)
        
        for db_type, url in NCBI_URLS.items():
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
                
                # Check file size
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
                    if line.startswith('# assembly_accession'):
                        # This is the header line
                        header = line.strip().lstrip('#').strip().split('\t')
                        break
                    elif line.startswith('##'):
                        # Comment line, skip
                        continue
                    elif line.startswith('#'):
                        # Might be header, try it
                        header = line.strip().lstrip('#').strip().split('\t')
                        break
                
                if not header:
                    print(f"  ✗ Could not find header line")
                    continue
                
                # Find column indices
                try:
                    # Note: NCBI uses '#assembly_accession' for first column (with #)
                    # But after split, the # is part of the column name
                    asm_name_idx = header.index('asm_name')
                    organism_idx = header.index('organism_name')
                    taxid_idx = header.index('taxid')
                    species_taxid_idx = header.index('species_taxid')
                    assembly_level_idx = header.index('assembly_level')
                    release_type_idx = header.index('release_type')
                    genome_rep_idx = header.index('genome_rep')
                    ftp_path_idx = header.index('ftp_path')
                except ValueError as e:
                    print(f"  ✗ Error finding columns: {e}")
                    print(f"  Available columns: {', '.join(header[:10])}...")
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
                        'assembly_level': fields[assembly_level_idx],
                        'release_type': fields[release_type_idx],
                        'genome_rep': fields[genome_rep_idx],
                        'ftp_path': fields[ftp_path_idx],
                        'db_type': db_type
                    }
                    
                    assemblies.append(assembly)
            
            self.assemblies[db_type] = assemblies
            print(f"  ✓ Parsed {len(assemblies):,} {db_type.upper()} assemblies")
        
        total = sum(len(v) for v in self.assemblies.values())
        print(f"\n✓ Total assemblies parsed: {total:,}")
    
    def fetch_taxonomy_from_ncbi(self, taxid: str) -> Dict[str, str]:
        """Fetch full taxonomic lineage from NCBI Taxonomy.
        
        Note: This is a placeholder. In practice, you'd want to:
        1. Download NCBI taxonomy dump files
        2. Parse taxdump files (names.dmp, nodes.dmp)
        3. Build taxonomy lookup table
        
        For now, we'll extract what we can from organism names.
        """
        # This is a simplified version - just parse organism name
        return self._parse_organism_name(taxid)
    
    def _parse_organism_name(self, organism_name: str) -> Dict[str, str]:
        """Extract taxonomic information from organism name.
        
        This is simplified - ideally use NCBI taxonomy database.
        """
        taxonomy = {}
        
        # For bacteria, organism name is typically: Genus species strain
        parts = organism_name.split()
        
        if len(parts) >= 1:
            taxonomy['genus'] = parts[0]
        
        if len(parts) >= 2:
            taxonomy['species'] = f"{parts[0]} {parts[1]}"
        
        # We can't reliably get higher taxonomy from name alone
        # Would need NCBI taxonomy database for full lineage
        
        return taxonomy
    
    def summarize_by_taxonomy(self):
        """Summarize assembly counts by taxonomic levels."""
        print("\n" + "=" * 70)
        print("SUMMARIZING BY TAXONOMY")
        print("=" * 70)
        
        print("\nNote: Full taxonomic lineages require NCBI taxonomy database.")
        print("For now, summarizing by genus and species from organism names.")
        
        for db_type in ['refseq', 'genbank']:
            if not self.assemblies[db_type]:
                continue
            
            print(f"\nProcessing {db_type.upper()} assemblies...")
            
            for assembly in self.assemblies[db_type]:
                organism = assembly['organism']
                
                # Parse genus and species from organism name
                taxonomy = self._parse_organism_name(organism)
                
                # Count by genus
                if 'genus' in taxonomy:
                    genus = taxonomy['genus']
                    self.taxonomy_counts[db_type]['genus'][genus] += 1
                    self.taxonomy_counts['combined']['genus'][genus] += 1
                
                # Count by species
                if 'species' in taxonomy:
                    species = taxonomy['species']
                    self.taxonomy_counts[db_type]['species'][species] += 1
                    self.taxonomy_counts['combined']['species'][species] += 1
            
            print(f"  ✓ Processed {len(self.assemblies[db_type]):,} assemblies")
        
        print("\n✓ Taxonomy summary complete")
    
    def generate_summary_tables(self, min_assemblies: int = 1):
        """Generate summary tables by taxonomic level."""
        print("\n" + "=" * 70)
        print("GENERATING SUMMARY TABLES")
        print("=" * 70)
        
        for rank in ['genus', 'species']:
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
            print(f"    {i:2d}. {row[rank]:<40s} RefSeq: {row['refseq']:>6,}  GenBank: {row['genbank']:>6,}  Total: {row['total']:>7,}")
    
    def generate_summary_report(self):
        """Generate overall summary report."""
        report_file = self.output_dir / "summary_report.txt"
        
        print(f"\nGenerating summary report...")
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("NCBI BACTERIAL ASSEMBLY STATISTICS\n")
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
            
            # Taxonomy counts
            for rank in ['genus', 'species']:
                f.write(f"\n{rank.upper()} COUNTS:\n")
                f.write("-" * 70 + "\n")
                
                refseq_taxa = len(self.taxonomy_counts['refseq'][rank])
                genbank_taxa = len(self.taxonomy_counts['genbank'][rank])
                combined_taxa = len(self.taxonomy_counts['combined'][rank])
                
                f.write(f"Unique {rank} in RefSeq:  {refseq_taxa:>6,}\n")
                f.write(f"Unique {rank} in GenBank: {genbank_taxa:>6,}\n")
                f.write(f"Unique {rank} combined:   {combined_taxa:>6,}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("See assembly_counts_by_*.tsv files for detailed breakdowns\n")
            f.write("=" * 70 + "\n")
        
        print(f"  ✓ Summary report: {report_file}")
    
    def run(self, download_only: bool = False, min_assemblies: int = 1):
        """Run complete analysis pipeline."""
        # Step 1: Download
        self.download_assembly_summaries()
        
        if download_only:
            print("\n✓ Download complete (skipping analysis)")
            return
        
        # Step 2: Parse
        self.parse_assembly_summaries()
        
        # Step 3: Summarize
        self.summarize_by_taxonomy()
        
        # Step 4: Generate tables
        self.generate_summary_tables(min_assemblies=min_assemblies)
        
        # Step 5: Generate report
        self.generate_summary_report()
        
        print("\n" + "=" * 70)
        print("ANALYSIS COMPLETE!")
        print("=" * 70)
        print(f"\nOutput directory: {self.output_dir}")
        print("\nGenerated files:")
        print("  - assembly_summary_refseq.txt")
        print("  - assembly_summary_genbank.txt")
        print("  - assembly_counts_by_genus.tsv")
        print("  - assembly_counts_by_species.tsv")
        print("  - summary_report.txt")


def main():
    parser = argparse.ArgumentParser(
        description="Download and analyze NCBI bacterial assembly statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python ncbi_assembly_stats.py --output-dir ncbi_stats/
  
  # Filter to taxa with at least 10 assemblies
  python ncbi_assembly_stats.py --output-dir ncbi_stats/ --min-assemblies 10
  
  # Download only (skip analysis)
  python ncbi_assembly_stats.py --output-dir ncbi_stats/ --download-only

Note:
  This script provides genus and species level summaries from organism names.
  For complete taxonomic lineages (phylum, class, order, family), you would need
  to download and parse NCBI's taxonomy database (taxdump files).
  
  See: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        """
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for downloads and results'
    )
    parser.add_argument(
        '--download-only',
        action='store_true',
        help='Download assembly summaries only (skip analysis)'
    )
    parser.add_argument(
        '--min-assemblies',
        type=int,
        default=1,
        help='Minimum number of assemblies to include in summary tables (default: 1)'
    )
    
    args = parser.parse_args()
    
    # Create analyzer
    analyzer = NCBIAssemblyStats(output_dir=args.output_dir)
    
    # Run analysis
    try:
        analyzer.run(
            download_only=args.download_only,
            min_assemblies=args.min_assemblies
        )
        return 0
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())