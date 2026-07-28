#!/usr/bin/env python3
"""Quick test to debug assembly summary parsing."""

import sys

summary_file = sys.argv[1] if len(sys.argv) > 1 else "ncbi_stats4/assembly_summary_refseq.txt"

print(f"Testing parsing of: {summary_file}\n")

with open(summary_file, 'r') as f:
    line_num = 0
    header = None
    
    for line in f:
        line_num += 1
        
        # Show first 10 lines
        if line_num <= 10:
            print(f"Line {line_num}: {line[:100]}")
        
        # Try to find header
        if header is None:
            # Skip double-hash comments
            if line.startswith('##'):
                print(f"  -> Skipping comment line {line_num}")
                continue
            
            # Header starts with single # and contains assembly_accession
            if line.startswith('#') and 'assembly_accession' in line:
                header = line.strip().lstrip('#').split('\t')
                print(f"\n✓ Found header at line {line_num}")
                print(f"  Columns: {len(header)}")
                print(f"  First 5 columns: {header[:5]}")
                
                # Find key columns
                try:
                    asm_name_idx = header.index('asm_name')
                    organism_idx = header.index('organism_name')
                    taxid_idx = header.index('taxid')
                    species_taxid_idx = header.index('species_taxid')
                    
                    print(f"\n  Column indices:")
                    print(f"    asm_name: {asm_name_idx}")
                    print(f"    organism_name: {organism_idx}")
                    print(f"    taxid: {taxid_idx}")
                    print(f"    species_taxid: {species_taxid_idx}")
                except ValueError as e:
                    print(f"\n  ✗ Error finding columns: {e}")
                    sys.exit(1)
                
                break
    
    if header is None:
        print("\n✗ Could not find header!")
        sys.exit(1)
    
    # Try to parse first data line
    print(f"\nParsing first data line...")
    for line in f:
        line_num += 1
        
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        print(f"\n✓ First data line (line {line_num}):")
        print(f"  Fields: {len(fields)}")
        print(f"  asm_name: {fields[asm_name_idx]}")
        print(f"  organism: {fields[organism_idx]}")
        print(f"  taxid: {fields[taxid_idx]}")
        print(f"  species_taxid: {fields[species_taxid_idx]}")
        
        break
    
    # Count total lines
    count = 1  # Already read one data line
    for line in f:
        if not line.startswith('#'):
            count += 1
    
    print(f"\n✓ Total data lines: {count:,}")
