#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_redundant_genomes.py

Remove redundant genome-assembly versions from a list of genome stems, without the
locale-dependent silent data loss of the bash `filter_for_redundancy` step it replaces
(unanchored `grep $I | sort | tail -n 1` over a `sed`-truncated accession number, which
collapses distinct user filenames that happen to share a numeric substring).

Only names that actually look like NCBI accessions (GCF_/GCA_ + digits + version) carry a
"redundant version" concept -- only those are de-versioned. Every other name (a user's own
filename stem) is passed through untouched: there is nothing to redundancy-filter about
`isolate_1` vs `isolate_11`, and treating them as related is exactly the bug being fixed.

Usage:
    filter_redundant_genomes.py --stems-file names.txt
    (or import partition_by_accession / dedupe_accessions directly for unit testing)
"""
import argparse
import re
import sys
from collections import defaultdict

# NCBI assembly accession: GCF_/GCA_ + 9 digits + ".version". Anything else is left alone.
ACCESSION_RE = re.compile(r'^(GC[AF])_(\d+)\.(\d+)$')


def partition_by_accession(names):
    """
    Split `names` into (accession_like, passthrough) based on ACCESSION_RE.

    accession_like: names matching the NCBI accession pattern.
    passthrough:    everything else, unchanged and in original order.
    """
    accession_like = [n for n in names if ACCESSION_RE.match(n)]
    passthrough = [n for n in names if not ACCESSION_RE.match(n)]
    return accession_like, passthrough


def dedupe_accessions(names):
    """
    Given accession-like names (GCF_/GCA_ + base number + version), keep exactly one
    per base number: prefer GCF over GCA, then the highest version number.

    Pure set/dict based -- no grep, no locale dependence.
    """
    by_base = defaultdict(list)
    for name in names:
        m = ACCESSION_RE.match(name)
        prefix, base, version = m.group(1), m.group(2), int(m.group(3))
        by_base[base].append((prefix, version, name))

    kept = []
    for base, group in by_base.items():
        refseq = [g for g in group if g[0] == 'GCF']
        genbank = [g for g in group if g[0] == 'GCA']
        pool = refseq if refseq else genbank
        pool.sort(key=lambda g: g[1], reverse=True)  # highest version first
        kept.append(pool[0][2])
    return kept


def filter_redundant(names):
    """
    Top-level: de-duplicate accession-like names by (base number, prefer GCF, highest
    version); pass every other name through untouched. Returns a sorted list.
    """
    accession_like, passthrough = partition_by_accession(names)
    kept = dedupe_accessions(accession_like) + passthrough
    return sorted(set(kept))


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Remove redundant NCBI-accession genome versions from a list of "
                    "genome stems (one per line), passing non-accession names through "
                    "unchanged.")
    p.add_argument("--stems-file", required=True,
                   help="File of genome stems, one per line (no extension).")
    args = p.parse_args(argv)

    with open(args.stems_file) as fh:
        names = [line.strip() for line in fh if line.strip()]

    kept = filter_redundant(names)
    for name in kept:
        print(name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
