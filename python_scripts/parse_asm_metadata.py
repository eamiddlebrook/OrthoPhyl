#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_asm_metadata.py

Fast replacements for the metadata join/filter steps of
utils/gather_filter_asms.sh. These were bash one-liners whose cost blows up on
large taxa (100k+ assemblies); each subcommand here is a single linear pass
using dicts/sets instead.

What was slow, and why:

  merge-geoloc   <- merge_metadata_geoloc()
      The bash version was a nested `while read` loop that re-`cat`ed the ENTIRE
      biosample geoloc file once per assembly row: O(n*m) shell-loop iterations
      plus one fork per row. At n=m=100k that is ~1e10 iterations. This is just a
      hash join on the BioSample accession -- O(n+m) with a dict.

  filter-by-acc  <- all_sample_metadata()
      The bash version ran `grep -f all_asm_acc.GCF` (and .GCA), i.e. ~100k
      patterns matched as REGEXES against every line, twice. Set membership on
      the accession column is O(1) per row.

  exclude-downloaded  <- get_non_datasets_assemblies()
      The bash version ran `grep -vf assemblies_datasets_uniq.names`, again ~100k
      regex patterns, to drop already-downloaded assemblies.

Behaviour is preserved bit-for-bit against the bash it replaces, including two
quirks that are load-bearing downstream:

  * merge-geoloc is an INNER join. A metadata row whose BioSample has no geoloc
    entry is DROPPED (the bash loop emitted nothing for it). A row matching
    several geoloc entries is emitted once per match, in file order.
  * merge-geoloc collapses whitespace. The bash `echo -e $BS'\\t'$BLAH'\\t'$BLAH1`
    left $BLAH unquoted, so xtract's tabs became single spaces. Downstream awk
    splits on any whitespace, so the collapsed form is what the column indices in
    all_sample_metadata()/filter_asm_by_taxCheck()/annotate_tree_by_country.py
    are counted against.

One deliberate tightening: accessions are matched by EXACT value in the accession
column rather than as a substring of the whole line. `grep -f` matched an
accession anywhere on the line (including inside the FTP path columns) and, being
a substring match, `GCF_000123456.1` also matched `GCF_000123456.10`. Exact
column matching is both faster and removes that false positive.

Field layout of All-<taxon>.assembly.BS_to_meta (from xtract, 1-based):
    1 BioSampleAccn  2 RefSeq  3 Genbank  4 SpeciesName  5 Sub_value
    6 FtpPath_GenBank  7 FtpPath_RefSeq  8 Taxid  9 taxonomy-check-status
    10 ExclFromRefSeq
...and after merge-geoloc appends the geolocation: 11 geo_loc_name.

Usage:
    parse_asm_metadata.py merge-geoloc --meta F --geoloc F --out F
    parse_asm_metadata.py filter-by-acc --meta F --acc-list F --out F \\
        [--out-gcf F] [--out-gca F]
    parse_asm_metadata.py exclude-downloaded --meta F --downloaded F --out F
"""

import argparse
import os
import sys


# --------------------------------------------------------------------------- #
# Readers                                                                     #
# --------------------------------------------------------------------------- #

def iter_rows(path):
    """Yield (raw_line, fields) for each non-blank line of a whitespace table.

    Fields are split on runs of any whitespace, matching how bash `read` and
    awk's default FS treat these files (xtract emits tabs, but earlier `sed`
    passes and unquoted echos mean tabs and spaces are used interchangeably).
    ``raw_line`` keeps the original separators for pass-through subcommands.
    """
    if not path or not os.path.exists(path):
        return
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            yield line, line.split()


def read_acc_set(path):
    """Read a one-accession-per-line file into a set. Missing file -> empty set."""
    accs = set()
    if not path or not os.path.exists(path):
        return accs
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            acc = line.strip()
            if acc:
                accs.add(acc)
    return accs


def read_geoloc_map(path):
    """Map BioSample accession -> [geolocation, ...] preserving file order.

    Input is All-<taxon>.biosample.BS_to_Geoloc: ``accession<TAB>geo_loc_name``.
    A list (not a single value) because the bash join emitted one output row per
    matching geoloc row.
    """
    geo = {}
    for _raw, fields in iter_rows(path):
        if not fields:
            continue
        key = fields[0]
        # Everything after the accession is the location; re-join with single
        #   spaces exactly as the unquoted `echo $BLAH1` did.
        value = " ".join(fields[1:])
        geo.setdefault(key, []).append(value)
    return geo


def split_accessions(accs):
    """Split an accession set into (refseq, genbank) by GCF_/GCA_ prefix."""
    gcf = {a for a in accs if a.startswith("GCF")}
    gca = {a for a in accs if a.startswith("GCA")}
    return gcf, gca


# --------------------------------------------------------------------------- #
# Pure transforms (unit-tested directly)                                      #
# --------------------------------------------------------------------------- #

def merge_geoloc(rows, geoloc_map):
    """Hash-join metadata rows against geoloc_map on the BioSample accession.

    ``rows`` is an iterable of (raw_line, fields). Yields joined output lines:
    ``BioSample<TAB>rest-of-metadata<TAB>geolocation``.

    INNER join, deliberately: rows with no geoloc match are dropped, matching the
    bash loop this replaces. Returns via generator; call merge_geoloc_stats for
    the drop count.
    """
    for _raw, fields in rows:
        if not fields:
            continue
        bs = fields[0]
        rest = " ".join(fields[1:])
        for loc in geoloc_map.get(bs, ()):
            yield "{}\t{}\t{}".format(bs, rest, loc)


def project_metadata(fields, acc_index):
    """Project a merged metadata row to the all_asm_acc_metadata column order.

    Mirrors the awk projection ``$N,$4,$4"."$5,$1,$8,$9,$10,$11`` where N is the
    accession column (2 = RefSeq for GCF rows, 3 = Genbank for GCA rows):

        accession  species  species.strain  biosample  taxid  taxcheck
        excl_from_refseq  geo_loc

    Short rows are padded with empty strings, matching awk's behaviour of
    printing an empty string for a field past NF.
    """
    def f(i):
        # 1-based, awk-style: out of range -> empty string.
        return fields[i - 1] if 0 < i <= len(fields) else ""

    return " ".join([
        f(acc_index),
        f(4),
        "{}.{}".format(f(4), f(5)),
        f(1),
        f(8),
        f(9),
        f(10),
        f(11),
    ])


def filter_by_acc(rows, acc_set, acc_index):
    """Yield projected lines for rows whose accession column is in acc_set.

    ``acc_index`` is the 1-based accession column: 2 (RefSeq) or 3 (Genbank).
    """
    for _raw, fields in rows:
        if not fields:
            continue
        acc = fields[acc_index - 1] if acc_index <= len(fields) else ""
        if acc in acc_set:
            yield project_metadata(fields, acc_index)


def exclude_downloaded(rows, downloaded):
    """Yield raw lines whose RefSeq (col 2) and Genbank (col 3) are both unseen.

    Replaces `grep -vf assemblies_datasets_uniq.names`: keep only assemblies that
    `datasets` did NOT already fetch. Lines pass through verbatim because the
    caller re-parses them with its own field layout.
    """
    for raw, fields in rows:
        refseq = fields[1] if len(fields) > 1 else ""
        genbank = fields[2] if len(fields) > 2 else ""
        if refseq in downloaded or genbank in downloaded:
            continue
        yield raw


# --------------------------------------------------------------------------- #
# Subcommands                                                                 #
# --------------------------------------------------------------------------- #

def cmd_merge_geoloc(args):
    geoloc_map = read_geoloc_map(args.geoloc)
    n_in = n_out = 0
    with open(args.out, "w") as out:
        for _raw, fields in iter_rows(args.meta):
            if not fields:
                continue
            n_in += 1
            bs = fields[0]
            rest = " ".join(fields[1:])
            locs = geoloc_map.get(bs)
            if not locs:
                continue
            for loc in locs:
                out.write("{}\t{}\t{}\n".format(bs, rest, loc))
                n_out += 1

    sys.stderr.write(
        "merge-geoloc: {} metadata rows, {} biosamples with geolocation, "
        "{} joined rows\n".format(n_in, len(geoloc_map), n_out)
    )
    dropped = n_in - n_out
    if n_in and not geoloc_map:
        sys.stderr.write(
            "WARNING: no biosample geolocation data; "
            "{} is empty (geolocation is decorative metadata).\n".format(args.out)
        )
    elif dropped > 0:
        sys.stderr.write(
            "NOTE: {} metadata row(s) had no matching biosample geolocation and "
            "were dropped (same as the previous bash join).\n".format(dropped)
        )
    return 0


def cmd_filter_by_acc(args):
    gcf, gca = split_accessions(read_acc_set(args.acc_list))

    # Preserve the intermediate accession lists the bash version wrote, so
    #   existing debugging habits (and any manual inspection) still work.
    if args.out_gcf:
        with open(args.out_gcf, "w") as fh:
            for acc in sorted(gcf):
                fh.write(acc + "\n")
    if args.out_gca:
        with open(args.out_gca, "w") as fh:
            for acc in sorted(gca):
                fh.write(acc + "\n")

    # Two passes so memory stays O(accessions), not O(rows): the bash version
    #   appended all GCF matches first, then all GCA matches, and downstream
    #   readers are order-insensitive but we keep the same layout.
    n_gcf = n_gca = 0
    with open(args.out, "w") as out:
        for line in filter_by_acc(iter_rows(args.meta), gcf, acc_index=2):
            out.write(line + "\n")
            n_gcf += 1
        for line in filter_by_acc(iter_rows(args.meta), gca, acc_index=3):
            out.write(line + "\n")
            n_gca += 1

    sys.stderr.write(
        "filter-by-acc: {} RefSeq + {} GenBank accessions requested -> "
        "{} + {} metadata rows written\n".format(len(gcf), len(gca), n_gcf, n_gca)
    )
    return 0


def cmd_exclude_downloaded(args):
    downloaded = read_acc_set(args.downloaded)
    n_in = n_out = 0
    with open(args.out, "w") as out:
        for raw, fields in iter_rows(args.meta):
            n_in += 1
            refseq = fields[1] if len(fields) > 1 else ""
            genbank = fields[2] if len(fields) > 2 else ""
            if refseq in downloaded or genbank in downloaded:
                continue
            out.write(raw + "\n")
            n_out += 1

    sys.stderr.write(
        "exclude-downloaded: {} metadata rows, {} already downloaded, "
        "{} candidates for supplementary download\n".format(
            n_in, len(downloaded), n_out
        )
    )
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description="Fast metadata joins/filters for gather_filter_asms.sh.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command")
    # Python 3.6 compatibility: `required=` on add_subparsers is 3.7+.
    sub.required = True

    p = sub.add_parser(
        "merge-geoloc",
        help="Join assembly metadata with biosample geolocation (hash join).",
    )
    p.add_argument("--meta", required=True,
                   help="All-<taxon>.assembly.BS_to_meta")
    p.add_argument("--geoloc", required=True,
                   help="All-<taxon>.biosample.BS_to_Geoloc")
    p.add_argument("--out", required=True,
                   help="Output All-<taxon>.BS_to_all_meta")
    p.set_defaults(func=cmd_merge_geoloc)

    p = sub.add_parser(
        "filter-by-acc",
        help="Keep metadata rows for downloaded accessions and project columns.",
    )
    p.add_argument("--meta", required=True,
                   help="All-<taxon>.BS_to_all_meta")
    p.add_argument("--acc-list", required=True,
                   help="all_asm_acc (mixed GCF_/GCA_ accessions)")
    p.add_argument("--out", required=True, help="Output all_asm_acc_metadata")
    p.add_argument("--out-gcf", default=None,
                   help="Optional: write the RefSeq accession subset here")
    p.add_argument("--out-gca", default=None,
                   help="Optional: write the GenBank accession subset here")
    p.set_defaults(func=cmd_filter_by_acc)

    p = sub.add_parser(
        "exclude-downloaded",
        help="Drop metadata rows whose assembly `datasets` already fetched.",
    )
    p.add_argument("--meta", required=True,
                   help="All-<taxon>.assembly.BS_to_meta")
    p.add_argument("--downloaded", required=True,
                   help="assemblies_datasets_uniq.names")
    p.add_argument("--out", required=True,
                   help="Output assemblies_not_in_datasets.BS_to_meta")
    p.set_defaults(func=cmd_exclude_downloaded)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
