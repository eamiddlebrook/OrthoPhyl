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

...and two subcommands that fix a CORRECTNESS bug rather than a speed one:

  parse-assembly-xml  <- get_asm_metadata()'s `xtract -element ...` pipeline
      xtract emits one whitespace-separated token per matched element, so a
      REPEATED element silently shifts every later column right. `Sub_value` is
      repeated whenever an assembly has more than one Infraspecie entry (e.g.
      both a culture-collection and a strain designation), which pushes
      taxonomy-check-status into the geolocation column -- the "geolocation is
      OK" bug. Measured on a real 75,155-row Pseudomonas run: 2,150 rows (2.9%)
      had extra columns. Parsing by ELEMENT NAME removes the whole bug class.

  parse-biosample-xml  <- get_biosample_GEOdata()'s xtract + 7-`sed` pipeline
      Same idea for the geolocation table, and it drops the fragile
      `sed 's/Missing.*$/Unknown/g'`-style normalisation in favour of explicit
      value matching.

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

Field layout of All-<taxon>.assembly.BS_to_meta (1-based). parse-assembly-xml
guarantees EXACTLY these 10 columns per row, one line per DocumentSummary:
    1 BioSampleAccn  2 RefSeq  3 Genbank  4 SpeciesName  5 Sub_value
    6 FtpPath_GenBank  7 FtpPath_RefSeq  8 Taxid  9 taxonomy-check-status
    10 ExclFromRefSeq
...and after merge-geoloc appends the geolocation: 11 geo_loc_name.

Usage:
    parse_asm_metadata.py parse-assembly-xml --xml F --out F
    parse_asm_metadata.py parse-biosample-xml --xml F --out F
    parse_asm_metadata.py merge-geoloc --meta F --geoloc F --out F
    parse_asm_metadata.py filter-by-acc --meta F --acc-list F --out F \\
        [--out-gcf F] [--out-gca F]
    parse_asm_metadata.py exclude-downloaded --meta F --downloaded F --out F
"""

import argparse
import os
import sys
import xml.etree.ElementTree as ET


# Value used for an absent field, matching xtract's `-def "NA"`.
MISSING = "NA"

# BS_to_meta must have exactly this many columns. Anything else means the row was
#   produced by the old xtract pipeline and a repeated element shifted the
#   columns; see the module docstring and the merge-geoloc recurrence guard.
#   Checking the width catches the shift directly, and is strictly better than
#   sniffing for taxcheck tokens in the geolocation column: "NA" is a legitimate
#   geolocation value, so a token-based guard would false-positive.
EXPECTED_META_COLUMNS = 10

# geo_loc_name values that mean "no real location". The bash pipeline matched
#   these with `sed 's/Missing.*$/Unknown/g'` and friends -- a prefix match on a
#   lowercased-or-not spelling. Compared case-insensitively against the value
#   prefix to reproduce that without the regex fragility.
UNKNOWN_GEOLOC_PREFIXES = (
    "missing",
    "unknown",
    "not collected",
    "not_collected",
    "not applicable",
    "not_applicable",
    "none",
    "not determined",
    "not_determined",
    "not provided",
    "not_provided",
)


# --------------------------------------------------------------------------- #
# XML parsing -- replaces the xtract pipelines (fixes the column-shift bug)    #
# --------------------------------------------------------------------------- #

def iter_document_summaries(path):
    """Stream <DocumentSummary> elements from an NCBI esummary XML file.

    This deliberately does NOT use a strict whole-document parser, because these
    files are NOT single XML documents. `esearch | esummary` appends one
    <DocumentSummarySet> root per fetched batch, so a real file is a CONCATENATION
    of many roots (the 75k-assembly Pseudomonas run has 76 of them, plus a repeated
    <?xml?>/<!DOCTYPE> preamble). ET.iterparse() stops at the end of the first root
    with "junk after document element" -- which silently yielded only the first
    1,000 of 75,155 records when this was first written.

    So: scan for <DocumentSummary>...</DocumentSummary> spans textually and parse
    each one on its own with ET.fromstring(). That is still real XML parsing per
    record (attributes, nesting and entities all handled by ET), it tolerates the
    multi-root layout and a truncated tail, and memory stays flat -- only one
    record is held at a time, which matters at 380-555 MB per file.

    A record that fails to parse is skipped with a warning rather than aborting
    the run; one bad record should not discard the other 75,154.
    """
    if not path or not os.path.exists(path):
        return

    open_tag = "<DocumentSummary>"
    open_tag_attr = "<DocumentSummary "   # e.g. <DocumentSummary uid="...">
    close_tag = "</DocumentSummary>"
    n_bad = 0

    with open(path, "r", errors="replace") as fh:
        buf = []
        inside = False
        for line in fh:
            if not inside and (open_tag in line or open_tag_attr in line):
                inside = True
                buf = []
            if inside:
                buf.append(line)
                if close_tag in line:
                    inside = False
                    chunk = "".join(buf)
                    buf = []
                    try:
                        yield ET.fromstring(chunk)
                    except ET.ParseError:
                        n_bad += 1

    if n_bad:
        sys.stderr.write(
            "WARNING: skipped {} malformed <DocumentSummary> record(s) in "
            "{}.\n".format(n_bad, path)
        )


def first_text(elem, tag, default=MISSING):
    """Return the text of the FIRST descendant <tag>, or default if absent/empty.

    Taking the first occurrence is what makes the output positionally stable:
    `Sub_value` repeats when an assembly has several Infraspecie entries, and
    xtract emitted one token per occurrence, shifting all later columns. One
    element name -> exactly one column, always.
    """
    for node in elem.iter(tag):
        text = (node.text or "").strip()
        if text:
            return text
        return default
    return default


def sanitize_field(value):
    """Collapse whitespace inside a field value to underscores.

    The bash pipeline ran `sed 's/ /_/g'` over the whole xtract output for exactly
    this reason: a species name like "Pseudomonas aeruginosa" would otherwise
    become two whitespace-separated columns. Tabs/newlines get the same treatment
    so a field can never split a row.
    """
    if value is None:
        return MISSING
    value = " ".join(str(value).split())
    if not value:
        return MISSING
    return value.replace(" ", "_")


def assembly_row(elem):
    """Extract the 10 BS_to_meta columns from one assembly <DocumentSummary>.

    Column order matches what get_asm_metadata()'s xtract emitted, so every
    downstream consumer keeps working -- but each column is now keyed by element
    name, so a repeated element cannot shift the row.
    """
    fields = [
        first_text(elem, "BioSampleAccn"),
        first_text(elem, "RefSeq"),
        first_text(elem, "Genbank"),
        first_text(elem, "SpeciesName"),
        first_text(elem, "Sub_value"),
        first_text(elem, "FtpPath_GenBank"),
        first_text(elem, "FtpPath_RefSeq"),
        first_text(elem, "Taxid"),
        first_text(elem, "taxonomy-check-status"),
        first_text(elem, "ExclFromRefSeq"),
    ]
    return [sanitize_field(f) for f in fields]


def normalize_geoloc(value):
    """Normalize a raw geo_loc_name to the form the old sed chain produced.

    - "India: Chhatrapati Sambhajinagar" -> "India"  (the bash split on ':' and
      kept the first field via awk; 75,609 values in the Pseudomonas run have a
      colon, so this matters).
    - missing/unknown/not-collected spellings -> "Unknown".
    - empty -> "NA" (the bash `awk '{if ($2 == "") print $1,"NA"}'`).
    - spaces -> underscores, so the value is one column.
    """
    if value is None:
        return MISSING
    value = " ".join(str(value).split())
    if not value:
        return MISSING

    low = value.lower()
    for prefix in UNKNOWN_GEOLOC_PREFIXES:
        if low.startswith(prefix):
            return "Unknown"

    # Keep only the country part, before the first colon.
    country = value.split(":", 1)[0].strip()
    if not country:
        return MISSING
    return country.replace(" ", "_")


def biosample_geoloc(elem):
    """Return (accession, normalized_geoloc) for a biosample <DocumentSummary>.

    Finds the Attribute whose harmonized_name is geo_loc_name, mirroring the
    xtract `-if Attribute@harmonized_name -equals geo_loc_name` selector. Returns
    ("", ...) if the record has no accession (skipped by the caller).
    """
    accession = first_text(elem, "Accession", default="")
    raw = None
    for attr in elem.iter("Attribute"):
        if attr.get("harmonized_name") == "geo_loc_name":
            raw = attr.text
            break
    return sanitize_field(accession) if accession else "", normalize_geoloc(raw)


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

def cmd_parse_assembly_xml(args):
    n_records = 0
    n_multi_subvalue = 0
    with open(args.out, "w") as out:
        for elem in iter_document_summaries(args.xml):
            # Count records that would have shifted columns under xtract, so the
            #   log shows how much the old pipeline was corrupting.
            if len(list(elem.iter("Sub_value"))) > 1:
                n_multi_subvalue += 1
            row = assembly_row(elem)
            out.write("\t".join(row) + "\n")
            n_records += 1

    sys.stderr.write(
        "parse-assembly-xml: {} records -> {} (10 columns each)\n".format(
            n_records, args.out
        )
    )
    if n_multi_subvalue:
        sys.stderr.write(
            "  {} record(s) had multiple <Sub_value> entries; under the old "
            "xtract pipeline these shifted every later column right (putting "
            "taxonomy-check-status where the geolocation belongs). Now keyed by "
            "element name, so columns stay aligned.\n".format(n_multi_subvalue)
        )
    if n_records == 0:
        sys.stderr.write(
            "WARNING: no DocumentSummary records parsed from {}; "
            "{} is empty.\n".format(args.xml, args.out)
        )
    return 0


def cmd_parse_biosample_xml(args):
    n_records = n_written = 0
    with open(args.out, "w") as out:
        for elem in iter_document_summaries(args.xml):
            n_records += 1
            accession, geoloc = biosample_geoloc(elem)
            if not accession:
                continue
            out.write("{}\t{}\n".format(accession, geoloc))
            n_written += 1

    sys.stderr.write(
        "parse-biosample-xml: {} records -> {} geolocation rows in {}\n".format(
            n_records, n_written, args.out
        )
    )
    if n_records == 0:
        sys.stderr.write(
            "WARNING: no DocumentSummary records parsed from {}; geolocation "
            "table is empty (geolocation is decorative metadata).\n".format(
                args.xml
            )
        )
    return 0


def cmd_merge_geoloc(args):
    geoloc_map = read_geoloc_map(args.geoloc)
    n_in = n_out = 0
    n_wrong_width = 0
    with open(args.out, "w") as out:
        for _raw, fields in iter_rows(args.meta):
            if not fields:
                continue
            n_in += 1
            # Recurrence guard for the column-shift bug: BS_to_meta must be
            #   exactly 10 columns, or the geolocation lands in the wrong place
            #   and taxonomy-check-status ("OK") gets read as a country.
            if len(fields) != EXPECTED_META_COLUMNS:
                n_wrong_width += 1
            bs = fields[0]
            rest = " ".join(fields[1:])
            locs = geoloc_map.get(bs)
            if not locs:
                continue
            for loc in locs:
                out.write("{}\t{}\t{}\n".format(bs, rest, loc))
                n_out += 1

    if n_wrong_width:
        sys.stderr.write(
            "WARNING: {} of {} metadata rows did not have the expected {} "
            "columns. Downstream column indices (taxonomy-check-status, "
            "geolocation) will be wrong for those rows. Regenerate {} with "
            "`parse_asm_metadata.py parse-assembly-xml` instead of the old "
            "xtract pipeline.\n".format(
                n_wrong_width, n_in, EXPECTED_META_COLUMNS, args.meta
            )
        )

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
        "parse-assembly-xml",
        help="Assembly esummary XML -> BS_to_meta (10 stable columns).",
    )
    p.add_argument("--xml", required=True,
                   help="All-<taxon>-info.assembly.xml")
    p.add_argument("--out", required=True,
                   help="Output All-<taxon>.assembly.BS_to_meta")
    p.set_defaults(func=cmd_parse_assembly_xml)

    p = sub.add_parser(
        "parse-biosample-xml",
        help="Biosample esummary XML -> accession<TAB>geolocation table.",
    )
    p.add_argument("--xml", required=True,
                   help="All-<taxon>-info.biosample.xml")
    p.add_argument("--out", required=True,
                   help="Output All-<taxon>.biosample.BS_to_Geoloc")
    p.set_defaults(func=cmd_parse_biosample_xml)

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
