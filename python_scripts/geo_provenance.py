#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geo_provenance.py

Flags cases where an assembly's "sampling" geolocation (all_asm_acc_metadata
column 8, from BioSample geo_loc_name) might actually be the location of the
lab that sequenced or submitted the genome, rather than where the sample was
collected. NCBI's schema keeps these concepts separate (geo_loc_name is
sample-collection metadata; sequenced_by / Assembly SubmitterOrganization /
BioSample Owner describe who did the work), but nothing stops a submitter
from typing their own institution's country into geo_loc_name -- which this
module cannot see directly, but CAN flag by comparing the two.

This is deliberately FLAG-ONLY: it never drops or relabels a genome. A
"differ" or "unresolved" flag just means a human should look closer.

Column contract (all_asm_acc_metadata, 1-based, whitespace-split):
    1 accession      5 taxid           9  collected_by
    2 species        6 taxcheck        10 sequenced_by
    3 species.strain 7 excl_refseq     11 biosample_owner
    4 biosample       8 geo_country     12 asm_submitter
Columns 9-12 are produced by parse_asm_metadata.py's project_metadata(); this
module reads a file with (at least) those 12 columns and appends:
    13 submitter_country   14 geo_vs_submitter

Usage:
    geo_provenance.py annotate --meta all_asm_acc_metadata --out F \\
        --report geo_provenance_report.tsv
"""

import argparse
import json
import os
import sys

MISSING = "NA"
UNKNOWN = "Unknown"

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
_DEFAULT_OVERRIDES = os.path.join(_DATA_DIR, "org_country_overrides.tsv")

# Values that mean "no usable text", mirroring parse_asm_metadata's MISSING/
# "Unknown" sentinels plus a plain empty string.
_EMPTY_VALUES = {"", MISSING.lower(), UNKNOWN.lower(), "none", "null", "-"}

# Country names that are too ambiguous to match from free text (also a very
# common non-country word/name, e.g. a US state or a person's name) --
# "Georgia" collides with the US state (and e.g. "University of Georgia").
_AMBIGUOUS_COUNTRIES = {"georgia"}

# candidate_lower -> tuple of preceding words that disqualify the match.
# "New Mexico" is a US state, not the country Mexico; "New Jersey"/"New York"
# would be similar but neither is a country name in our index.
_PRECEDING_WORD_GUARDS = {
    "mexico": ("new",),
}


# --------------------------------------------------------------------------- #
# Country index: NE 110m country names + INSDC/NCBI alias table             #
# --------------------------------------------------------------------------- #

class CountryIndex(object):
    """Lower-cased name/alias -> canonical NE country name.

    `names` is the full list of canonical country names (for text scanning,
    longest-first so multi-word names win over short substrings).
    """

    def __init__(self, name_lookup, alias_lookup, names):
        self.name_lookup = name_lookup
        self.alias_lookup = alias_lookup
        self.names = names

    def canonicalize(self, raw):
        """Return the canonical NE country name for a geo_loc_name-style
        string (already country-only, underscore-joined), or None."""
        if raw is None:
            return None
        text = raw.replace("_", " ").strip()
        if not text or text.lower() in _EMPTY_VALUES:
            return None
        low = text.lower()
        if low in self.name_lookup:
            return self.name_lookup[low]
        if low in self.alias_lookup:
            return self.alias_lookup[low]
        return None

    def match_candidates(self):
        """Candidate (lower, canonical) pairs for free-text scanning, longest
        surface form first so e.g. "United States of America" is tried
        before any shorter substring it contains."""
        seen = {}
        for low, name in self.name_lookup.items():
            if low in _AMBIGUOUS_COUNTRIES:
                continue
            seen[low] = name
        for low, name in self.alias_lookup.items():
            if low in _AMBIGUOUS_COUNTRIES:
                continue
            seen.setdefault(low, name)
        return sorted(seen.items(), key=lambda kv: len(kv[0]), reverse=True)


def _load_ne_country_names():
    path = os.path.join(_DATA_DIR, "ne_110m_countries.json")
    with open(path) as fh:
        data = json.load(fh)
    return [c["name"] for c in data["countries"]]


def _load_raw_aliases():
    path = os.path.join(_DATA_DIR, "country_aliases.tsv")
    aliases = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 3:
                continue
            insdc_name, kind, value = parts
            aliases[insdc_name.lower()] = (kind, value)
    return aliases


def load_country_index():
    """Build a CountryIndex from the bundled Natural Earth country list and
    the INSDC/NCBI alias table (same data files annotate_tree_by_country.py
    uses for map rendering, loaded independently here to avoid a matplotlib
    import). Only "country"-kind aliases are usable for country matching;
    "latlon"-kind entries (seas, oceans, regions) have no canonical country
    and are skipped.
    """
    names = _load_ne_country_names()
    name_lookup = {n.lower(): n for n in names}
    raw_aliases = _load_raw_aliases()

    alias_lookup = {}
    # Fixed-point resolution: an alias can point at another alias that is
    # itself a country-kind entry (not a direct NE name), so a single pass
    # over the dict can miss a chain depending on iteration order -- same
    # issue annotate_tree_by_country.load_country_coords() guards against.
    pending = {
        alias: value for alias, (kind, value) in raw_aliases.items()
        if kind == "country"
    }
    changed = True
    while changed and pending:
        changed = False
        for alias, value in list(pending.items()):
            vl = value.lower()
            if vl in name_lookup:
                alias_lookup[alias] = name_lookup[vl]
                del pending[alias]
                changed = True
            elif vl in alias_lookup:
                alias_lookup[alias] = alias_lookup[vl]
                del pending[alias]
                changed = True

    return CountryIndex(name_lookup, alias_lookup, names)


def canonical_country(raw, index):
    """Canonicalize a country-only string (e.g. all_asm_acc_metadata col 8,
    or a resolved submitter country) so equal countries compare equal
    regardless of spelling ("USA" / "United_States" / "United States of
    America" -> the same NE name)."""
    return index.canonicalize(raw)


# --------------------------------------------------------------------------- #
# Free-text country extraction (for submitter/sequencer/owner strings)       #
# --------------------------------------------------------------------------- #

def country_from_text(text, index):
    """Scan free text (an institution/submitter name) for a country name or
    known alias, preferring the LAST occurrence (affiliations conventionally
    end with the country, e.g. "Dept of X, University of Y, Country").

    Guards against well-known false positives (_AMBIGUOUS_COUNTRIES,
    _PRECEDING_WORD_GUARDS) rather than matching every substring blindly.
    Returns the canonical NE country name, or None if nothing matched.
    """
    if not text:
        return None
    low = text.replace("_", " ").lower()
    if low.strip() in _EMPTY_VALUES:
        return None

    import re

    best = None  # (start_index, canonical_name)
    for candidate, canonical in index.match_candidates():
        pattern = r"(?<!\w)" + re.escape(candidate) + r"(?!\w)"
        for m in re.finditer(pattern, low):
            guard_words = _PRECEDING_WORD_GUARDS.get(candidate)
            if guard_words:
                prefix = low[: m.start()].rstrip()
                prev_word = prefix.split()[-1] if prefix.split() else ""
                if prev_word in guard_words:
                    continue
            if best is None or m.start() >= best[0]:
                best = (m.start(), canonical)

    return best[1] if best else None


def load_overrides(path=None):
    """Read org_country_overrides.tsv: substring (case-insensitive) ->
    country name. Returns a list of (lower_substring, country) sorted
    longest-substring-first so the most specific match wins."""
    path = path or _DEFAULT_OVERRIDES
    overrides = []
    if not path or not os.path.exists(path):
        return overrides
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            substring, country = parts
            overrides.append((substring.lower(), country))
    overrides.sort(key=lambda kv: len(kv[0]), reverse=True)
    return overrides


def match_overrides(text, overrides):
    """First (longest-substring) override whose key appears in text."""
    if not text:
        return None
    low = text.replace("_", " ").lower()
    for substring, country in overrides:
        if substring in low:
            return country
    return None


# --------------------------------------------------------------------------- #
# Submitter country resolution + the geo-vs-submitter flag                   #
# --------------------------------------------------------------------------- #

def is_usable_text(value):
    if not value:
        return False
    return value.strip().lower() not in _EMPTY_VALUES


def resolve_submitter_country(sequenced_by, asm_submitter, biosample_owner,
                               index, overrides):
    """Try, in order, sequenced_by -> asm_submitter -> biosample_owner: the
    override table first (for submitters whose name mentions no place),
    then free-text country scanning. Returns (country, source_field) or
    (None, None) if nothing resolved.
    """
    for source_name, value in (
        ("sequenced_by", sequenced_by),
        ("asm_submitter", asm_submitter),
        ("biosample_owner", biosample_owner),
    ):
        if not is_usable_text(value):
            continue
        text = value.replace("_", " ")
        country = match_overrides(text, overrides)
        if country:
            canon = canonical_country(country, index) or country
            return canon, source_name
        country = country_from_text(text, index)
        if country:
            return country, source_name
    return None, None


def geo_vs_submitter(geo_country, submitter_country, index):
    """Compare all_asm_acc_metadata's col-8 geo country against the resolved
    submitter country. Flag values:
      match      - same country (the leak signal this module exists to find)
      differ     - both resolved, and they disagree (no leak signal)
      no_geo     - no usable geo_loc_name to compare against
      unresolved - geo is usable but no submitter country could be resolved
    """
    geo_canon = canonical_country(geo_country, index)
    if geo_canon is None:
        return "no_geo"
    if submitter_country is None:
        return "unresolved"
    return "match" if geo_canon == submitter_country else "differ"


# --------------------------------------------------------------------------- #
# Row IO                                                                      #
# --------------------------------------------------------------------------- #

# project_metadata() in parse_asm_metadata.py appends these four columns
# (9-12) ahead of this module's two (13-14).
_META_FIELD_NAMES = (
    "accession", "species", "species_strain", "biosample", "taxid",
    "taxcheck", "excl_refseq", "geo_country", "collected_by", "sequenced_by",
    "biosample_owner", "asm_submitter",
)
N_META_FIELDS = len(_META_FIELD_NAMES)


def iter_rows(path):
    if not path or not os.path.exists(path):
        return
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            yield line.split()


def _pad(fields, n):
    fields = list(fields[:n])
    fields += [MISSING] * (n - len(fields))
    return fields


def annotate_rows(rows, index, overrides):
    """Yield (extended_fields, record) for each input row.

    extended_fields: the 12 input fields plus submitter_country and
    geo_vs_submitter (14 total). record: a dict of every named field, for
    the report writer and for the stderr summary.
    """
    for fields in rows:
        row = _pad(fields, N_META_FIELDS)
        rec = dict(zip(_META_FIELD_NAMES, row))
        sub_country, source = resolve_submitter_country(
            rec["sequenced_by"], rec["asm_submitter"], rec["biosample_owner"],
            index, overrides,
        )
        flag = geo_vs_submitter(rec["geo_country"], sub_country, index)
        # Canonical NE country names can contain spaces ("United States of
        # America"); the output row is whitespace-split downstream (same
        # contract as parse_asm_metadata.sanitize_field), so this column
        # cannot be allowed to contain one.
        rec["submitter_country"] = (sub_country or MISSING).replace(" ", "_")
        rec["submitter_source"] = source or MISSING
        rec["geo_vs_submitter"] = flag
        extended = row + [rec["submitter_country"], rec["geo_vs_submitter"]]
        yield extended, rec


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #

def cmd_annotate(args):
    index = load_country_index()
    overrides = load_overrides(args.overrides)

    n_total = 0
    flag_counts = {}
    # submitter_source value -> [n_match, n_total], for the "lab-default"
    # fingerprint summary: an org whose samples are overwhelmingly "match"
    # is a candidate for a systematic geo_loc_name-is-really-the-lab issue.
    org_counts = {}

    report_header = _META_FIELD_NAMES + (
        "submitter_country", "submitter_source", "geo_vs_submitter",
    )

    with open(args.out, "w") as out, open(args.report, "w") as report:
        report.write("\t".join(report_header) + "\n")
        for extended, rec in annotate_rows(iter_rows(args.meta), index, overrides):
            n_total += 1
            out.write(" ".join(extended) + "\n")
            report.write("\t".join(rec[name] for name in report_header) + "\n")

            flag_counts[rec["geo_vs_submitter"]] = (
                flag_counts.get(rec["geo_vs_submitter"], 0) + 1
            )

            if rec["submitter_source"] != MISSING:
                org_key = (rec["submitter_source"], rec[rec["submitter_source"]])
                counts = org_counts.setdefault(org_key, [0, 0])
                counts[1] += 1
                if rec["geo_vs_submitter"] == "match":
                    counts[0] += 1

    sys.stderr.write(
        "geo-provenance annotate: {} rows -> {} (14 columns)\n".format(
            n_total, args.out
        )
    )
    for flag in ("match", "differ", "no_geo", "unresolved"):
        if flag in flag_counts:
            sys.stderr.write("  {}: {}\n".format(flag, flag_counts[flag]))

    lab_default = [
        (org_key, counts) for org_key, counts in org_counts.items()
        if counts[1] >= 5 and counts[0] / float(counts[1]) >= 0.9
    ]
    if lab_default:
        lab_default.sort(key=lambda kv: kv[1][1], reverse=True)
        sys.stderr.write(
            "  {} submitter(s) have >=90% geo==submitter-country match "
            "(possible lab-default geolocation):\n".format(len(lab_default))
        )
        for (source, org), (n_match, n_tot) in lab_default[:10]:
            sys.stderr.write(
                "    [{}] {}: {}/{} match\n".format(source, org, n_match, n_tot)
            )

    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description="Flag NCBI geo_loc_name values that may actually be "
                     "the sequencing/submitting lab's location.",
    )
    sub = parser.add_subparsers(dest="command")
    sub.required = True

    p = sub.add_parser(
        "annotate",
        help="Append submitter_country/geo_vs_submitter to "
             "all_asm_acc_metadata and write a per-genome report.",
    )
    p.add_argument("--meta", required=True,
                    help="all_asm_acc_metadata (>=12 columns)")
    p.add_argument("--out", required=True,
                    help="Output file (14 columns)")
    p.add_argument("--report", required=True,
                    help="Output geo_provenance_report.tsv")
    p.add_argument("--overrides", default=None,
                    help="Override org->country TSV "
                         "(default: data/org_country_overrides.tsv)")
    p.set_defaults(func=cmd_annotate)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
