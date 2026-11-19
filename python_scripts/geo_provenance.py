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
    1 accession      5 taxid           9  collected_by    13 lat_lon
    2 species        6 taxcheck        10 sequenced_by
    3 species.strain 7 excl_refseq     11 biosample_owner
    4 biosample       8 geo_country     12 asm_submitter
Columns 9-13 are produced by parse_asm_metadata.py's project_metadata(); this
module reads a file with (at least) those 13 columns and appends:
    14 submitter_country   15 geo_vs_submitter
    16 latlon_country      17 latlon_vs_geo
    18 submitter_km        19 latlon_vs_submitter

latlon_* columns are NA/no_country/unresolved whenever column 13 (lat_lon)
has no usable value -- most genomes have no BioSample lat_lon, so these are
best-effort, not universal. submitter_km/latlon_vs_submitter additionally
require --ror (opt-in network lookup of the submitter org's own coordinates);
without it they stay "NA"/"unresolved".

Usage:
    geo_provenance.py annotate --meta all_asm_acc_metadata --out F \\
        --report geo_provenance_report.tsv [--ror --ror-cache ror_cache.json]
        [--at-submitter-km 25]
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

# A point that misses every polygon (common right at a simplified 110m
# coastline) is assigned to the nearest boundary vertex if it is within this
# distance -- beyond it we call it open ocean / no_country rather than
# guessing.
_COASTAL_FALLBACK_KM = 50.0

# Default threshold for latlon_vs_submitter's "at_submitter" flag: a sample's
# coordinates landing within this many km of the submitter org's own
# coordinates (from ROR) is the strongest available leak signal -- the
# "sample" location IS the lab.
DEFAULT_AT_SUBMITTER_KM = 25.0

_ROR_BASE_URL = "https://api.ror.org/v2/organizations"


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
# lat_lon parsing, point-in-polygon country resolution, and distance          #
# --------------------------------------------------------------------------- #

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance between two lat/lon points, in km."""
    import math

    r = 6371.0088  # IUGG mean Earth radius
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = (
        math.sin(dphi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    )
    return 2.0 * r * math.asin(min(1.0, math.sqrt(a)))


def parse_lat_lon(raw):
    """Parse a BioSample lat_lon value into (lat, lon) floats, or None.

    ``raw`` is already sanitize_field()-clean (underscore-joined, no raw
    whitespace), so the first step undoes that. Accepts INSDC format
    ("38.98 N 77.11 W") and plain decimal ("38.98, -77.11" / "38.98 -77.11").
    Out-of-range values (|lat|>90 or |lon|>180) and the missing/unknown
    sentinels are rejected, returning None.
    """
    if not raw:
        return None
    text = raw.replace("_", " ").strip()
    if not text or text.lower() in _EMPTY_VALUES:
        return None

    import re

    m = re.match(
        r"^(-?\d+(?:\.\d+)?)\s*([NSns])\s*,?\s+(-?\d+(?:\.\d+)?)\s*([EWew])$",
        text,
    )
    if m:
        lat = float(m.group(1))
        if m.group(2).upper() == "S":
            lat = -lat
        lon = float(m.group(3))
        if m.group(4).upper() == "W":
            lon = -lon
    else:
        parts = [p for p in re.split(r"[,\s]+", text) if p]
        if len(parts) != 2:
            return None
        try:
            lat, lon = float(parts[0]), float(parts[1])
        except ValueError:
            return None

    if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
        return None
    return (lat, lon)


_COUNTRY_POLYGONS_CACHE = None


def load_country_polygons():
    """Load (and cache) the bundled NE 110m country polygons.

    Each entry: {"name": ..., "polygons": [polygon, ...]}, where a polygon
    is a list of rings and a ring is a list of [lon, lat] points -- the raw
    shape of ne_110m_countries.json, used directly rather than through
    CountryIndex (which only keeps names, not geometry).
    """
    global _COUNTRY_POLYGONS_CACHE
    if _COUNTRY_POLYGONS_CACHE is None:
        path = os.path.join(_DATA_DIR, "ne_110m_countries.json")
        with open(path) as fh:
            data = json.load(fh)
        _COUNTRY_POLYGONS_CACHE = data["countries"]
    return _COUNTRY_POLYGONS_CACHE


def _point_in_ring(lat, lon, ring):
    inside = False
    n = len(ring)
    j = n - 1
    for i in range(n):
        xi, yi = ring[i]
        xj, yj = ring[j]
        if (yi > lat) != (yj > lat):
            x_at_lat = (xj - xi) * (lat - yi) / (yj - yi) + xi
            if lon < x_at_lat:
                inside = not inside
        j = i
    return inside


def _point_in_polygon(lat, lon, polygon):
    """polygon: a list of rings. Rings are XOR-ed together, which is the
    standard even-odd way of handling holes regardless of winding order."""
    inside = False
    for ring in polygon:
        if _point_in_ring(lat, lon, ring):
            inside = not inside
    return inside


def country_at(lat, lon, countries=None):
    """Resolve (lat, lon) to a canonical NE country name via point-in-polygon.

    Falls back to the nearest boundary vertex across all countries if the
    point misses every polygon (common for a coastal point against the
    simplified 110m outline) and that vertex is within
    _COASTAL_FALLBACK_KM. Returns None for open ocean / no match (the
    caller maps that to the "no_country" flag).
    """
    countries = countries if countries is not None else load_country_polygons()

    for country in countries:
        for polygon in country["polygons"]:
            if _point_in_polygon(lat, lon, polygon):
                return country["name"]

    best = None  # (distance_km, name)
    for country in countries:
        for polygon in country["polygons"]:
            for ring in polygon:
                for plon, plat in ring:
                    d = haversine_km(lat, lon, plat, plon)
                    if best is None or d < best[0]:
                        best = (d, country["name"])
    if best is not None and best[0] <= _COASTAL_FALLBACK_KM:
        return best[1]
    return None


# --------------------------------------------------------------------------- #
# ROR lookup (opt-in, network): submitter org -> country/lat/lng             #
# --------------------------------------------------------------------------- #

def _ror_http_fetch(org, timeout=10):
    """Default `fetch` for ror_lookup(): a plain urllib GET. Tests inject
    their own `fetch` so this never actually runs under pytest."""
    import urllib.parse
    import urllib.request

    url = _ROR_BASE_URL + "?affiliation=" + urllib.parse.quote(org)
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def load_ror_cache(path):
    if path and os.path.exists(path):
        try:
            with open(path) as fh:
                return json.load(fh)
        except ValueError:
            return {}
    return {}


def save_ror_cache(path, cache):
    if not path:
        return
    with open(path, "w") as fh:
        json.dump(cache, fh, indent=2, sort_keys=True)


def ror_lookup(org, cache, fetch=None):
    """Resolve a submitter org string to its own (country, lat, lng) via the
    ROR API, accepting only a `chosen: true` match (ROR's own best-match
    flag) with usable geonames_details coordinates.

    ``cache`` is a plain dict the caller persists as JSON (load_ror_cache /
    save_ror_cache); a hit avoids the network call entirely, so a submitter
    seen across thousands of genomes costs one real lookup. ``fetch`` is
    injectable for tests -- defaults to a real HTTP GET. Any fetch error
    (timeout, bad JSON, DNS, ...) is swallowed: this check is a nice-to-have,
    never worth failing a pipeline run over. Returns (country, lat, lng), or
    (None, None, None) if nothing resolved.
    """
    if not is_usable_text(org):
        return (None, None, None)
    if org in cache:
        cached = cache[org]
        return tuple(cached) if cached else (None, None, None)

    fetch = fetch or _ror_http_fetch
    try:
        payload = fetch(org)
    except Exception:
        payload = None

    result = (None, None, None)
    if payload:
        for item in payload.get("items", []):
            if not item.get("chosen"):
                continue
            locations = (item.get("organization") or {}).get("locations") or []
            if not locations:
                break
            geo = locations[0].get("geonames_details") or {}
            country, lat, lng = (
                geo.get("country_name"), geo.get("lat"), geo.get("lng"),
            )
            if country is not None and lat is not None and lng is not None:
                result = (country, float(lat), float(lng))
            break  # only the chosen match is authoritative

    cache[org] = list(result)
    return result


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


def resolve_latlon_country(lat_lon_raw, countries=None):
    """Parse all_asm_acc_metadata col 13 and resolve it to a country via
    point-in-polygon. Returns (lat, lon, country) -- lat/lon are None if the
    value is missing/unparseable, country is None if a value parsed but
    resolved to no country (ocean / too far from the coastal fallback)."""
    parsed = parse_lat_lon(lat_lon_raw)
    if parsed is None:
        return None, None, None
    lat, lon = parsed
    return lat, lon, country_at(lat, lon, countries)


def latlon_vs_geo(lat_lon_raw, geo_country, index, countries=None):
    """Compare the BioSample lat_lon's resolved country against col-8's geo
    country. Flag values:
      NA         - no usable lat_lon to check
      no_country - lat_lon parsed but didn't resolve to any country
      no_geo     - lat_lon resolved to a country, but geo_loc_name is unusable
      match/differ - both resolved; same country or not
    """
    lat, lon, latlon_country = resolve_latlon_country(lat_lon_raw, countries)
    if lat is None:
        return "NA"
    if latlon_country is None:
        return "no_country"
    geo_canon = canonical_country(geo_country, index)
    if geo_canon is None:
        return "no_geo"
    return "match" if geo_canon == latlon_country else "differ"


def latlon_vs_submitter(lat, lon, latlon_country, submitter_lat, submitter_lon,
                         submitter_country, at_submitter_km):
    """Compare the BioSample lat_lon against the submitter org's own
    coordinates (from an opt-in ROR lookup). Returns (submitter_km, flag):
      submitter_km - great-circle distance in km, or None if unresolved
      at_submitter - within at_submitter_km of the submitter org (the
                     strongest leak signal: the "sample" IS the lab)
      same_country - farther than that, but the same country
      differ       - resolved and in different countries
      unresolved   - lat_lon and/or the submitter's own coordinates (ROR)
                     are unavailable
    """
    if lat is None or lon is None or submitter_lat is None or submitter_lon is None:
        return None, "unresolved"
    dist = haversine_km(lat, lon, submitter_lat, submitter_lon)
    if dist <= at_submitter_km:
        return dist, "at_submitter"
    if latlon_country and submitter_country and latlon_country == submitter_country:
        return dist, "same_country"
    return dist, "differ"


def resolve_submitter_ror(sequenced_by, asm_submitter, biosample_owner,
                           cache, fetch=None):
    """Try, in order, asm_submitter -> biosample_owner -> sequenced_by: ROR's
    affiliation matcher works best against a formal institution name, which
    asm_submitter/biosample_owner are more likely to be than sequenced_by
    (often a generic surveillance network like "PulseNet"). Returns
    (country, lat, lng, source_field) or (None, None, None, None).
    """
    for source_name, value in (
        ("asm_submitter", asm_submitter),
        ("biosample_owner", biosample_owner),
        ("sequenced_by", sequenced_by),
    ):
        if not is_usable_text(value):
            continue
        country, lat, lng = ror_lookup(value.replace("_", " "), cache, fetch)
        if country is not None:
            return country, lat, lng, source_name
    return None, None, None, None


# --------------------------------------------------------------------------- #
# Row IO                                                                      #
# --------------------------------------------------------------------------- #

# project_metadata() in parse_asm_metadata.py appends these five columns
# (9-13) ahead of this module's six (14-19).
_META_FIELD_NAMES = (
    "accession", "species", "species_strain", "biosample", "taxid",
    "taxcheck", "excl_refseq", "geo_country", "collected_by", "sequenced_by",
    "biosample_owner", "asm_submitter", "lat_lon",
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


def annotate_rows(rows, index, overrides, ror_enabled=False, ror_cache=None,
                   ror_fetch=None, at_submitter_km=DEFAULT_AT_SUBMITTER_KM,
                   countries=None):
    """Yield (extended_fields, record) for each input row.

    extended_fields: the 13 input fields plus submitter_country,
    geo_vs_submitter, latlon_country, latlon_vs_geo, submitter_km, and
    latlon_vs_submitter (19 total). record: a dict of every named field, for
    the report writer and for the stderr summary.

    ror_enabled gates the (opt-in, network) ROR lookup that resolves
    submitter_km/latlon_vs_submitter; without it those two columns are
    always NA/"unresolved". ror_cache is a plain dict the caller persists
    (load_ror_cache/save_ror_cache) -- repeated submitters across many rows
    cost at most one real lookup each.
    """
    ror_cache = ror_cache if ror_cache is not None else {}
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

        lat, lon, latlon_country = resolve_latlon_country(
            rec["lat_lon"], countries)
        rec["latlon_country"] = (latlon_country or MISSING).replace(" ", "_")
        rec["latlon_vs_geo"] = latlon_vs_geo(
            rec["lat_lon"], rec["geo_country"], index, countries)

        ror_country = ror_lat = ror_lng = None
        if ror_enabled:
            ror_country, ror_lat, ror_lng, _ror_source = resolve_submitter_ror(
                rec["sequenced_by"], rec["asm_submitter"],
                rec["biosample_owner"], ror_cache, ror_fetch,
            )
        dist, sub_flag = latlon_vs_submitter(
            lat, lon, latlon_country, ror_lat, ror_lng,
            (canonical_country(ror_country, index) if ror_country else None),
            at_submitter_km,
        )
        rec["submitter_km"] = MISSING if dist is None else "{:.1f}".format(dist)
        rec["latlon_vs_submitter"] = sub_flag

        extended = row + [
            rec["submitter_country"], rec["geo_vs_submitter"],
            rec["latlon_country"], rec["latlon_vs_geo"],
            rec["submitter_km"], rec["latlon_vs_submitter"],
        ]
        yield extended, rec


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #

def cmd_annotate(args):
    index = load_country_index()
    overrides = load_overrides(args.overrides)
    ror_cache = load_ror_cache(args.ror_cache) if args.ror else {}

    n_total = 0
    flag_counts = {}
    latlon_flag_counts = {}
    submitter_flag_counts = {}
    # submitter_source value -> [n_match, n_total], for the "lab-default"
    # fingerprint summary: an org whose samples are overwhelmingly "match"
    # is a candidate for a systematic geo_loc_name-is-really-the-lab issue.
    org_counts = {}

    report_header = _META_FIELD_NAMES + (
        "submitter_country", "submitter_source", "geo_vs_submitter",
        "latlon_country", "latlon_vs_geo", "submitter_km",
        "latlon_vs_submitter",
    )

    with open(args.out, "w") as out, open(args.report, "w") as report:
        report.write("\t".join(report_header) + "\n")
        for extended, rec in annotate_rows(
            iter_rows(args.meta), index, overrides,
            ror_enabled=args.ror, ror_cache=ror_cache,
            at_submitter_km=args.at_submitter_km,
        ):
            n_total += 1
            out.write(" ".join(extended) + "\n")
            report.write("\t".join(rec[name] for name in report_header) + "\n")

            flag_counts[rec["geo_vs_submitter"]] = (
                flag_counts.get(rec["geo_vs_submitter"], 0) + 1
            )
            latlon_flag_counts[rec["latlon_vs_geo"]] = (
                latlon_flag_counts.get(rec["latlon_vs_geo"], 0) + 1
            )
            submitter_flag_counts[rec["latlon_vs_submitter"]] = (
                submitter_flag_counts.get(rec["latlon_vs_submitter"], 0) + 1
            )

            if rec["submitter_source"] != MISSING:
                org_key = (rec["submitter_source"], rec[rec["submitter_source"]])
                counts = org_counts.setdefault(org_key, [0, 0])
                counts[1] += 1
                if rec["geo_vs_submitter"] == "match":
                    counts[0] += 1

    if args.ror:
        save_ror_cache(args.ror_cache, ror_cache)

    sys.stderr.write(
        "geo-provenance annotate: {} rows -> {} (19 columns)\n".format(
            n_total, args.out
        )
    )
    for flag in ("match", "differ", "no_geo", "unresolved"):
        if flag in flag_counts:
            sys.stderr.write("  geo_vs_submitter {}: {}\n".format(flag, flag_counts[flag]))
    for flag in ("match", "differ", "no_country", "no_geo", "NA"):
        if flag in latlon_flag_counts:
            sys.stderr.write("  latlon_vs_geo {}: {}\n".format(flag, latlon_flag_counts[flag]))
    if args.ror:
        for flag in ("at_submitter", "same_country", "differ", "unresolved"):
            if flag in submitter_flag_counts:
                sys.stderr.write(
                    "  latlon_vs_submitter {}: {}\n".format(
                        flag, submitter_flag_counts[flag])
                )

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
                    help="all_asm_acc_metadata (>=8 columns; 13 recommended)")
    p.add_argument("--out", required=True,
                    help="Output file (19 columns)")
    p.add_argument("--report", required=True,
                    help="Output geo_provenance_report.tsv")
    p.add_argument("--overrides", default=None,
                    help="Override org->country TSV "
                         "(default: data/org_country_overrides.tsv)")
    p.add_argument("--ror", action="store_true",
                    help="Opt-in: look up the submitter org's own "
                         "coordinates via the ROR API (network) to compute "
                         "submitter_km/latlon_vs_submitter")
    p.add_argument("--ror-cache", default=None,
                    help="JSON cache file for ROR lookups (persists across "
                         "runs; without --ror this is unused)")
    p.add_argument("--at-submitter-km", type=float,
                    default=DEFAULT_AT_SUBMITTER_KM,
                    help="Distance threshold (km) for the latlon_vs_submitter "
                         "'at_submitter' flag (default: {})".format(
                             DEFAULT_AT_SUBMITTER_KM))
    p.set_defaults(func=cmd_annotate)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
