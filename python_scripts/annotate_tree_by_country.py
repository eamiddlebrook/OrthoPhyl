#!/usr/bin/env python3
"""
Annotate OrthoPhyl species trees with the geographic origin of each assembly.

OrthoPhyl writes its final Newick trees to FINAL_SPECIES_TREES/ (or, via the
Python wrapper, to 03_results/trees/.../*.nwk). When the input assemblies were
gathered with utils/gather_filter_asms.sh, each tree leaf label is the assembly
accession (GCF_/GCA_...), which is column 1 of the gather output file
`all_asm_acc_metadata`. Column 8 of that same file is the country/geo-location
of origin (e.g. "USA:_Washington", underscore-joined, colon-separated region).

For genomes with no NCBI metadata (query genomes placed by ReLeaf, or genomes
ingested via --genome-dir), a user-supplied two-column TSV
(genome_id<TAB>location) can be passed with --geo-tsv; it takes precedence
over --metadata for any leaf present in both.

This script renders, per tree:
  - a country-annotated tree figure (colored tip strip + legend)
  - a world map of sample origins, with lines connecting the --top-k closest
    cross-country genome pairs by patristic (tree) distance. Line opacity is
    scaled so the single most distant pair in that top-K set is fully
    transparent and the closest pair is fully opaque.
  - <label>_geo_pairs.tsv and <label>_geo_countries.tsv data exports.

Tree rendering requires the OrthoPhyl conda env (ships ete3) and renders
headlessly via offscreen Qt -- no $DISPLAY needed. The map requires
matplotlib and numpy; if matplotlib is unavailable, pass --no-map to skip it
(the tree figure and data exports still work).

Example:
    python annotate_tree_by_country.py \\
        --trees-dir  /path/to/FINAL_SPECIES_TREES \\
        --metadata   /path/to/all_asm_acc_metadata \\
        --geo-tsv    /path/to/query_locations.tsv \\
        --out-dir    /path/to/03_results/figures \\
        --top-k 100
"""

import argparse
import colorsys
import glob
import hashlib
import json
import math
import os
import re
import sys

# ete3 needs a Qt backend; force offscreen so it works with no $DISPLAY.
# Must be set before ete3 (and therefore Qt) is imported.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

UNKNOWN = "Unknown"
# values in the metadata's country column that mean "no real location"
_UNKNOWN_VALUES = {"", "na", "unknown", "none", "null", "-"}

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

_ACCESSION_VERSION_RE = re.compile(r"^(GC[AF]_\d+)\.\d+$")
_FASTA_EXTS = (".fna", ".fasta", ".fa")


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate OrthoPhyl species trees with assembly geographic origin.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tree", default=None,
                   help="Path to a single Newick tree (.nwk/.tree/.treefile). "
                        "Mutually exclusive with --trees-dir.")
    p.add_argument("--trees-dir", default=None,
                   help="Directory containing Newick files (*.nwk/*.tree/*.treefile), "
                        "e.g. a FINAL_SPECIES_TREES directory. Mutually exclusive with --tree.")
    p.add_argument("--metadata", default=None,
                   help="Path to all_asm_acc_metadata from gather_filter_asms.sh. "
                        "At least one of --metadata/--geo-tsv is required.")
    p.add_argument("--geo-tsv", default=None,
                   help="Optional TSV of genome_id<TAB>location (header/# comments OK). "
                        "Overrides --metadata for any leaf present in both; the only "
                        "location source in --genome-dir mode (no NCBI metadata).")
    p.add_argument("--out-dir", default=None,
                   help="Directory for output figures "
                        "(default: <trees-dir>/geo_annotated, or alongside --tree).")
    p.add_argument("--format", default="png", choices=["png", "svg", "pdf"],
                   help="Output image format.")
    p.add_argument("--label", default=None,
                   help="Base name for output files when --tree is used "
                        "(default: the tree file's stem). Ignored for --trees-dir, "
                        "where each tree uses its own stem.")
    p.add_argument("--meta-acc-col", type=int, default=1,
                   help="1-based column in metadata holding the accession / leaf name.")
    p.add_argument("--meta-country-col", type=int, default=8,
                   help="1-based column in metadata holding the country of origin.")
    p.add_argument("--label-with-species", action="store_true",
                   help="Append the species.strain field (metadata col 3) to each tip label.")
    p.add_argument("--species-col", type=int, default=3,
                   help="1-based column used for --label-with-species.")
    p.add_argument("--dpi", type=int, default=300, help="Raster DPI for png output.")
    p.add_argument("--top-k", type=int, default=100,
                   help="Draw lines for only the K closest cross-country genome pairs "
                        "by patristic distance (0 disables pair lines).")
    p.add_argument("--no-map", action="store_true",
                   help="Skip the world-map figure (tree figure and TSV exports still run).")
    args = p.parse_args()
    if bool(args.tree) == bool(args.trees_dir):
        p.error("exactly one of --tree or --trees-dir is required")
    if not args.metadata and not args.geo_tsv:
        p.error("at least one of --metadata or --geo-tsv is required")
    return args


def clean_value(value):
    """Strip surrounding quotes/whitespace from a metadata token."""
    return value.strip().strip('"').strip("'")


def normalize_country(value):
    """Map a raw country token to a display bucket; missing -> Unknown.

    Handles the gather_filter_asms.sh encoding (underscores for spaces,
    "Country:_Region" for sub-national detail) as well as plain
    human-typed --geo-tsv values ("USA: Montana", "Viet Nam").
    """
    v = clean_value(value).replace("_", " ")
    v = v.split(":")[0].strip()
    if v.lower() in _UNKNOWN_VALUES:
        return UNKNOWN
    return v


def parse_metadata(path, acc_col, country_col, species_col, want_species):
    """Return (acc->country, acc->species) dicts. Columns are 1-based."""
    acc_i = acc_col - 1
    country_i = country_col - 1
    species_i = species_col - 1
    need = max(acc_i, country_i, species_i if want_species else 0)

    acc_country = {}
    acc_species = {}
    n_lines = 0
    n_dup = 0
    n_short = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) <= need:
                n_short += 1
                continue
            n_lines += 1
            acc = clean_value(fields[acc_i])
            if acc in acc_country:
                n_dup += 1
            acc_country[acc] = normalize_country(fields[country_i])
            if want_species:
                acc_species[acc] = clean_value(fields[species_i])

    eprint(f"[metadata] parsed {n_lines} records from {path}")
    if n_short:
        eprint(f"[metadata] WARNING: skipped {n_short} line(s) with too few columns")
    if n_dup:
        eprint(f"[metadata] WARNING: {n_dup} duplicate accession(s); last value kept")
    return acc_country, acc_species


def parse_geo_tsv(path):
    """Parse a user-supplied genome_id<TAB>location TSV into acc->country.

    A header row and '#'-prefixed comment lines are tolerated. Blank
    locations map to Unknown.
    """
    acc_country = {}
    n_lines = 0
    n_skipped = 0
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                n_skipped += 1
                continue
            genome_id, location = clean_value(parts[0]), parts[1]
            if genome_id.lower() in ("genome_id", "id", "accession", "genome"):
                # tolerate an optional header row
                continue
            acc_country[genome_id] = normalize_country(location)
            n_lines += 1

    eprint(f"[geo-tsv] parsed {n_lines} record(s) from {path}")
    if n_skipped:
        eprint(f"[geo-tsv] WARNING: skipped {n_skipped} malformed line(s)")
    return acc_country


def _candidate_keys(name):
    """Alternate forms of a leaf/accession name worth trying when matching
    against a metadata or geo-tsv key: the bare name, its stem without a
    FASTA extension, and its accession without an NCBI version suffix."""
    keys = [name]
    for ext in _FASTA_EXTS:
        if name.endswith(ext):
            keys.append(name[: -len(ext)])
    m = _ACCESSION_VERSION_RE.match(name)
    if m:
        keys.append(m.group(1))
    return keys


def _expand_keys(mapping):
    """Index a name->value mapping under every alternate key form too, so
    lookups succeed regardless of which side (leaf vs. metadata record)
    carries the extension/version suffix."""
    expanded = dict(mapping)
    for k, v in mapping.items():
        for ck in _candidate_keys(k):
            expanded.setdefault(ck, v)
    return expanded


def build_leaf_country_map(leaf_names, acc_country, geo_tsv_country):
    """Merge metadata + geo-tsv into a leaf_name -> country map.
    geo_tsv entries take precedence over metadata for the same leaf."""
    expanded_meta = _expand_keys(acc_country)
    expanded_tsv = _expand_keys(geo_tsv_country)

    leaf_country = {}
    missing = []
    for name in leaf_names:
        country = None
        for ck in _candidate_keys(name):
            if ck in expanded_meta:
                country = expanded_meta[ck]
                break
        for ck in _candidate_keys(name):
            if ck in expanded_tsv:
                country = expanded_tsv[ck]
                break
        if country is None:
            country = UNKNOWN
            missing.append(name)
        leaf_country[name] = country
    return leaf_country, missing


def build_color_map(countries):
    """Assign a stable, distinct color per country. Unknown -> light grey."""
    ordered = sorted(c for c in countries if c != UNKNOWN)
    color_map = {}

    # Try matplotlib's qualitative palettes first for nice distinct colors.
    palette = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.colors as mcolors
        for cmap_name in ("tab20", "tab20b", "tab20c"):
            cmap = matplotlib.colormaps[cmap_name]
            for i in range(cmap.N):
                palette.append(mcolors.to_hex(cmap(i)))
    except Exception as exc:  # matplotlib missing/old -> HSV fallback
        eprint(f"[color] matplotlib palette unavailable ({exc}); using HSV fallback")

    n = len(ordered)
    for idx, country in enumerate(ordered):
        if idx < len(palette):
            color_map[country] = palette[idx]
        else:
            # Evenly spaced hues so we never run out of distinct colors.
            hue = (idx / max(n, 1)) % 1.0
            r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.85)
            color_map[country] = "#{:02x}{:02x}{:02x}".format(
                int(r * 255), int(g * 255), int(b * 255))

    color_map[UNKNOWN] = "#d0d0d0"
    return color_map


def load_tree(path):
    """Load a Newick tree, tolerating internal support labels. Returns Tree or None."""
    from ete3 import Tree
    for fmt in (1, 0, 5):
        try:
            return Tree(path, format=fmt)
        except Exception:
            continue
    eprint(f"[tree] ERROR: could not parse {path}; skipping")
    return None


# --------------------------------------------------------------------------- #
# Patristic top-K cross-country pairs                                        #
# --------------------------------------------------------------------------- #

def top_k_cross_country_pairs(tree, leaf_country, k, row_chunk=2000):
    """Return the K closest (by patristic distance) genome pairs whose
    leaves have two different, known countries, as a list of
    (name_a, name_b, distance) sorted ascending by distance.

    For every internal node, the pairs whose MRCA is exactly that node are
    the outer product of leaves under each pair of its children (since any
    two leaves under the SAME child have a deeper, more specific MRCA).
    d(a, b) = depth(a) + depth(b) - 2*depth(mrca), and depth(mrca) is a
    single scalar per node, so each node's contribution is one vectorized
    numpy outer-sum. Candidate pairs are restricted to known, differing
    countries before ranking, and only the K best survive per chunk via
    argpartition, so memory never holds a full all-pairs matrix.
    """
    import heapq
    import itertools

    import numpy as np

    if k <= 0:
        return []

    leaves = tree.get_leaves()
    if len(leaves) < 2:
        return []
    names = [leaf.name for leaf in leaves]
    name_to_idx = {name: i for i, name in enumerate(names)}
    countries = np.array([leaf_country.get(n, UNKNOWN) for n in names], dtype=object)

    depth = {tree: 0.0}
    for node in tree.traverse("preorder"):
        if node.is_root():
            depth[node] = 0.0
        else:
            depth[node] = depth[node.up] + (node.dist or 0.0)
    leaf_depth = np.array([depth[leaf] for leaf in leaves], dtype=float)

    # max-heap of size <=k on distance, via negation (heap[0] = largest d kept)
    heap = []

    def _consider_block(idx_i, idx_j, node_depth):
        if idx_i.size == 0 or idx_j.size == 0:
            return
        for start in range(0, idx_i.size, row_chunk):
            sub_i = idx_i[start:start + row_chunk]
            d_block = (leaf_depth[sub_i][:, None] + leaf_depth[idx_j][None, :]
                       - 2.0 * node_depth)
            mismatch = countries[sub_i][:, None] != countries[idx_j][None, :]
            d_block = np.where(mismatch, d_block, np.inf)
            flat = d_block.ravel()
            if flat.size == 0:
                continue
            kk = min(k, flat.size)
            part = np.argpartition(flat, kk - 1)[:kk]
            for fi in part:
                d = flat[fi]
                if not np.isfinite(d):
                    continue
                r, c = divmod(int(fi), d_block.shape[1])
                a_idx, b_idx = int(sub_i[r]), int(idx_j[c])
                if len(heap) < k:
                    heapq.heappush(heap, (-d, a_idx, b_idx))
                elif -d > heap[0][0]:
                    heapq.heapreplace(heap, (-d, a_idx, b_idx))

    known_mask = countries != UNKNOWN
    for node in tree.traverse("postorder"):
        children = node.children
        if len(children) < 2:
            continue
        child_known_idx = []
        for child in children:
            idx = np.array([name_to_idx[l.name] for l in child.get_leaves()], dtype=int)
            idx = idx[known_mask[idx]] if idx.size else idx
            child_known_idx.append(idx)
        node_depth = depth[node]
        for ci, cj in itertools.combinations(range(len(children)), 2):
            _consider_block(child_known_idx[ci], child_known_idx[cj], node_depth)

    results = sorted(((-negd, a, b) for negd, a, b in heap), key=lambda t: t[0])
    return [(names[a], names[b], d) for d, a, b in results]


def pair_alphas(dists):
    """alpha = 1 - d/d_max over the given (already top-K) distances, so the
    most distant pair is fully transparent and distance 0 is fully opaque.
    All-zero distances (e.g. a single pair, or all ties at 0) -> alpha 1."""
    if not dists:
        return []
    d_max = max(dists)
    if d_max <= 0:
        return [1.0 for _ in dists]
    return [max(0.0, 1.0 - d / d_max) for d in dists]


# --------------------------------------------------------------------------- #
# Country coordinate lookup (bundled Natural Earth + alias table)            #
# --------------------------------------------------------------------------- #

def _load_ne_countries():
    path = os.path.join(_DATA_DIR, "ne_110m_countries.json")
    with open(path) as fh:
        data = json.load(fh)
    return data["countries"]


def _load_aliases():
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
            if kind == "latlon":
                lat_s, lon_s = value.split(",")
                aliases[insdc_name.lower()] = ("latlon", (float(lat_s), float(lon_s)))
            else:
                aliases[insdc_name.lower()] = ("country", value)
    return aliases


def load_country_coords():
    """Return (land_polygons, name_to_latlon).

    land_polygons: list of rings, each a list of [lon, lat] pairs, for a
    plain background landmass fill.
    name_to_latlon: lower-cased country/location name -> (lat, lon).
    """
    countries = _load_ne_countries()
    aliases = _load_aliases()

    land_polygons = []
    name_to_latlon = {}
    for c in countries:
        for poly in c["polygons"]:
            for ring in poly:
                land_polygons.append(ring)
        if c.get("label_lat") is not None and c.get("label_lon") is not None:
            name_to_latlon[c["name"].lower()] = (c["label_lat"], c["label_lon"])

    # Resolve in repeated passes: a "country" alias may point at another
    # alias (e.g. "Cape Verde" -> "Cabo Verde", itself only a latlon alias,
    # not a direct Natural Earth name) rather than straight at an NE country,
    # so a single pass can miss it depending on dict iteration order.
    pending = dict(aliases)
    for alias_name in list(pending):
        if alias_name in name_to_latlon:
            del pending[alias_name]
    while pending:
        progressed = False
        for alias_name, (kind, value) in list(pending.items()):
            if kind == "latlon":
                name_to_latlon[alias_name] = value
                del pending[alias_name]
                progressed = True
            else:  # kind == "country"
                target = name_to_latlon.get(value.lower())
                if target:
                    name_to_latlon[alias_name] = target
                    del pending[alias_name]
                    progressed = True
        if not progressed:
            break  # remaining entries point at unresolvable names

    return land_polygons, name_to_latlon


def resolve_country_latlon(country_name, name_to_latlon):
    if not country_name or country_name == UNKNOWN:
        return None
    return name_to_latlon.get(country_name.strip().lower())


def _jitter(name, lat, lon, max_radius=1.5):
    """Deterministic small offset so genomes from the same country don't
    all stack on one point. Samples uniformly within a disk of radius
    max_radius (degrees) -- not a square -- so offsets never land farther
    from the country point than max_radius in any direction."""
    h = hashlib.sha256(name.encode("utf-8")).digest()
    u = h[0] / 255.0
    angle = (h[1] / 255.0) * 2 * math.pi
    r = max_radius * math.sqrt(u)
    dlat = r * math.sin(angle)
    dlon = r * math.cos(angle)
    return lat + dlat, lon + dlon


# --------------------------------------------------------------------------- #
# Rendering                                                                   #
# --------------------------------------------------------------------------- #

def render_tree(tree_path, out_path, leaf_country, acc_species, color_map,
                want_species, img_format, dpi, title=None):
    """Render one annotated tree. Returns (country_counts, missing_accessions)."""
    from ete3 import TreeStyle, NodeStyle, RectFace, TextFace

    tree = load_tree(tree_path)
    if tree is None:
        return {}, []

    tree_name = title or os.path.basename(tree_path)
    country_counts = {}
    missing = []
    countries_here = set()

    for leaf in tree.get_leaves():
        acc = leaf.name
        country = leaf_country.get(acc, UNKNOWN)
        if country == UNKNOWN:
            missing.append(acc)
        countries_here.add(country)
        country_counts[country] = country_counts.get(country, 0) + 1

        color = color_map.get(country, color_map[UNKNOWN])

        ns = NodeStyle()
        ns["size"] = 0
        leaf.set_style(ns)

        # optional richer tip label
        if want_species and acc in acc_species and acc_species[acc]:
            leaf.add_face(TextFace(f"{acc}  {acc_species[acc]}", fsize=9),
                          column=0, position="branch-right")

        # aligned color strip keyed to country
        strip = RectFace(width=18, height=14, fgcolor=color, bgcolor=color)
        strip.margin_left = 6
        leaf.add_face(strip, column=1, position="aligned")

    ts = TreeStyle()
    ts.show_leaf_name = not want_species  # avoid double labels when species shown
    ts.show_branch_support = False
    ts.mode = "r"  # rectangular
    ts.scale = None
    ts.title.add_face(TextFace(tree_name, fsize=12, bold=True), column=0)

    # legend: one colored square + country name per country present in THIS tree
    ts.legend_position = 4
    legend_order = sorted(c for c in countries_here if c != UNKNOWN)
    if UNKNOWN in countries_here:
        legend_order.append(UNKNOWN)
    for i, country in enumerate(legend_order):
        color = color_map.get(country, color_map[UNKNOWN])
        swatch = RectFace(width=16, height=16, fgcolor=color, bgcolor=color)
        swatch.margin_right = 4
        swatch.margin_top = 2
        swatch.margin_bottom = 2
        ts.legend.add_face(swatch, column=0)
        label = f"{country} (n={country_counts.get(country, 0)})"
        ts.legend.add_face(TextFace(label, fsize=10), column=1)

    render_kwargs = {"tree_style": ts}
    if img_format == "png":
        render_kwargs["dpi"] = dpi
    tree.render(out_path, **render_kwargs)
    return country_counts, missing


def _dot_area_points2(n):
    """Scatter marker area (points^2) for a country with n samples."""
    return 25 * (n ** 0.5)


def render_map(out_path, leaf_country, country_counts, pairs, alphas, color_map,
               land_polygons, name_to_latlon, img_format, dpi, title=None):
    """Render the world map: country points sized by sample count, plus
    top-K closest cross-country genome-pair lines. Returns the set of
    countries that could not be geocoded."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PolyCollection

    unresolved = set()
    resolved = {}
    for country in country_counts:
        if country == UNKNOWN:
            continue
        coord = resolve_country_latlon(country, name_to_latlon)
        if coord is None:
            unresolved.add(country)
        else:
            resolved[country] = coord

    fig, ax = plt.subplots(figsize=(14, 7))

    if land_polygons:
        pc = PolyCollection(land_polygons, facecolor="#e8e8e0",
                            edgecolor="#b0b0a8", linewidths=0.3, zorder=1)
        ax.add_collection(pc)

    # Axis limits/aspect must be finalized before converting between data
    # (degrees) and display (points) units below, so jitter radii can be
    # expressed in degrees but bounded by each country's actual on-screen
    # dot radius -- otherwise jittered line endpoints can land outside the
    # dot they belong to.
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("equal")
    ax.apply_aspect()

    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((1, 0))
    p2 = ax.transData.transform((0, 1))
    points_per_degree = (
        (abs(p1[0] - p0[0]) + abs(p2[1] - p0[1])) / 2.0
    ) * 72.0 / fig.dpi

    def dot_radius_degrees(n):
        radius_points = math.sqrt(_dot_area_points2(n) / math.pi)
        return radius_points / points_per_degree

    jittered = {}
    for name, country in leaf_country.items():
        if country in resolved:
            lat, lon = resolved[country]
            n = country_counts.get(country, 0)
            # Keep every jittered sample point strictly inside its
            # country's rendered dot.
            jittered[name] = _jitter(name, lat, lon, max_radius=0.6 * dot_radius_degrees(n))

    if pairs:
        segments = []
        colors = []
        for (a, b, _d), alpha in zip(pairs, alphas):
            if a not in jittered or b not in jittered:
                continue
            lat_a, lon_a = jittered[a]
            lat_b, lon_b = jittered[b]
            segments.append([(lon_a, lat_a), (lon_b, lat_b)])
            colors.append((0.05, 0.05, 0.05, alpha))
        if segments:
            lc = LineCollection(segments, colors=colors, linewidths=0.9, zorder=2)
            ax.add_collection(lc)

    legend_order = sorted(resolved, key=lambda c: -country_counts.get(c, 0))
    for country in legend_order:
        lat, lon = resolved[country]
        n = country_counts.get(country, 0)
        ax.scatter([lon], [lat], s=_dot_area_points2(n),
                   color=color_map.get(country, "#999999"),
                   edgecolor="black", linewidths=0.5, zorder=3)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=8,
                    color=color_map.get(c, "#999999"),
                    label=f"{c} (n={country_counts.get(c, 0)})")
        for c in legend_order
    ]
    if handles:
        # A single lower-left column works fine for a handful of countries,
        # but at dense-map scale (dozens of countries) that column grows
        # taller than the map itself and covers data. Past that point, lay
        # the legend out in several columns below the map instead.
        if len(handles) > 20:
            ncol = math.ceil(len(handles) / 12)
            ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02),
                       fontsize=8, framealpha=0.9, ncol=ncol)
        else:
            ax.legend(handles=handles, loc="lower left", fontsize=8,
                       framealpha=0.9, ncol=1)

    save_kwargs = {"bbox_inches": "tight"}
    if img_format == "png":
        save_kwargs["dpi"] = dpi
    fig.savefig(out_path, **save_kwargs)
    plt.close(fig)
    return unresolved


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def _collect_tree_paths(trees_dir):
    patterns = ("*.nwk", "*.tree", "*.treefile")
    seen = set()
    paths = []
    for pat in patterns:
        for p in sorted(glob.glob(os.path.join(trees_dir, pat))):
            if p not in seen:
                seen.add(p)
                paths.append(p)
    return sorted(paths)


def main():
    args = parse_args()

    if args.tree:
        tree_path = os.path.abspath(args.tree)
        if not os.path.isfile(tree_path):
            eprint(f"ERROR: --tree file not found: {tree_path}")
            sys.exit(1)
        tree_paths = [tree_path]
        out_dir = args.out_dir or os.path.join(os.path.dirname(tree_path), "geo_annotated")
        labels = {tree_path: args.label or os.path.splitext(os.path.basename(tree_path))[0]}
    else:
        trees_dir = os.path.abspath(args.trees_dir)
        if not os.path.isdir(trees_dir):
            eprint(f"ERROR: --trees-dir is not a directory: {trees_dir}")
            sys.exit(1)
        tree_paths = _collect_tree_paths(trees_dir)
        if not tree_paths:
            eprint(f"ERROR: no *.nwk/*.tree/*.treefile files found in {trees_dir}")
            sys.exit(1)
        out_dir = args.out_dir or os.path.join(trees_dir, "geo_annotated")
        labels = {tp: os.path.splitext(os.path.basename(tp))[0] for tp in tree_paths}

    if args.metadata and not os.path.isfile(args.metadata):
        eprint(f"ERROR: --metadata file not found: {args.metadata}")
        sys.exit(1)
    if args.geo_tsv and not os.path.isfile(args.geo_tsv):
        eprint(f"ERROR: --geo-tsv file not found: {args.geo_tsv}")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)

    # ete3 is heavy and needs Qt; import here so --help works without it.
    try:
        import ete3  # noqa: F401
    except Exception as exc:
        eprint(f"ERROR: could not import ete3 ({exc}).")
        eprint("       Activate the OrthoPhyl conda env (ships ete3=3.1.3) and retry.")
        sys.exit(1)

    acc_country = {}
    acc_species = {}
    if args.metadata:
        acc_country, acc_species = parse_metadata(
            args.metadata, args.meta_acc_col, args.meta_country_col,
            args.species_col, args.label_with_species)

    geo_tsv_country = {}
    if args.geo_tsv:
        geo_tsv_country = parse_geo_tsv(args.geo_tsv)

    land_polygons, name_to_latlon = ([], {})
    if not args.no_map:
        try:
            land_polygons, name_to_latlon = load_country_coords()
        except Exception as exc:
            eprint(f"WARNING: could not load bundled country coordinate data ({exc}); "
                   f"disabling map rendering")
            args.no_map = True

    # First pass: collect every country that actually appears on a tip across
    # all trees, so colors stay consistent from figure to figure.
    trees_by_path = {}
    present = set()
    for tp in tree_paths:
        t = load_tree(tp)
        trees_by_path[tp] = t
        if t is None:
            continue
        leaf_names = [leaf.name for leaf in t.get_leaves()]
        leaf_country, _ = build_leaf_country_map(leaf_names, acc_country, geo_tsv_country)
        present.update(leaf_country.values())
    color_map = build_color_map(present)

    total_counts = {}
    all_missing = set()
    all_unresolved = set()
    rendered = 0
    for tp in tree_paths:
        tree = trees_by_path.get(tp)
        if tree is None:
            continue
        label = labels[tp]
        leaf_names = [leaf.name for leaf in tree.get_leaves()]
        leaf_country, missing = build_leaf_country_map(leaf_names, acc_country, geo_tsv_country)

        tree_out = os.path.join(out_dir, f"{label}_geo_tree.{args.format}")
        counts, _ = render_tree(
            tp, tree_out, leaf_country, acc_species, color_map,
            args.label_with_species, args.format, args.dpi, title=label)
        if not counts:
            continue
        rendered += 1
        print(f"[rendered] {label} -> {tree_out}")
        for c, n in counts.items():
            total_counts[c] = total_counts.get(c, 0) + n
        all_missing.update(missing)

        pairs = top_k_cross_country_pairs(tree, leaf_country, args.top_k)
        alphas = pair_alphas([d for _, _, d in pairs])

        pairs_path = os.path.join(out_dir, f"{label}_geo_pairs.tsv")
        with open(pairs_path, "w") as fh:
            fh.write("genome_a\tgenome_b\tcountry_a\tcountry_b\tdistance\talpha\n")
            for (a, b, d), alpha in zip(pairs, alphas):
                fh.write(f"{a}\t{b}\t{leaf_country.get(a, UNKNOWN)}\t"
                         f"{leaf_country.get(b, UNKNOWN)}\t{d:.6g}\t{alpha:.4f}\n")

        unresolved = set()
        if not args.no_map:
            map_out = os.path.join(out_dir, f"{label}_geo_map.{args.format}")
            unresolved = render_map(
                map_out, leaf_country, counts, pairs, alphas, color_map,
                land_polygons, name_to_latlon, args.format, args.dpi, title=label)
            print(f"[rendered] {label} -> {map_out}")
            all_unresolved.update(unresolved)

        countries_path = os.path.join(out_dir, f"{label}_geo_countries.tsv")
        with open(countries_path, "w") as fh:
            fh.write("country\tn_genomes\tlat\tlon\n")
            for country, n in sorted(counts.items(), key=lambda kv: (kv[0] == UNKNOWN, -kv[1])):
                coord = None if country == UNKNOWN else resolve_country_latlon(country, name_to_latlon)
                if coord:
                    fh.write(f"{country}\t{n}\t{coord[0]:.4f}\t{coord[1]:.4f}\n")
                else:
                    fh.write(f"{country}\t{n}\tunresolved\tunresolved\n")

    # Summary report
    print("\n===== summary =====")
    print(f"trees rendered: {rendered}/{len(tree_paths)}")
    print(f"output dir:     {out_dir}")
    print("country -> tip count (across all trees):")
    for c in sorted(total_counts, key=lambda k: (k == UNKNOWN, k)):
        print(f"  {c}: {total_counts[c]}")
    if all_missing:
        print(f"\nWARNING: {len(all_missing)} tip label(s) had no metadata/geo-tsv match "
              f"(bucketed as {UNKNOWN}):")
        for acc in sorted(all_missing):
            print(f"  {acc}")
        print(f"  -> check that tree leaf names match column {args.meta_acc_col} of "
              "the metadata, or add them to --geo-tsv.")
    if all_unresolved:
        print(f"\nWARNING: {len(all_unresolved)} country name(s) could not be geocoded "
              "(shown in the tree but omitted from the map):")
        for c in sorted(all_unresolved):
            print(f"  {c}")
        print("  -> add a row to python_scripts/data/country_aliases.tsv.")


if __name__ == "__main__":
    main()
