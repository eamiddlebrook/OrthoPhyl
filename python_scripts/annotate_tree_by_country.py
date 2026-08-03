#!/usr/bin/env python3
"""
Annotate OrthoPhyl species trees with the country of origin of each assembly.

OrthoPhyl writes its final Newick trees to $store/FINAL_SPECIES_TREES/*.tree .
When the input assemblies were gathered with utils/gather_filter_asms.sh, each
tree leaf label is the assembly accession (GCF_/GCA_...), which is column 1 of
the gather output file `all_asm_acc_metadata`.  Column 8 of that same file is the
country / geo-location of origin (underscore-joined, e.g. United_Kingdom).

This script reads the metadata and every *.tree in a directory (typically
FINAL_SPECIES_TREES) and renders a figure per tree where each tip carries a
colored strip keyed to its country of origin, plus a legend.

Requires the OrthoPhyl conda env (ships ete3 + matplotlib).  Renders headlessly
via offscreen Qt -- no $DISPLAY needed.

Example:
    python annotate_tree_by_country.py \\
        --trees-dir  /path/to/FINAL_SPECIES_TREES \\
        --metadata   /path/to/all_asm_acc_metadata \\
        --out-dir    /path/to/FINAL_SPECIES_TREES/country_annotated
"""

import argparse
import colorsys
import glob
import os
import sys

# ete3 needs a Qt backend; force offscreen so it works with no $DISPLAY.
# Must be set before ete3 (and therefore Qt) is imported.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

UNKNOWN = "Unknown"
# values in the metadata's country column that mean "no real location"
_UNKNOWN_VALUES = {"", "na", "unknown", "none", "null", "-"}


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate OrthoPhyl species trees with assembly country of origin.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--trees-dir", required=True,
                   help="Directory containing *.tree Newick files "
                        "(e.g. a FINAL_SPECIES_TREES directory).")
    p.add_argument("--metadata", required=True,
                   help="Path to all_asm_acc_metadata from gather_filter_asms.sh.")
    p.add_argument("--out-dir", default=None,
                   help="Directory for output figures "
                        "(default: <trees-dir>/country_annotated).")
    p.add_argument("--format", default="png", choices=["png", "svg", "pdf"],
                   help="Output image format.")
    p.add_argument("--meta-acc-col", type=int, default=1,
                   help="1-based column in metadata holding the accession / leaf name.")
    p.add_argument("--meta-country-col", type=int, default=8,
                   help="1-based column in metadata holding the country of origin.")
    p.add_argument("--label-with-species", action="store_true",
                   help="Append the species.strain field (metadata col 3) to each tip label.")
    p.add_argument("--species-col", type=int, default=3,
                   help="1-based column used for --label-with-species.")
    p.add_argument("--dpi", type=int, default=300, help="Raster DPI for png output.")
    return p.parse_args()


def clean_value(value):
    """Strip surrounding quotes/whitespace from a metadata token."""
    return value.strip().strip('"').strip("'")


def normalize_country(value):
    """Map raw metadata country token to a display bucket; missing -> Unknown."""
    v = clean_value(value)
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
        from matplotlib import cm
        for cmap_name in ("tab20", "tab20b", "tab20c"):
            cmap = cm.get_cmap(cmap_name)
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


def render_tree(tree_path, out_path, acc_country, acc_species, color_map,
                want_species, img_format, dpi):
    """Render one annotated tree. Returns (country_counts, missing_accessions)."""
    from ete3 import TreeStyle, NodeStyle, RectFace, TextFace

    tree = load_tree(tree_path)
    if tree is None:
        return {}, []

    tree_name = os.path.basename(tree_path)
    country_counts = {}
    missing = []
    countries_here = set()

    for leaf in tree.get_leaves():
        acc = leaf.name
        if acc in acc_country:
            country = acc_country[acc]
        else:
            country = UNKNOWN
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


def main():
    args = parse_args()

    trees_dir = os.path.abspath(args.trees_dir)
    if not os.path.isdir(trees_dir):
        eprint(f"ERROR: --trees-dir is not a directory: {trees_dir}")
        sys.exit(1)

    tree_paths = sorted(glob.glob(os.path.join(trees_dir, "*.tree")))
    if not tree_paths:
        eprint(f"ERROR: no *.tree files found in {trees_dir}")
        sys.exit(1)

    if not os.path.isfile(args.metadata):
        eprint(f"ERROR: --metadata file not found: {args.metadata}")
        sys.exit(1)

    out_dir = args.out_dir or os.path.join(trees_dir, "country_annotated")
    os.makedirs(out_dir, exist_ok=True)

    # ete3 is heavy and needs Qt; import here so --help works without it.
    try:
        import ete3  # noqa: F401
    except Exception as exc:
        eprint(f"ERROR: could not import ete3 ({exc}).")
        eprint("       Activate the OrthoPhyl conda env (ships ete3=3.1.3) and retry.")
        sys.exit(1)

    acc_country, acc_species = parse_metadata(
        args.metadata, args.meta_acc_col, args.meta_country_col,
        args.species_col, args.label_with_species)

    # First pass: collect every country that actually appears on a tip across all
    # trees, so colors are consistent from figure to figure.
    from ete3 import Tree  # noqa: F401  (ensures ete3 fully importable)
    present = set()
    for tp in tree_paths:
        t = load_tree(tp)
        if t is None:
            continue
        for leaf in t.get_leaves():
            present.add(acc_country.get(leaf.name, UNKNOWN))
    color_map = build_color_map(present)

    total_counts = {}
    all_missing = set()
    rendered = 0
    for tp in tree_paths:
        base = os.path.basename(tp)
        out_path = os.path.join(out_dir, f"{base}.country.{args.format}")
        counts, missing = render_tree(
            tp, out_path, acc_country, acc_species, color_map,
            args.label_with_species, args.format, args.dpi)
        if counts:
            rendered += 1
            print(f"[rendered] {base} -> {out_path}")
            for c, n in counts.items():
                total_counts[c] = total_counts.get(c, 0) + n
            all_missing.update(missing)

    # Summary report
    print("\n===== summary =====")
    print(f"trees rendered: {rendered}/{len(tree_paths)}")
    print(f"output dir:     {out_dir}")
    print("country -> tip count (across all trees):")
    for c in sorted(total_counts, key=lambda k: (k == UNKNOWN, k)):
        print(f"  {c}: {total_counts[c]}")
    if all_missing:
        print(f"\nWARNING: {len(all_missing)} tip label(s) had no metadata match "
              f"(bucketed as {UNKNOWN}):")
        for acc in sorted(all_missing):
            print(f"  {acc}")
        print("  -> check that tree leaf names match column "
              f"{args.meta_acc_col} of the metadata.")


if __name__ == "__main__":
    main()
