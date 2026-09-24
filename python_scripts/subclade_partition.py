#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
subclade_partition.py

Partition an oversized set of genomes into size-bounded "subclades" using MASH
sequence distances, so a huge taxon can be split into several manageable trees
instead of one intractable one.

Given a directory of genome FASTAs and a target maximum subclade size N, this:
  1. runs `mash triangle` (all-vs-all) to get pairwise distances,
  2. builds a distance matrix and clusters with average-linkage (UPGMA),
  3. recursively splits the tree so every subclade has <= N genomes,
  4. merges any too-small (< min-size) cluster into its nearest neighbor,
  5. numbers subclades deterministically as <Taxon>_1, <Taxon>_2, ...,
  6. writes one combined MASH sketch (.msh) per subclade for the router to
     compare future queries against,
  7. records which subclade each --query genome landed in,
  8. emits partition_manifest.json describing the result.

MASH parameters MUST match those used elsewhere in OrthoPhyl (script_lib/
functions.sh: `mash triangle -k 17 -s 5000`) so the per-subclade sketches this
writes are directly comparable to a query sketched with the same params by the
assembly router.

The distance matrix is kept in scipy's CONDENSED (1-D, upper-triangle-only) form
throughout -- parse_mash_edges builds it directly and partition_matrix/
_merge_tiny_clusters consume it via condensed_index -- rather than ever
materializing a dense n x n matrix. This roughly halves peak memory at large n
(measured: ~45 GB dense vs ~22 GB condensed for a real 74,707-genome taxon).
`mash triangle` itself is still O(n^2) in TIME regardless (unaffected by this).

The clustering logic is separated from the mash calls (see partition_matrix) so
it can be unit-tested by injecting a pre-computed condensed distance array
without running mash.

Usage:
    subclade_partition.py --genome-dir DIR --taxon NAME --out-dir DIR \
        --max-size 150 [--min-size 4] [--threads 1] [--query STEM ...]
"""
import argparse
import glob
import json
import os
import subprocess
import sys

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree

# MASH sketch parameters -- MUST match script_lib/functions.sh so distances and
#   sketches are comparable across OrthoPhyl, this tool, and the router.
MASH_K = "17"
MASH_S = "5000"

# mash reports distances in [0, 1]; use the max (1.0) for pairs it could not
#   compare (too divergent to share hashes). NO artificial floor (unlike the old
#   ANI_genome_picking.py, which used 50 for ANI% -- wrong axis here).
MASH_MAX_DIST = 1.0

# FASTA extensions to consider as genomes.
FASTA_GLOBS = ("*.fna", "*.fasta", "*.fa")


def condensed_index(i, j, n):
    """
    Index into a scipy condensed distance array for pair (i, j), 0 <= i,j < n,
    i != j. Matches scipy.spatial.distance.squareform's own ordering exactly --
    the condensed array here is consumed directly by scipy.cluster.hierarchy.
    linkage(), which requires that exact convention.
    """
    if i == j:
        raise ValueError("condensed_index has no entry for the diagonal (i == j)")
    if i > j:
        i, j = j, i
    return n * i - i * (i + 1) // 2 + (j - i - 1)


def _pairwise_dist(condensed, i, j, n):
    """Look up the distance between genome indices i and j in a condensed array."""
    return condensed[condensed_index(i, j, n)]


def find_genomes(genome_dir):
    """Return a sorted list of genome FASTA paths in genome_dir (deterministic)."""
    files = []
    for pat in FASTA_GLOBS:
        files.extend(glob.glob(os.path.join(genome_dir, pat)))
    # de-dup (a file matching two globs) and sort for determinism
    return sorted(set(files))


def _run_mash(cmd, stdout_path=None):
    """Run a mash command (shell=False). If stdout_path is given, redirect stdout there."""
    if stdout_path is not None:
        with open(stdout_path, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh)
    else:
        subprocess.run(cmd, check=True)


def run_mash_triangle(genome_files, out_path, threads):
    """Run all-vs-all `mash triangle -E` on genome_files, writing edges to out_path."""
    cmd = [
        "mash", "triangle",
        "-k", MASH_K, "-s", MASH_S,
        "-E", "-p", str(threads),
    ] + list(genome_files)
    _run_mash(cmd, stdout_path=out_path)


def parse_mash_edges(edge_path, names):
    """
    Parse `mash triangle -E` output into a condensed distance array (scipy's
    1-D upper-triangle-only form -- see condensed_index).

    Each edge line is: <seqA> <seqB> <dist> <p-value> <shared-hashes>
    (seqA/seqB are the paths mash was given). Returns a 1-D numpy array of
    length n*(n-1)/2, indexed via condensed_index(i, j, n) where i/j are
    positions in `names` (basenames). Pairs absent from the edge list default
    to MASH_MAX_DIST. There is no diagonal in condensed form (always 0,
    never looked up).
    """
    idx = {name: i for i, name in enumerate(names)}
    n = len(names)
    condensed = np.full(n * (n - 1) // 2, MASH_MAX_DIST, dtype=float)
    with open(edge_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            a = os.path.basename(parts[0])
            b = os.path.basename(parts[1])
            if a not in idx or b not in idx:
                continue
            try:
                dist = float(parts[2])
            except ValueError:
                continue
            i, j = idx[a], idx[b]
            if i == j:
                continue
            condensed[condensed_index(i, j, n)] = dist
    return condensed


def _split_tree(node, max_size):
    """
    Recursively descend a scipy cluster tree, emitting a node's leaf-id list as a
    subclade as soon as its leaf count <= max_size. Returns a list of leaf-id lists.
    """
    leaves = node.pre_order(lambda x: x.id)
    if len(leaves) <= max_size:
        return [leaves]
    # Too big: recurse into both children (a non-leaf node always has both).
    return _split_tree(node.left, max_size) + _split_tree(node.right, max_size)


def _merge_tiny_clusters(clusters, condensed, n, min_size):
    """
    Merge clusters smaller than min_size into the (non-tiny) cluster holding
    their single nearest member (min pairwise MASH distance). Iterates until no
    tiny cluster remains that can be merged. `clusters` is a list of index-lists.
    `condensed`/`n` are the condensed distance array and total genome count (see
    condensed_index) used to look up pairwise distances. Returns a new list of
    index-lists.
    """
    clusters = [list(c) for c in clusters]
    while True:
        big = [c for c in clusters if len(c) >= min_size]
        tiny = [c for c in clusters if len(c) < min_size]
        if not tiny:
            break
        if not big:
            # Nothing to merge into (every cluster is tiny) -- give up, caller
            #   handles the "can't partition" fallback.
            break
        # Merge the single smallest tiny cluster into its nearest big cluster,
        #   then re-loop (a merge can turn a big cluster even bigger, fine).
        tiny.sort(key=lambda c: (len(c), min(c)))
        t = tiny[0]
        best_big = None
        best_dist = None
        for c in big:
            d = min(_pairwise_dist(condensed, i, j, n) for i in t for j in c)
            if best_dist is None or d < best_dist:
                best_dist = d
                best_big = c
        best_big.extend(t)
        clusters = [c for c in clusters if c is not t and c is not best_big]
        clusters.append(best_big)
    return clusters


def partition_matrix(names, condensed, max_size, min_size):
    """
    Cluster `names` (given their condensed distance array -- see condensed_index)
    into subclades of <= max_size, merging clusters < min_size. Returns a list of
    subclades, each a sorted list of names, ordered deterministically (size desc,
    then min name).

    Pure function of (names, condensed) -- no mash, so this is directly
    unit-testable. `condensed` is scipy's condensed (1-D, upper-triangle-only)
    form, consumed directly by linkage() -- never materialized as a dense
    matrix, which would roughly double peak memory at large n.
    """
    n = len(names)
    if n <= max_size:
        return [sorted(names)]

    Z = linkage(condensed, method="average")
    tree = to_tree(Z)
    id_clusters = _split_tree(tree, max_size)
    id_clusters = _merge_tiny_clusters(id_clusters, condensed, n, min_size)

    name_clusters = [sorted(names[i] for i in c) for c in id_clusters]
    # Deterministic ordering: largest first, ties broken by lexical min member.
    name_clusters.sort(key=lambda c: (-len(c), c[0]))
    return name_clusters


def write_sketch(subclade_name, member_paths, out_dir, threads):
    """Write one combined MASH sketch for a subclade. Returns the .msh path."""
    prefix = os.path.join(out_dir, subclade_name)
    cmd = [
        "mash", "sketch",
        "-k", MASH_K, "-s", MASH_S,
        "-p", str(threads),
        "-o", prefix,
    ] + list(member_paths)
    _run_mash(cmd)
    return prefix + ".msh"


def partition(genome_dir, taxon, out_dir, max_size, min_size, threads, queries):
    """
    Top-level orchestration: glob genomes, run mash, cluster, write per-subclade
    sketches + member files, and return the manifest dict (also written to disk).
    """
    os.makedirs(out_dir, exist_ok=True)
    genome_files = find_genomes(genome_dir)
    if not genome_files:
        raise SystemExit("ERROR: no genome FASTAs (%s) found in %s"
                         % ("/".join(FASTA_GLOBS), genome_dir))

    names = [os.path.basename(f) for f in genome_files]
    path_by_name = dict(zip(names, genome_files))
    queries = queries or []

    if len(names) <= max_size:
        # No partitioning needed: one un-suffixed "subclade" covering everything.
        manifest = _build_manifest(
            partitioned=False, taxon=taxon, max_size=max_size,
            subclades=[(taxon, sorted(names))],
            out_dir=out_dir, path_by_name=path_by_name,
            threads=threads, queries=queries, write_sketches=False,
        )
        _write_manifest(manifest, out_dir)
        return manifest

    mash_out = os.path.join(out_dir, "MASH_out")
    run_mash_triangle(genome_files, mash_out, threads)
    condensed = parse_mash_edges(mash_out, names)

    clusters = partition_matrix(names, condensed, max_size, min_size)

    # If clustering could not produce any subclade meeting the min-size floor
    #   (e.g. a single genome set that all merged tiny), fall back to unpartitioned.
    if len(clusters) == 1 and len(clusters[0]) < min_size:
        manifest = _build_manifest(
            partitioned=False, taxon=taxon, max_size=max_size,
            subclades=[(taxon, sorted(names))],
            out_dir=out_dir, path_by_name=path_by_name,
            threads=threads, queries=queries, write_sketches=False,
        )
        _write_manifest(manifest, out_dir)
        return manifest

    named = [("%s_%d" % (taxon, i + 1), members) for i, members in enumerate(clusters)]
    manifest = _build_manifest(
        partitioned=True, taxon=taxon, max_size=max_size,
        subclades=named, out_dir=out_dir, path_by_name=path_by_name,
        threads=threads, queries=queries, write_sketches=True,
    )
    _write_manifest(manifest, out_dir)
    return manifest


def _build_manifest(partitioned, taxon, max_size, subclades, out_dir,
                    path_by_name, threads, queries, write_sketches):
    """Assemble the manifest dict, writing member files (+ sketches if requested)."""
    query_set = set(queries)
    query_assignments = {}
    subclade_entries = []
    for i, (name, members) in enumerate(subclades):
        members_file = os.path.join(out_dir, name + ".members.txt")
        with open(members_file, "w") as fh:
            for m in members:
                fh.write(m + "\n")
        sketch_file = None
        if write_sketches:
            member_paths = [path_by_name[m] for m in members]
            sketch_file = write_sketch(name, member_paths, out_dir, threads)
        for q in query_set:
            if _matches_member(q, members):
                query_assignments[q] = name
        subclade_entries.append({
            "subclade_id": i + 1,
            "name": name,
            "n_genomes": len(members),
            "members_file": members_file,
            "sketch_file": sketch_file,
        })
    return {
        "partitioned": partitioned,
        "parent_taxon": taxon,
        "max_size": max_size,
        "n_subclades": len(subclade_entries),
        "subclades": subclade_entries,
        "query_assignments": query_assignments,
    }


def _matches_member(query_stem, members):
    """
    True if a --query stem corresponds to one of `members` (which are basenames
    like GCF_x.fna). Matches on the stem (extension stripped) either way so the
    caller can pass a bare accession or a filename.
    """
    q = _strip_ext(query_stem)
    for m in members:
        if m == query_stem or _strip_ext(m) == q:
            return True
    return False


def _strip_ext(name):
    base = os.path.basename(name)
    for ext in (".fna", ".fasta", ".fa"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def _write_manifest(manifest, out_dir):
    path = os.path.join(out_dir, "partition_manifest.json")
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Partition a genome set into size-bounded MASH subclades.")
    p.add_argument("--genome-dir", required=True,
                   help="Directory of genome FASTAs (.fna/.fasta/.fa) to partition.")
    p.add_argument("--taxon", required=True,
                   help="Parent taxon name; subclades are named <taxon>_1, _2, ...")
    p.add_argument("--out-dir", required=True,
                   help="Directory for MASH_out, per-subclade .msh/.members.txt, "
                        "and partition_manifest.json.")
    p.add_argument("--max-size", type=int, default=150,
                   help="Target maximum genomes per subclade (default 150).")
    p.add_argument("--min-size", type=int, default=4,
                   help="Minimum genomes per subclade; smaller clusters are merged "
                        "into their nearest neighbor (default 4, OrthoPhyl's floor).")
    p.add_argument("--threads", type=int, default=1,
                   help="Threads for mash (default 1).")
    p.add_argument("--query", action="append", default=[],
                   help="A query genome stem/basename to locate in the partition; "
                        "repeatable. Recorded in the manifest's query_assignments.")
    args = p.parse_args(argv)

    manifest = partition(
        genome_dir=args.genome_dir, taxon=args.taxon, out_dir=args.out_dir,
        max_size=args.max_size, min_size=args.min_size, threads=args.threads,
        queries=args.query,
    )
    print("Partitioned=%s: %d subclade(s) from taxon '%s'"
          % (manifest["partitioned"], manifest["n_subclades"], args.taxon))
    for sc in manifest["subclades"]:
        print("  %s: %d genomes" % (sc["name"], sc["n_genomes"]))
    if args.query:
        for q, name in sorted(manifest["query_assignments"].items()):
            print("  query %s -> %s" % (q, name))
    return 0


if __name__ == "__main__":
    sys.exit(main())
