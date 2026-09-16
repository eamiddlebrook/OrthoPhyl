#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
subsample_genomes.py

Pick a diverse, size-bounded subset of a large genome set using MASH sequence
distances, so a huge taxon can be reduced to ONE manageable tree instead of an
intractably large one (or an O(n^2) partition).

Given a directory of genome FASTAs and a target count N, this:
  1. sketches ALL genomes once with `mash sketch` (linear in the number of
     genomes -- no all-vs-all matrix),
  2. greedily selects a maximally-diverse subset via farthest-point sampling
     (max-min): repeatedly add the genome whose nearest already-selected genome
     is the most distant,
  3. optionally seeds the selection with must-keep genomes (e.g. query genomes)
     so they are always retained,
  4. writes subsample_members.txt (chosen basenames) and subsample_manifest.json.

Why greedy max-min instead of clustering: farthest-point sampling needs only one
`mash dist <candidate> combined.msh` call per pick (a single row of distances),
giving O(n*N) time and O(n) memory. subclade_partition.py's approach builds a
dense N_total x N_total distance matrix (O(n^2) memory) which OOMs on very large
taxa -- exactly the case this tool exists to handle.

MASH parameters MUST match those used elsewhere in OrthoPhyl (script_lib/
functions.sh: `mash ... -k 17 -s 5000`) so sketches/distances are comparable.

The selection logic is separated from the mash calls (see greedy_maxmin) so it
can be unit-tested by injecting a distance oracle without running mash.

Usage:
    subsample_genomes.py --genome-dir DIR --out-dir DIR --n 500 \
        [--must-keep STEM ...] [--threads 1]
"""
import argparse
import glob
import json
import math
import os
import subprocess
import sys

# MASH sketch parameters -- MUST match script_lib/functions.sh and
#   subclade_partition.py so distances/sketches are comparable across tools.
MASH_K = "17"
MASH_S = "5000"

# FASTA extensions to consider as genomes.
FASTA_GLOBS = ("*.fna", "*.fasta", "*.fa")


def find_genomes(genome_dir):
    """Return a sorted list of genome FASTA paths in genome_dir (deterministic)."""
    files = []
    for pat in FASTA_GLOBS:
        files.extend(glob.glob(os.path.join(genome_dir, pat)))
    # de-dup (a file matching two globs) and sort for determinism
    return sorted(set(files))


def _strip_ext(name):
    base = os.path.basename(name)
    for ext in (".fna", ".fasta", ".fa"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def _run_mash(cmd, stdout_path=None):
    """Run a mash command (shell=False). If stdout_path is given, redirect stdout there."""
    if stdout_path is not None:
        with open(stdout_path, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh)
    else:
        subprocess.run(cmd, check=True)


def greedy_maxmin(names, dist_fn, target, seed_names=None):
    """
    Select `target` maximally-diverse names via farthest-point (max-min) sampling.

    Pure function of (names, dist_fn) -- no mash -- so it is directly unit-testable
    by injecting a distance oracle.

    Args:
      names:      list of candidate names (basenames). Order does not matter; the
                  algorithm sorts internally for determinism.
      dist_fn:    callable name -> {other_name: distance}. Distances for pairs the
                  oracle omits default to +inf (treated as maximally divergent).
      target:     desired number of selected names.
      seed_names: names that MUST be selected first (e.g. query genomes). Any not
                  in `names` are ignored. If None/empty, the lexicographically
                  first name seeds the selection (deterministic).

    Returns a list of selected names in selection order (seeds first). If
    target >= len(names) all names are returned (sorted).
    """
    names = sorted(set(names))
    n = len(names)
    if target >= n:
        return list(names)
    if target <= 0:
        return []

    name_set = set(names)
    seeds = [s for s in (seed_names or []) if s in name_set]
    # De-dup seeds preserving order.
    seen = set()
    seeds = [s for s in seeds if not (s in seen or seen.add(s))]

    if not seeds:
        seeds = [names[0]]

    selected = []
    selected_set = set()
    # min_dist[name] = distance from `name` to its nearest selected genome.
    min_dist = {nm: math.inf for nm in names}

    def _select(name):
        selected.append(name)
        selected_set.add(name)
        row = dist_fn(name)
        for other in names:
            if other in selected_set:
                continue
            d = row.get(other, math.inf)
            if d < min_dist[other]:
                min_dist[other] = d

    for s in seeds:
        if s not in selected_set and len(selected) < target:
            _select(s)

    while len(selected) < target:
        # Pick the unselected genome farthest from the current selection. `names`
        # is sorted and we compare with strict `>`, so on a distance tie the
        # lexicographically-first name wins -- deterministic tie-break.
        best = None
        best_dist = None
        for nm in names:
            if nm in selected_set:
                continue
            if best is None or min_dist[nm] > best_dist:
                best = nm
                best_dist = min_dist[nm]
        if best is None:
            break
        _select(best)

    return selected


def _mash_dist_row(query_path, combined_msh):
    """
    Run `mash dist <query_path> <combined_msh>` and return {member_name: dist}.

    mash dist output lines: <ref> <query> <dist> <p-value> <shared-hashes>, one
    per sketched reference in combined_msh. We map each ref path back to its
    basename so keys align with `names`.
    """
    proc = subprocess.run(
        ["mash", "dist", str(query_path), str(combined_msh)],
        check=True, capture_output=True, text=True,
    )
    row = {}
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        ref = os.path.basename(parts[0])
        try:
            dist = float(parts[2])
        except ValueError:
            continue
        row[ref] = dist
    return row


def subsample(genome_dir, out_dir, target, threads, must_keep=None):
    """
    Top-level orchestration: glob genomes, sketch once, greedily pick `target`
    diverse genomes (seeded by must_keep), write members + manifest.

    Returns the manifest dict (also written to disk).
    """
    os.makedirs(out_dir, exist_ok=True)
    genome_files = find_genomes(genome_dir)
    if not genome_files:
        raise SystemExit("ERROR: no genome FASTAs (%s) found in %s"
                         % ("/".join(FASTA_GLOBS), genome_dir))

    names = [os.path.basename(f) for f in genome_files]
    path_by_name = dict(zip(names, genome_files))
    n_total = len(names)

    # Normalize must-keep (accepts bare accession or filename) to member basenames.
    stem_to_name = {_strip_ext(nm): nm for nm in names}
    seed_names = []
    for mk in (must_keep or []):
        cand = stem_to_name.get(_strip_ext(mk))
        if cand:
            seed_names.append(cand)

    if n_total <= target:
        # Nothing to do: keep everything.
        manifest = _build_manifest(
            subsampled=False, out_dir=out_dir, target=target,
            n_total=n_total, members=sorted(names), seeds=sorted(set(seed_names)))
        _write_manifest(manifest, out_dir)
        return manifest

    combined = os.path.join(out_dir, "combined")
    # Pass genome paths via a file-of-filenames (-l), not as bare argv: at tens
    # of thousands of genomes the combined path list overflows execve's
    # ARG_MAX (confirmed: 74707 paths fails with "Argument list too long").
    filelist = os.path.join(out_dir, "sketch_input.txt")
    with open(filelist, "w") as fh:
        for f in genome_files:
            fh.write(f + "\n")
    cmd = [
        "mash", "sketch",
        "-k", MASH_K, "-s", MASH_S,
        "-p", str(threads),
        "-o", combined,
        "-l", filelist,
    ]
    _run_mash(cmd)
    combined_msh = combined + ".msh"

    def dist_fn(name):
        return _mash_dist_row(path_by_name[name], combined_msh)

    selected = greedy_maxmin(names, dist_fn, target, seed_names=seed_names)

    manifest = _build_manifest(
        subsampled=True, out_dir=out_dir, target=target,
        n_total=n_total, members=sorted(selected), seeds=sorted(set(seed_names)))
    _write_manifest(manifest, out_dir)
    return manifest


def _build_manifest(subsampled, out_dir, target, n_total, members, seeds):
    """Assemble the manifest dict and write subsample_members.txt."""
    members_file = os.path.join(out_dir, "subsample_members.txt")
    with open(members_file, "w") as fh:
        for m in members:
            fh.write(m + "\n")
    return {
        "subsampled": subsampled,
        "target": target,
        "n_total": n_total,
        "n_selected": len(members),
        "members_file": members_file,
        "members": members,
        "seeds": seeds,
    }


def _write_manifest(manifest, out_dir):
    path = os.path.join(out_dir, "subsample_manifest.json")
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Pick a diverse size-bounded subset of a genome set via MASH "
                    "greedy max-min sampling (no O(n^2) matrix).")
    p.add_argument("--genome-dir", required=True,
                   help="Directory of genome FASTAs (.fna/.fasta/.fa) to subsample.")
    p.add_argument("--out-dir", required=True,
                   help="Directory for combined.msh, subsample_members.txt, and "
                        "subsample_manifest.json.")
    p.add_argument("--n", type=int, required=True, dest="target",
                   help="Target number of genomes to select.")
    p.add_argument("--must-keep", action="append", default=[],
                   help="A genome stem/basename that MUST be selected (seeds the "
                        "greedy pick); repeatable. Missing ones are ignored.")
    p.add_argument("--threads", type=int, default=1,
                   help="Threads for mash (default 1).")
    args = p.parse_args(argv)

    manifest = subsample(
        genome_dir=args.genome_dir, out_dir=args.out_dir, target=args.target,
        threads=args.threads, must_keep=args.must_keep,
    )
    print("Subsampled=%s: selected %d of %d genome(s) (target %d)"
          % (manifest["subsampled"], manifest["n_selected"],
             manifest["n_total"], manifest["target"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
