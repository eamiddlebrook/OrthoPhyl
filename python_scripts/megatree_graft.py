#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
megatree_graft.py

Assemble ONE "megatree" from a set of per-subclade trees and a small backbone
tree, for the wrapper's opt-in --megatree large-taxon path.

A huge taxon is partitioned into size-bounded subclades (subclade_partition.py);
a tree is built for each subclade over ALL its genomes; and a small BACKBONE tree
is built from a few diverse representatives of every subclade (so each subclade is
guaranteed several anchor leaves in the backbone). This module then grafts each
subclade's full tree onto its representatives in the backbone, producing a single
tree that contains every genome while keeping the backbone's inter-subclade
topology.

Graft rule (MRCA-replace with monophyly check):
  For each subclade, find the MRCA of its backbone representatives. In the
  expected case that clade is monophyletic (its reps come from one MASH cluster),
  so the whole clade is replaced by the subclade's full tree. If the reps are NOT
  monophyletic in the backbone, we do not discard the intervening foreign leaves:
  we anchor-graft the subclade tree and prune only that subclade's own reps, and
  we flag the subclade as non-monophyletic in the conflict report.

Conflict flagging (NOT resolution):
  Independently, we compare the backbone's topology restricted to a subclade's
  reps against the subclade tree's topology restricted to the same reps, and
  record any bipartition that is strongly supported (>= --min-support) on one side
  but incompatible with a strongly-supported bipartition on the other. These are
  written to a JSON report for a future reconciliation pass; this module does NOT
  attempt to resolve them.

Only ete3 is required (already an OrthoPhyl dependency). The graft/conflict logic
is separated from any file/subprocess I/O (see graft / induced_conflicts) so it is
directly unit-testable by injecting small newick strings.

Usage:
    megatree_graft.py --backbone BACKBONE.nwk \\
        --subclade NAME:SUBCLADE.nwk:REP1,REP2,REP3 \\
        [--subclade ...] \\
        --out-tree MERGED.nwk --out-report conflicts.json [--min-support 90]
"""
import argparse
import json
import sys

from ete3 import Tree


def _leaf_labels(tree):
    """Return the set of leaf names of an ete3 tree."""
    return set(tree.get_leaf_names())


def _canon_split(side, leaves):
    """Canonical key for a bipartition: the frozenset of the lexically-smaller
    side (ties broken by min element), so the same split hashes identically no
    matter which side a tree happened to hang it from."""
    side = frozenset(side)
    other = frozenset(leaves) - side
    if len(side) < len(other):
        return side
    if len(other) < len(side):
        return other
    # Equal-size sides: pick the one whose min element sorts first.
    return side if min(side) <= min(other) else other


def _bipartitions(tree):
    """Map every non-trivial internal edge of `tree` to its support.

    Returns {canonical_split (frozenset of one side): support}. Trivial splits
    (a single leaf, or all-but-one) are skipped -- they carry no topological
    information to conflict over.
    """
    leaves = _leaf_labels(tree)
    splits = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        side = frozenset(node.get_leaf_names())
        other = leaves - side
        if len(side) < 2 or len(other) < 2:
            continue
        splits[_canon_split(side, leaves)] = float(node.support)
    return splits


def _compatible(a, leaves, b):
    """True if bipartitions `a|~a` and `b|~b` (over leaf set `leaves`) can coexist
    on one tree. Two splits are compatible iff at least one of the four
    side-intersections is empty."""
    leaves = set(leaves)
    a = set(a)
    na = leaves - a
    b = set(b)
    nb = leaves - b
    return not (a & b and a & nb and na & b and na & nb)


def induced_conflicts(backbone_tree, subclade_tree, rep_labels, min_support):
    """Flag strongly-supported topological disagreements between the backbone and
    a subclade tree, restricted to that subclade's representative leaves.

    Both trees are pruned to `rep_labels` (their shared leaf set), then each
    strongly-supported (>= min_support) bipartition of the subclade tree is
    checked for compatibility with each strongly-supported bipartition of the
    backbone. Incompatible high-support pairs are returned as conflicts.

    Pure function of the two trees + reps -- no I/O.

    Returns a list of dicts:
        {"subclade_split": [...], "backbone_split": [...],
         "subclade_support": float, "backbone_support": float}
    """
    reps = [r for r in rep_labels if r in _leaf_labels(backbone_tree)
            and r in _leaf_labels(subclade_tree)]
    if len(reps) < 4:
        # Need at least 4 leaves for a non-trivial bipartition on each side.
        return []

    bb = backbone_tree.copy()
    sc = subclade_tree.copy()
    bb.prune(reps)
    sc.prune(reps)

    leaves = set(reps)
    bb_splits = {s: sup for s, sup in _bipartitions(bb).items() if sup >= min_support}
    sc_splits = {s: sup for s, sup in _bipartitions(sc).items() if sup >= min_support}

    conflicts = []
    for s_split, s_sup in sorted(sc_splits.items(), key=lambda kv: sorted(kv[0])):
        for b_split, b_sup in sorted(bb_splits.items(), key=lambda kv: sorted(kv[0])):
            if not _compatible(s_split, leaves, b_split):
                conflicts.append({
                    "subclade_split": sorted(s_split),
                    "backbone_split": sorted(b_split),
                    "subclade_support": s_sup,
                    "backbone_support": b_sup,
                })
    return conflicts


def _find_leaf(tree, name):
    hits = tree.search_nodes(name=name)
    return hits[0] if hits else None


def graft(backbone_tree, subclade_trees, rep_map, min_support=90):
    """Graft each subclade's full tree onto its representatives in the backbone.

    Pure-ish transform (mutates a copy of backbone_tree, never the inputs):

      backbone_tree : ete3.Tree whose leaves are subclade representatives.
      subclade_trees: {subclade_name: ete3.Tree} full per-subclade trees.
      rep_map       : {subclade_name: [rep_leaf_label, ...]} that subclade's
                      backbone anchors.
      min_support   : support threshold for conflict flagging.

    For each subclade (processed in sorted name order for determinism):
      * locate its reps' MRCA in the (current) backbone,
      * if the reps are monophyletic, replace that whole clade with the subclade
        tree; otherwise anchor-graft the subclade tree and prune only that
        subclade's reps (preserving any foreign leaves) and flag it,
      * record induced high-support bipartition conflicts.

    Returns (merged_tree, report) where report is a list of per-subclade dicts:
        {"subclade", "monophyletic", "n_reps", "conflicts": [...]}
    """
    merged = backbone_tree.copy()
    report = []

    for name in sorted(subclade_trees):
        sc_tree = subclade_trees[name]
        reps = rep_map.get(name, [])
        present = [r for r in reps if _find_leaf(merged, r) is not None]

        entry = {"subclade": name, "n_reps": len(present),
                 "monophyletic": None, "conflicts": []}

        if not present:
            # No anchor to graft onto -- cannot place this subclade.
            entry["monophyletic"] = False
            entry["error"] = "no representatives found in backbone"
            report.append(entry)
            continue

        # Conflict flagging is independent of the graft mechanics.
        entry["conflicts"] = induced_conflicts(
            merged, sc_tree, present, min_support)

        graft_subtree = sc_tree.copy()

        if len(present) == 1:
            leaf = _find_leaf(merged, present[0])
            entry["monophyletic"] = True  # a single leaf is trivially monophyletic
            parent = leaf.up
            graft_subtree.dist = leaf.dist
            if parent is None:
                merged = graft_subtree
            else:
                parent.add_child(graft_subtree)
                leaf.detach()
            report.append(entry)
            continue

        is_mono, _, _ = merged.check_monophyly(values=present, target_attr="name")
        entry["monophyletic"] = bool(is_mono)
        mrca = merged.get_common_ancestor(present)

        if is_mono and mrca.up is not None:
            # Clean replacement: the reps' clade is exactly this subclade.
            graft_subtree.dist = mrca.dist
            mrca.up.add_child(graft_subtree)
            mrca.detach()
        else:
            # Non-monophyletic (or reps span the whole backbone): anchor-graft and
            # remove only this subclade's own reps so foreign leaves survive.
            anchor = _find_leaf(merged, present[0])
            parent = anchor.up or mrca
            parent.add_child(graft_subtree)
            for r in present:
                node = _find_leaf(merged, r)
                if node is not None and node.up is not None:
                    node.detach()

        report.append(entry)

    _collapse_unifurcations(merged)
    return merged, report


def _collapse_unifurcations(tree):
    """Remove internal nodes left with a single child after grafting/pruning."""
    for node in list(tree.traverse()):
        if node.is_root():
            continue
        if not node.is_leaf() and len(node.children) == 1:
            node.delete()  # reconnect the lone child to node's parent


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #

def _parse_subclade_arg(spec):
    """Parse a --subclade NAME:TREE_PATH:REP1,REP2,... spec into
    (name, tree_path, [reps]). The name and path must not contain ':'."""
    parts = spec.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "expected NAME:TREE_PATH:REP1,REP2,... got %r" % spec)
    name, tree_path, reps = parts
    rep_list = [r for r in reps.split(",") if r]
    if not rep_list:
        raise argparse.ArgumentTypeError(
            "subclade %s has no representatives" % name)
    return name, tree_path, rep_list


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Graft per-subclade trees onto a backbone into one megatree.")
    p.add_argument("--backbone", required=True,
                   help="Backbone newick (leaves are subclade representatives).")
    p.add_argument("--subclade", action="append", default=[], required=True,
                   type=_parse_subclade_arg, metavar="NAME:TREE:REP,REP,...",
                   help="A subclade's full tree and its backbone representatives; "
                        "repeatable.")
    p.add_argument("--out-tree", required=True,
                   help="Path to write the merged megatree (newick).")
    p.add_argument("--out-report", required=True,
                   help="Path to write the JSON conflict report.")
    p.add_argument("--min-support", type=float, default=90,
                   help="Support threshold for flagging bipartition conflicts "
                        "(assumes a 0-100 support scale, e.g. IQ-TREE UFBoot). "
                        "Default 90.")
    args = p.parse_args(argv)

    # format=0: flexible parse that reads internal-node labels as support
    # (IQ-TREE UFBoot / FastTree support), defaulting to 1.0 when absent.
    backbone = Tree(args.backbone, format=0)
    subclade_trees = {}
    rep_map = {}
    for name, tree_path, reps in args.subclade:
        subclade_trees[name] = Tree(tree_path, format=0)
        rep_map[name] = reps

    merged, report = graft(backbone, subclade_trees, rep_map,
                           min_support=args.min_support)

    merged.write(outfile=args.out_tree, format=0)
    with open(args.out_report, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)

    n_conflicts = sum(len(e["conflicts"]) for e in report)
    n_nonmono = sum(1 for e in report if e.get("monophyletic") is False)
    print("Grafted %d subclade(s) onto backbone; %d non-monophyletic, "
          "%d high-support bipartition conflict(s) flagged."
          % (len(report), n_nonmono, n_conflicts))
    return 0


if __name__ == "__main__":
    sys.exit(main())
