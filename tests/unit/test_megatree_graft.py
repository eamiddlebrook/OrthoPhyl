"""Tests for python_scripts/megatree_graft.py.

The graft / conflict-detection core takes ete3 Tree objects, so it is tested by
injecting small newick strings -- no mash, OrthoPhyl, or subprocess involved.
Trees are parsed with ete3 format=0 so internal-node labels are read as support
(matching what the module's CLI does).
"""

import pytest

from ete3 import Tree


@pytest.fixture
def mg(megatree_graft_module):
    return megatree_graft_module


def _tree(nw):
    return Tree(nw, format=0)


# --------------------------------------------------------------------------- #
# graft -- MRCA-replace                                                        #
# --------------------------------------------------------------------------- #

class TestGraft:
    def test_replaces_monophyletic_clade(self, mg):
        # Backbone A|B split; subclade A adds A3, subclade B adds B4.
        backbone = _tree("((A1,A2),(B1,B2));")
        subs = {"A": _tree("(A1,(A2,A3));"), "B": _tree("(B1,(B2,B4));")}
        reps = {"A": ["A1", "A2"], "B": ["B1", "B2"]}

        merged, report = mg.graft(backbone, subs, reps, min_support=90)

        assert set(merged.get_leaf_names()) == {"A1", "A2", "A3", "B1", "B2", "B4"}
        # Backbone A|B bipartition preserved: all A leaves share an ancestor that
        # excludes every B leaf.
        a_anc = merged.get_common_ancestor(["A1", "A2", "A3"])
        assert set(a_anc.get_leaf_names()) == {"A1", "A2", "A3"}
        assert all(e["monophyletic"] for e in report)
        assert all(e["conflicts"] == [] for e in report)

    def test_inputs_not_mutated(self, mg):
        backbone = _tree("((A1,A2),(B1,B2));")
        subA = _tree("(A1,(A2,A3));")
        subs = {"A": subA, "B": _tree("(B1,(B2,B4));")}
        reps = {"A": ["A1", "A2"], "B": ["B1", "B2"]}

        mg.graft(backbone, subs, reps, min_support=90)
        # Original backbone and subclade trees are untouched (graft copies them).
        assert set(backbone.get_leaf_names()) == {"A1", "A2", "B1", "B2"}
        assert set(subA.get_leaf_names()) == {"A1", "A2", "A3"}

    def test_single_rep_subclade_grafts(self, mg):
        # Subclade C anchored by a single backbone leaf.
        backbone = _tree("((A1,A2),(C1,B1));")
        subs = {"C": _tree("(C1,(C2,C3));")}
        reps = {"C": ["C1"]}

        merged, report = mg.graft(backbone, subs, reps, min_support=90)
        assert {"C1", "C2", "C3"}.issubset(set(merged.get_leaf_names()))
        # A1, A2, B1 survive untouched.
        assert {"A1", "A2", "B1"}.issubset(set(merged.get_leaf_names()))
        c = [e for e in report if e["subclade"] == "C"][0]
        assert c["monophyletic"] is True
        assert c["n_reps"] == 1

    def test_missing_reps_recorded_as_error(self, mg):
        backbone = _tree("((A1,A2),(B1,B2));")
        subs = {"Z": _tree("(Z1,(Z2,Z3));")}
        reps = {"Z": ["Q1", "Q2"]}  # none present in backbone

        merged, report = mg.graft(backbone, subs, reps, min_support=90)
        z = [e for e in report if e["subclade"] == "Z"][0]
        assert z["monophyletic"] is False
        assert "error" in z
        # Backbone unchanged, Z not grafted.
        assert "Z1" not in merged.get_leaf_names()

    def test_non_monophyletic_reps_flagged_but_leaves_survive(self, mg):
        # Reps A1 and A3 are interleaved with B leaves -> not monophyletic.
        backbone = _tree("(((A1,B1),(A3,B2)),X1);")
        subs = {"A": _tree("(A1,(A3,A9));")}
        reps = {"A": ["A1", "A3"]}

        merged, report = mg.graft(backbone, subs, reps, min_support=90)
        a = [e for e in report if e["subclade"] == "A"][0]
        assert a["monophyletic"] is False
        # Foreign leaves (B1, B2, X1) preserved; subclade leaves added.
        leaves = set(merged.get_leaf_names())
        assert {"B1", "B2", "X1", "A1", "A3", "A9"}.issubset(leaves)

    def test_deterministic(self, mg):
        backbone = _tree("((A1,A2),(B1,B2));")
        subs = {"A": _tree("(A1,(A2,A3));"), "B": _tree("(B1,(B2,B4));")}
        reps = {"A": ["A1", "A2"], "B": ["B1", "B2"]}
        m1, _ = mg.graft(backbone, subs, reps, min_support=90)
        m2, _ = mg.graft(backbone, subs, reps, min_support=90)
        assert m1.write(format=9) == m2.write(format=9)


# --------------------------------------------------------------------------- #
# induced_conflicts -- bipartition disagreement                                #
# --------------------------------------------------------------------------- #

class TestInducedConflicts:
    def test_finds_high_support_conflict(self, mg):
        # Over reps {A1..A4}: backbone groups (A1,A2) vs subclade groups (A1,A3),
        # both at support 100 -> incompatible high-support conflict.
        backbone = _tree("((A1:1,A2:1)100:1,(A3:1,A4:1)100:1);")
        subclade = _tree("((A1:1,A3:1)100:1,(A2:1,A4:1)100:1);")
        reps = ["A1", "A2", "A3", "A4"]

        conflicts = mg.induced_conflicts(backbone, subclade, reps, min_support=90)
        assert len(conflicts) >= 1
        c = conflicts[0]
        assert c["subclade_support"] >= 90
        assert c["backbone_support"] >= 90

    def test_ignores_low_support_disagreement(self, mg):
        # Same disagreement but subclade split is weakly supported (50) -> ignored.
        backbone = _tree("((A1:1,A2:1)100:1,(A3:1,A4:1)100:1);")
        subclade = _tree("((A1:1,A3:1)50:1,(A2:1,A4:1)50:1);")
        reps = ["A1", "A2", "A3", "A4"]

        conflicts = mg.induced_conflicts(backbone, subclade, reps, min_support=90)
        assert conflicts == []

    def test_agreeing_topologies_no_conflict(self, mg):
        backbone = _tree("((A1:1,A2:1)100:1,(A3:1,A4:1)100:1);")
        subclade = _tree("((A1:1,A2:1)100:1,(A3:1,A4:1)100:1);")
        reps = ["A1", "A2", "A3", "A4"]
        assert mg.induced_conflicts(backbone, subclade, reps, min_support=90) == []

    def test_too_few_reps_returns_empty(self, mg):
        backbone = _tree("((A1,A2),(A3,B1));")
        subclade = _tree("(A1,(A2,A3));")
        # Only 3 shared reps -> no non-trivial bipartition possible.
        assert mg.induced_conflicts(backbone, subclade, ["A1", "A2", "A3"],
                                    min_support=90) == []


# --------------------------------------------------------------------------- #
# helpers                                                                      #
# --------------------------------------------------------------------------- #

class TestHelpers:
    def test_leaf_labels(self, mg):
        assert mg._leaf_labels(_tree("((A1,A2),B1);")) == {"A1", "A2", "B1"}

    def test_bipartitions_skips_trivial(self, mg):
        # 4-taxon balanced tree -> exactly one non-trivial internal split.
        splits = mg._bipartitions(_tree("((A1:1,A2:1)70:1,(A3:1,A4:1)80:1);"))
        assert len(splits) == 1
        (side, sup), = splits.items()
        assert sup in (70.0, 80.0)
        assert side in (frozenset({"A1", "A2"}), frozenset({"A3", "A4"}))

    def test_compatible(self, mg):
        leaves = {"a", "b", "c", "d"}
        # Nested splits are compatible; crossing splits are not.
        assert mg._compatible({"a", "b"}, leaves, {"a", "b"})
        assert not mg._compatible({"a", "b"}, leaves, {"a", "c"})


# --------------------------------------------------------------------------- #
# CLI arg parsing                                                              #
# --------------------------------------------------------------------------- #

class TestParseSubcladeArg:
    def test_parses_spec(self, mg):
        name, path, reps = mg._parse_subclade_arg("Andr_1:/x/tree.nwk:A1,A2,A3")
        assert name == "Andr_1"
        assert path == "/x/tree.nwk"
        assert reps == ["A1", "A2", "A3"]

    def test_rejects_missing_reps(self, mg):
        import argparse
        with pytest.raises(argparse.ArgumentTypeError):
            mg._parse_subclade_arg("Andr_1:/x/tree.nwk:")

    def test_rejects_malformed(self, mg):
        import argparse
        with pytest.raises(argparse.ArgumentTypeError):
            mg._parse_subclade_arg("nocolons")
