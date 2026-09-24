"""Tests for python_scripts/subclade_partition.py (MASH subclade partitioning).

Clustering (partition_matrix) is a pure function of (names, condensed distance
array), so it is tested directly by injecting a synthetic array -- no mash
needed. partition_matrix consumes scipy's CONDENSED (1-D, upper-triangle-only)
form directly (never a dense matrix -- see condensed_index), so tests build a
dense matrix for readability (block patterns are easy to eyeball) and convert
via squareform() immediately before calling partition_matrix. The orchestration
layer (partition/find_genomes/_matches_member) is tested with the mash calls
monkeypatched.
"""

import numpy as np
import pytest
from scipy.spatial.distance import squareform


@pytest.fixture
def sp(subclade_partition_module):
    return subclade_partition_module


# --------------------------------------------------------------------------- #
# partition_matrix -- pure clustering                                          #
# --------------------------------------------------------------------------- #

class TestPartitionMatrix:
    def _two_group_matrix(self, names):
        """Build a block distance matrix: group A {0,1,2,3}, group B {4,5,6,7}.

        Within-group distances small (0.05), between-group large (0.8).
        """
        n = len(names)
        D = np.full((n, n), 0.8)
        np.fill_diagonal(D, 0.0)
        groups = [range(0, 4), range(4, 8)]
        for g in groups:
            for i in g:
                for j in g:
                    if i != j:
                        D[i, j] = 0.05
        return D

    def test_passthrough_when_under_max(self, sp):
        names = [f"g{i}.fna" for i in range(5)]
        D = np.zeros((5, 5))
        condensed = squareform(D, checks=False)
        clusters = sp.partition_matrix(names, condensed, max_size=10, min_size=4)
        assert len(clusters) == 1
        assert clusters[0] == sorted(names)

    def test_splits_two_clear_groups(self, sp):
        names = [f"g{i}.fna" for i in range(8)]
        D = self._two_group_matrix(names)
        condensed = squareform(D, checks=False)
        clusters = sp.partition_matrix(names, condensed, max_size=4, min_size=4)
        assert len(clusters) == 2
        # Each cluster is one of the two blocks.
        as_sets = sorted([tuple(sorted(c)) for c in clusters])
        assert as_sets == [
            tuple(sorted(f"g{i}.fna" for i in range(0, 4))),
            tuple(sorted(f"g{i}.fna" for i in range(4, 8))),
        ]

    def test_every_cluster_within_max(self, sp):
        names = [f"g{i}.fna" for i in range(20)]
        rng = np.random.RandomState(0)
        D = rng.uniform(0.1, 0.9, size=(20, 20))
        D = (D + D.T) / 2
        np.fill_diagonal(D, 0.0)
        condensed = squareform(D, checks=False)
        clusters = sp.partition_matrix(names, condensed, max_size=5, min_size=1)
        assert all(len(c) <= 5 for c in clusters)
        # Every genome appears exactly once.
        flat = sorted(x for c in clusters for x in c)
        assert flat == sorted(names)

    def test_tiny_cluster_merges_into_nearest(self, sp):
        # Group A {0,1,2,3,4} tight; lone outlier g5 closest to A.
        names = [f"g{i}.fna" for i in range(6)]
        n = 6
        D = np.full((n, n), 0.9)
        np.fill_diagonal(D, 0.0)
        for i in range(5):
            for j in range(5):
                if i != j:
                    D[i, j] = 0.05
        # g5 is somewhat close to group A but not part of it.
        for i in range(5):
            D[5, i] = D[i, 5] = 0.3
        condensed = squareform(D, checks=False)
        clusters = sp.partition_matrix(names, condensed, max_size=5, min_size=4)
        # The lone outlier (< min_size) must be merged, not left alone.
        assert all(len(c) >= 4 for c in clusters)
        assert sum(len(c) for c in clusters) == 6

    def test_deterministic(self, sp):
        names = [f"g{i}.fna" for i in range(8)]
        D = self._two_group_matrix(names)
        condensed = squareform(D, checks=False)
        a = sp.partition_matrix(names, condensed, max_size=4, min_size=4)
        b = sp.partition_matrix(names, condensed, max_size=4, min_size=4)
        assert a == b

    def test_ordering_size_desc(self, sp):
        # Group A has 5, group B has 3; A must come first (size desc).
        names = [f"g{i}.fna" for i in range(8)]
        n = 8
        D = np.full((n, n), 0.9)
        np.fill_diagonal(D, 0.0)
        for g in (range(0, 5), range(5, 8)):
            for i in g:
                for j in g:
                    if i != j:
                        D[i, j] = 0.05
        condensed = squareform(D, checks=False)
        clusters = sp.partition_matrix(names, condensed, max_size=5, min_size=1)
        assert len(clusters[0]) >= len(clusters[-1])


# --------------------------------------------------------------------------- #
# condensed_index -- cross-validated against scipy's own squareform            #
# --------------------------------------------------------------------------- #

class TestCondensedIndex:
    def test_matches_scipy_squareform(self, sp):
        n = 7
        rng = np.random.RandomState(0)
        D = rng.uniform(0.1, 0.9, size=(n, n))
        D = (D + D.T) / 2
        np.fill_diagonal(D, 0.0)
        condensed = squareform(D, checks=False)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                assert condensed[sp.condensed_index(i, j, n)] == pytest.approx(D[i, j])

    def test_symmetric_in_i_j(self, sp):
        assert sp.condensed_index(2, 5, 8) == sp.condensed_index(5, 2, 8)

    def test_diagonal_rejected(self, sp):
        with pytest.raises(ValueError):
            sp.condensed_index(3, 3, 8)


# --------------------------------------------------------------------------- #
# parse_mash_edges                                                             #
# --------------------------------------------------------------------------- #

class TestParseMashEdges:
    def test_symmetric_missing_default_max(self, sp, tmp_path):
        names = ["a.fna", "b.fna", "c.fna"]
        edges = tmp_path / "MASH_out"
        # Only a-b given; a-c and b-c missing -> default MASH_MAX_DIST.
        edges.write_text("/path/a.fna\t/path/b.fna\t0.10\t0\t100/1000\n")
        condensed = sp.parse_mash_edges(str(edges), names)
        n = len(names)
        assert condensed[sp.condensed_index(0, 1, n)] == pytest.approx(0.10)
        assert condensed[sp.condensed_index(0, 2, n)] == pytest.approx(sp.MASH_MAX_DIST)
        assert condensed[sp.condensed_index(1, 2, n)] == pytest.approx(sp.MASH_MAX_DIST)


# --------------------------------------------------------------------------- #
# _matches_member / _strip_ext                                                 #
# --------------------------------------------------------------------------- #

class TestMemberMatch:
    def test_matches_with_and_without_ext(self, sp):
        members = ["GCF_001.fna", "GCF_002.fna"]
        assert sp._matches_member("GCF_001", members)
        assert sp._matches_member("GCF_001.fna", members)
        assert not sp._matches_member("GCF_999", members)


# --------------------------------------------------------------------------- #
# partition orchestration (mash mocked)                                        #
# --------------------------------------------------------------------------- #

class TestPartitionOrchestration:
    def _make_genomes(self, d, n):
        d.mkdir(parents=True, exist_ok=True)
        paths = []
        for i in range(n):
            p = d / f"g{i}.fna"
            p.write_text(">c\nACGT\n")
            paths.append(p)
        return paths

    def test_passthrough_no_mash(self, sp, tmp_path, monkeypatch):
        gdir = tmp_path / "genomes"
        self._make_genomes(gdir, 3)
        out = tmp_path / "out"

        called = {"triangle": 0, "sketch": 0}
        monkeypatch.setattr(sp, "run_mash_triangle",
                            lambda *a, **k: called.__setitem__("triangle", called["triangle"] + 1))
        monkeypatch.setattr(sp, "write_sketch",
                            lambda *a, **k: called.__setitem__("sketch", called["sketch"] + 1))

        manifest = sp.partition(str(gdir), "Taxon", str(out), max_size=10,
                                min_size=4, threads=1, queries=[])
        assert manifest["partitioned"] is False
        assert manifest["n_subclades"] == 1
        assert called["triangle"] == 0
        assert called["sketch"] == 0

    def test_partition_writes_sketches_and_assigns_query(self, sp, tmp_path, monkeypatch):
        gdir = tmp_path / "genomes"
        self._make_genomes(gdir, 8)
        out = tmp_path / "out"
        out.mkdir()
        names = [f"g{i}.fna" for i in range(8)]

        # Fake mash triangle: write a two-block edge file.
        def fake_triangle(genome_files, out_path, threads):
            lines = []
            gs = [range(0, 4), range(4, 8)]
            for gi, g in enumerate(gs):
                members = list(g)
                for a in range(len(members)):
                    for b in range(a):
                        i, j = members[a], members[b]
                        lines.append(f"{gdir}/g{i}.fna\t{gdir}/g{j}.fna\t0.05\t0\t900/1000")
            # cross-block distances
            for i in range(4):
                for j in range(4, 8):
                    lines.append(f"{gdir}/g{i}.fna\t{gdir}/g{j}.fna\t0.8\t0\t10/1000")
            with open(out_path, "w") as fh:
                fh.write("\n".join(lines) + "\n")
        monkeypatch.setattr(sp, "run_mash_triangle", fake_triangle)

        sketch_args = []
        def fake_sketch(name, member_paths, out_dir, threads):
            sketch_args.append((name, list(member_paths)))
            p = f"{out_dir}/{name}.msh"
            open(p, "w").close()
            return p
        monkeypatch.setattr(sp, "write_sketch", fake_sketch)

        manifest = sp.partition(str(gdir), "Taxon", str(out), max_size=4,
                                min_size=4, threads=2, queries=["g0.fna"])
        assert manifest["partitioned"] is True
        assert manifest["n_subclades"] == 2
        # One sketch per subclade.
        assert len(sketch_args) == 2
        # Query g0 assigned to exactly one subclade.
        assert "g0.fna" in manifest["query_assignments"]

    def test_run_mash_triangle_args(self, sp, monkeypatch):
        captured = {}
        def fake_run(cmd, stdout_path=None):
            captured["cmd"] = cmd
            captured["stdout"] = stdout_path
        monkeypatch.setattr(sp, "_run_mash", fake_run)
        sp.run_mash_triangle(["a.fna", "b.fna"], "/tmp/OUT", threads=4)
        cmd = captured["cmd"]
        assert cmd[:2] == ["mash", "triangle"]
        assert "-k" in cmd and cmd[cmd.index("-k") + 1] == "17"
        assert "-s" in cmd and cmd[cmd.index("-s") + 1] == "5000"
        assert "-p" in cmd and cmd[cmd.index("-p") + 1] == "4"
        assert "-E" in cmd
        assert captured["stdout"] == "/tmp/OUT"

    def test_write_sketch_args(self, sp, monkeypatch):
        captured = {}
        monkeypatch.setattr(sp, "_run_mash", lambda cmd, **k: captured.__setitem__("cmd", cmd))
        path = sp.write_sketch("Taxon_1", ["a.fna", "b.fna"], "/tmp/o", threads=3)
        cmd = captured["cmd"]
        assert cmd[:2] == ["mash", "sketch"]
        assert cmd[cmd.index("-k") + 1] == "17"
        assert cmd[cmd.index("-s") + 1] == "5000"
        assert cmd[cmd.index("-p") + 1] == "3"
        assert path.endswith("Taxon_1.msh")
