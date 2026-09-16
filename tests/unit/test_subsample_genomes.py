"""Tests for python_scripts/subsample_genomes.py (MASH greedy max-min subsample).

The selection (greedy_maxmin) is a pure function of (names, distance oracle), so
it is tested directly by injecting distances -- no mash needed. The orchestration
layer (subsample/find_genomes/manifest) is tested with the mash calls mocked.
"""

import json
import math
import os

import pytest


@pytest.fixture
def ss(subsample_module):
    return subsample_module


# --------------------------------------------------------------------------- #
# greedy_maxmin -- pure farthest-point selection                              #
# --------------------------------------------------------------------------- #

class TestGreedyMaxMin:
    def _oracle(self, dist):
        """Build a dist_fn from a {(a,b): d} dict (symmetric); missing -> inf."""
        def dist_fn(name):
            row = {}
            for (a, b), d in dist.items():
                if a == name:
                    row[b] = d
                elif b == name:
                    row[a] = d
            return row
        return dist_fn

    def test_returns_all_when_target_ge_n(self, ss):
        names = ["a", "b", "c"]
        got = ss.greedy_maxmin(names, self._oracle({}), target=5)
        assert got == ["a", "b", "c"]

    def test_empty_when_target_zero(self, ss):
        assert ss.greedy_maxmin(["a", "b"], self._oracle({}), target=0) == []

    def test_picks_farthest_first(self, ss):
        # a-b very close, a-c and b-c far. Seeded at 'a' (lexical first), the next
        # pick must be the farthest genome = c, not b.
        names = ["a", "b", "c"]
        dist = {("a", "b"): 0.01, ("a", "c"): 0.9, ("b", "c"): 0.9}
        got = ss.greedy_maxmin(names, self._oracle(dist), target=2)
        assert got == ["a", "c"]

    def test_honors_seeds(self, ss):
        names = ["a", "b", "c", "d"]
        dist = {
            ("a", "b"): 0.5, ("a", "c"): 0.5, ("a", "d"): 0.5,
            ("b", "c"): 0.5, ("b", "d"): 0.5, ("c", "d"): 0.5,
        }
        got = ss.greedy_maxmin(names, self._oracle(dist), target=2, seed_names=["d"])
        assert got[0] == "d"
        assert len(got) == 2

    def test_ignores_unknown_seeds(self, ss):
        names = ["a", "b"]
        got = ss.greedy_maxmin(names, self._oracle({("a", "b"): 0.5}),
                               target=1, seed_names=["zzz"])
        # Unknown seed dropped -> falls back to lexical-first seed 'a'.
        assert got == ["a"]

    def test_deterministic_tiebreak_by_name(self, ss):
        # All equidistant: after seed 'a', b/c/d tie -> lexical order b then c.
        names = ["a", "b", "c", "d"]
        dist = {
            ("a", "b"): 0.5, ("a", "c"): 0.5, ("a", "d"): 0.5,
            ("b", "c"): 0.5, ("b", "d"): 0.5, ("c", "d"): 0.5,
        }
        a = ss.greedy_maxmin(names, self._oracle(dist), target=3)
        b = ss.greedy_maxmin(names, self._oracle(dist), target=3)
        assert a == b
        assert a[0] == "a"

    def test_selection_count(self, ss):
        names = [f"g{i}" for i in range(10)]
        # Random-ish but fixed distances via a simple formula.
        dist = {}
        for i in range(10):
            for j in range(i):
                dist[(f"g{i}", f"g{j}")] = 0.1 + 0.01 * (i + j)
        got = ss.greedy_maxmin(names, self._oracle(dist), target=4)
        assert len(got) == 4
        assert len(set(got)) == 4


# --------------------------------------------------------------------------- #
# subsample orchestration (mash mocked)                                        #
# --------------------------------------------------------------------------- #

class TestSubsampleOrchestration:
    def _make_genomes(self, d, n):
        d.mkdir(parents=True, exist_ok=True)
        for i in range(n):
            (d / f"g{i}.fna").write_text(">c\nACGT\n")

    def test_passthrough_when_under_target(self, ss, tmp_path, monkeypatch):
        gdir = tmp_path / "genomes"
        self._make_genomes(gdir, 3)
        out = tmp_path / "out"

        called = {"n": 0}
        monkeypatch.setattr(ss, "_run_mash",
                            lambda *a, **k: called.__setitem__("n", called["n"] + 1))
        manifest = ss.subsample(str(gdir), str(out), target=10, threads=1)
        assert manifest["subsampled"] is False
        assert manifest["n_selected"] == 3
        assert manifest["n_total"] == 3
        assert called["n"] == 0  # no sketch when nothing to subsample

    def test_subsamples_and_writes_manifest(self, ss, tmp_path, monkeypatch):
        gdir = tmp_path / "genomes"
        self._make_genomes(gdir, 6)
        out = tmp_path / "out"

        # mash sketch is a no-op; mash dist returns a fixed distance to every ref.
        monkeypatch.setattr(ss, "_run_mash", lambda *a, **k: None)

        def fake_row(query_path, combined_msh):
            # Distance = 0.1 to everything (ties -> deterministic lexical picks).
            row = {}
            for i in range(6):
                row[f"g{i}.fna"] = 0.1
            return row
        monkeypatch.setattr(ss, "_mash_dist_row", fake_row)

        manifest = ss.subsample(str(gdir), str(out), target=3, threads=1)
        assert manifest["subsampled"] is True
        assert manifest["n_selected"] == 3
        assert manifest["n_total"] == 6
        # Members file + manifest written.
        assert (out / "subsample_members.txt").exists()
        assert (out / "subsample_manifest.json").exists()
        on_disk = json.loads((out / "subsample_manifest.json").read_text())
        assert on_disk["n_selected"] == 3

    def test_must_keep_seeds_retained(self, ss, tmp_path, monkeypatch):
        gdir = tmp_path / "genomes"
        self._make_genomes(gdir, 6)
        out = tmp_path / "out"
        monkeypatch.setattr(ss, "_run_mash", lambda *a, **k: None)
        monkeypatch.setattr(ss, "_mash_dist_row",
                            lambda q, m: {f"g{i}.fna": 0.1 for i in range(6)})

        # Seed with g5 by bare accession (no extension) -- must survive.
        manifest = ss.subsample(str(gdir), str(out), target=2, threads=1,
                                must_keep=["g5"])
        assert "g5.fna" in manifest["members"]
        assert "g5.fna" in manifest["seeds"]

    def test_no_genomes_raises(self, ss, tmp_path):
        gdir = tmp_path / "empty"
        gdir.mkdir()
        with pytest.raises(SystemExit):
            ss.subsample(str(gdir), str(tmp_path / "out"), target=5, threads=1)

    def test_sketch_uses_file_of_filenames_not_bare_argv(self, ss, tmp_path, monkeypatch):
        # mash sketch must take -l <filelist>, not one argv entry per genome --
        # bare argv overflows execve's ARG_MAX at tens of thousands of genomes.
        gdir = tmp_path / "genomes"
        self._make_genomes(gdir, 6)
        out = tmp_path / "out"

        seen_cmd = {}
        def fake_run_mash(cmd, stdout_path=None):
            seen_cmd["cmd"] = cmd
        monkeypatch.setattr(ss, "_run_mash", fake_run_mash)
        monkeypatch.setattr(ss, "_mash_dist_row",
                            lambda q, m: {f"g{i}.fna": 0.1 for i in range(6)})

        ss.subsample(str(gdir), str(out), target=3, threads=1)

        cmd = seen_cmd["cmd"]
        assert "-l" in cmd
        filelist = cmd[cmd.index("-l") + 1]
        assert os.path.exists(filelist)
        listed = set(open(filelist).read().split())
        assert len(listed) == 6
        # None of the genome paths should appear as bare positional argv entries.
        for i in range(6):
            assert str(gdir / f"g{i}.fna") not in cmd
