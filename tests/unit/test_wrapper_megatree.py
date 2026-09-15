"""Tests for the opt-in --megatree large-taxon path in the wrapper.

The megatree path (only when self.megatree AND raw_count > max_tree_genomes)
is: enforce the total-genome ceiling -> partition into size-bounded subclades ->
build a full tree per subclade -> pick backbone reps per subclade -> build one
backbone tree -> shell to megatree_graft.py -> create the taxon DB from the
backbone run. Every heavy step (partition, subsample, build, OrthoPhyl, grafter
subprocess, DB creation) is mocked; these tests assert the wiring/dispatch, not
the bioinformatics.
"""

from pathlib import Path

import pytest


@pytest.fixture
def Wrapper(wrapper_module):
    return wrapper_module.PipelineWrapper


def _make_wrapper(Wrapper, tmp_path, **overrides):
    db_dir = tmp_path / "databases"
    db_dir.mkdir(exist_ok=True)
    fake_gather = tmp_path / "fake_gather.sh"
    fake_gather.write_text("#!/bin/bash\necho fake\n")
    fake_gather.chmod(0o755)
    kwargs = dict(
        input_file=tmp_path / "assemblies.tsv",
        database_dir=db_dir,
        output_dir=tmp_path / "out",
        threads=4,
        gather_script=fake_gather,
        max_tree_genomes=5,
    )
    kwargs.update(overrides)
    w = Wrapper(**kwargs)
    w.logs_dir.mkdir(parents=True, exist_ok=True)
    w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
    return w


def _query(tmp_path, stem="GCF_query"):
    p = tmp_path / f"{stem}.fna"
    p.write_text(">c\nACGT\n")
    return {"assembly_id": stem, "assembly_path": str(p),
            "download_taxonomy": "d__Bacteria;g__Andreesenella"}


def _stub_download(w, monkeypatch, n):
    def fake_dl_raw(taxon, output_dir, query_assemblies=None):
        d = output_dir / "assemblies_all.TMP"
        d.mkdir(parents=True, exist_ok=True)
        for i in range(n):
            (d / f"g{i}.fna").write_text(">c\nAC\n")
        return d
    monkeypatch.setattr(w, "_download_raw", fake_dl_raw)


# --------------------------------------------------------------------------- #
# Dispatch: --megatree routes to _run_megatree, default subsamples             #
# --------------------------------------------------------------------------- #

class TestMegatreeDispatch:
    def test_over_ceiling_with_megatree_calls_run_megatree(
            self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True)
        _stub_download(w, monkeypatch, 8)  # 8 > ceiling (5)

        mega_calls = []
        monkeypatch.setattr(w, "_run_megatree",
                            lambda **k: mega_calls.append(k))
        monkeypatch.setattr(w, "_subsample_genomes",
                            lambda *a, **k: pytest.fail("megatree must not subsample"))
        monkeypatch.setattr(w, "_build_subclade",
                            lambda **k: pytest.fail("dispatch must go via _run_megatree"))

        w._process_orthophyl_taxon("Andreesenella", [_query(tmp_path)])
        assert len(mega_calls) == 1
        assert mega_calls[0]["taxon_name"] == "Andreesenella"

    def test_under_ceiling_with_megatree_single_build(
            self, Wrapper, tmp_path, monkeypatch):
        # Even with --megatree, a set under the ceiling is a plain single tree.
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=10, megatree=True)
        _stub_download(w, monkeypatch, 3)
        monkeypatch.setattr(w, "_run_megatree",
                            lambda **k: pytest.fail("under ceiling must not megatree"))
        built = []
        monkeypatch.setattr(w, "_build_subclade", lambda **k: built.append(k))

        w._process_orthophyl_taxon("Andreesenella", [_query(tmp_path)])
        assert len(built) == 1
        assert built[0]["is_subclade"] is False

    def test_over_ceiling_without_megatree_subsamples(
            self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=False)
        _stub_download(w, monkeypatch, 8)
        monkeypatch.setattr(w, "_run_megatree",
                            lambda **k: pytest.fail("default path must not megatree"))
        sub_calls = []
        def fake_subsample(taxon, raw_dir, target, must_keep_stems=None):
            sub_calls.append(target)
            sel = w.orthophyl_dir / "subsample" / taxon / "selected"
            sel.mkdir(parents=True, exist_ok=True)
            (sel / "g0.fna").write_text(">c\nAC\n")
            return sel
        monkeypatch.setattr(w, "_subsample_genomes", fake_subsample)
        monkeypatch.setattr(w, "_build_subclade", lambda **k: None)

        w._process_orthophyl_taxon("Andreesenella", [_query(tmp_path)])
        assert len(sub_calls) == 1


class TestCreateModeMegatreeDispatch:
    def _stub_gatherer(self, monkeypatch):
        import sys
        import types
        class FakeGatherer:
            taxid = "123"
            taxon_rank = "genus"
            def __init__(self, **k): pass
            def get_taxonomy_string(self): return "d__Bacteria;g__Andreesenella"
        fake_mod = types.ModuleType("taxon_assembly_gatherer")
        fake_mod.TaxonAssemblyGatherer = FakeGatherer
        monkeypatch.setitem(sys.modules, "taxon_assembly_gatherer", fake_mod)

    def test_create_mode_over_ceiling_megatree(
            self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, input_file=None, taxon="Andreesenella",
                          max_tree_genomes=5, megatree=True)
        self._stub_gatherer(monkeypatch)
        _stub_download(w, monkeypatch, 8)

        mega_calls = []
        monkeypatch.setattr(w, "_run_megatree", lambda **k: mega_calls.append(k))
        monkeypatch.setattr(w, "_subsample_genomes",
                            lambda *a, **k: pytest.fail("megatree must not subsample"))
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        rc = w._run_taxon_create_mode()
        assert rc == 0
        assert len(mega_calls) == 1
        # Create mode has no query.
        assert mega_calls[0]["query_assemblies"] == []
        assert mega_calls[0]["taxon_name"] == "Andreesenella"


# --------------------------------------------------------------------------- #
# _run_megatree internals: ceiling -> partition -> per-subclade build ->       #
#   backbone build -> grafter subprocess -> DB entry                           #
# --------------------------------------------------------------------------- #

class TestRunMegatree:
    def _wire(self, w, wrapper_module, monkeypatch, n_subclades=2):
        """Stub every heavy dependency of _run_megatree and return recorders."""
        rec = {"ceiling": [], "partition": [], "build": [], "subsample": [],
               "orthophyl": [], "graft_cmd": [], "db": []}

        monkeypatch.setattr(w, "_enforce_total_genome_ceiling",
                            lambda taxon, n: rec["ceiling"].append((taxon, n)))

        def fake_partition(taxon, raw_dir, queries, max_size=None):
            rec["partition"].append({"taxon": taxon, "max_size": max_size})
            subs = [{"subclade_id": i + 1, "name": f"{taxon}_{i + 1}",
                     "n_genomes": 3, "members_file": None, "sketch_file": None}
                    for i in range(n_subclades)]
            assignments = {}
            for asm in queries:
                assignments[Path(asm["assembly_path"]).name] = subs[0]["name"]
            return {"partitioned": True, "subclades": subs,
                    "query_assignments": assignments}
        monkeypatch.setattr(w, "_partition_genomes", fake_partition)

        monkeypatch.setattr(w, "_build_subclade",
                            lambda **k: rec["build"].append(k))

        def fake_subsample(taxon, raw_dir, target, must_keep_stems=None):
            rec["subsample"].append({"taxon": taxon, "target": target,
                                     "seeds": must_keep_stems})
            sel = w.orthophyl_dir / "sub" / taxon
            sel.mkdir(parents=True, exist_ok=True)
            (sel / f"{taxon}_rep.fna").write_text(">c\nAC\n")
            return sel
        monkeypatch.setattr(w, "_subsample_genomes", fake_subsample)

        monkeypatch.setattr(w, "_run_orthophyl",
                            lambda **k: rec["orthophyl"].append(k))
        monkeypatch.setattr(w, "_locate_species_tree",
                            lambda out: out / "tree.nwk")
        monkeypatch.setattr(w, "_create_database_entry",
                            lambda **k: rec["db"].append(k))

        class FakeCompleted:
            returncode = 0
        def fake_run(cmd, *a, **k):
            rec["graft_cmd"].append([str(x) for x in cmd])
            return FakeCompleted()
        monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
        return rec

    def test_full_wiring(self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True,
                          subclade_size=150, backbone_reps=5)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=2)

        raw = w.orthophyl_dir / "downloads" / "Andreesenella" / "assemblies_all.TMP"
        raw.mkdir(parents=True, exist_ok=True)
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")
        # subclade genomes_to_keep dirs so backbone rep collection has something.
        for sc in ("Andreesenella_1", "Andreesenella_2"):
            gtk = w.orthophyl_dir / "downloads" / sc / "genomes_to_keep"
            gtk.mkdir(parents=True, exist_ok=True)
            (gtk / f"{sc}_a.fna").write_text(">c\nAC\n")

        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[_query(tmp_path)],
                        taxonomy="d__Bacteria;g__Andreesenella")

        # (1) ceiling enforced with the raw count.
        assert rec["ceiling"] == [("Andreesenella", 8)]
        # (2) partitioned with --max-size = subclade_size.
        assert len(rec["partition"]) == 1
        assert rec["partition"][0]["max_size"] == 150
        # (3) one full-tree build per subclade, all is_subclade=True.
        assert len(rec["build"]) == 2
        assert all(b["is_subclade"] is True for b in rec["build"])
        # (4) backbone reps picked per subclade (target = backbone_reps).
        assert len(rec["subsample"]) == 2
        assert all(s["target"] == 5 for s in rec["subsample"])
        # (5) exactly one backbone OrthoPhyl run.
        assert len(rec["orthophyl"]) == 1
        assert rec["orthophyl"][0]["taxon_name"] == "Andreesenella_backbone"
        # (6) grafter invoked and taxon DB created from the backbone run.
        assert len(rec["graft_cmd"]) == 1
        graft = rec["graft_cmd"][0]
        assert str(w.megatree_grafter) in graft
        assert "--backbone" in graft
        assert graft.count("--subclade") == 2  # one per subclade with reps
        assert len(rec["db"]) == 1
        assert rec["db"][0]["taxon_name"] == "Andreesenella"

    def test_backbone_reps_threads_through(self, Wrapper, wrapper_module,
                                           tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True,
                          backbone_reps=3)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=1)

        raw = tmp_path / "raw"
        raw.mkdir()
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")
        gtk = w.orthophyl_dir / "downloads" / "Andreesenella_1" / "genomes_to_keep"
        gtk.mkdir(parents=True, exist_ok=True)
        (gtk / "a.fna").write_text(">c\nAC\n")

        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[], taxonomy="d__Bacteria;g__Andreesenella")
        assert all(s["target"] == 3 for s in rec["subsample"])
