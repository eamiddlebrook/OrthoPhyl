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
        # The backbone DB entry is marked is_backbone so the router can tell
        # it apart from a dense subclade sharing the same parent taxonomy.
        assert rec["db"][0]["subclade_meta"]["is_backbone"] is True

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

    def test_megatree_lazy_registers_subclades_with_no_query(
            self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        """--megatree-lazy: a subclade with NO query is only registered
        (built=false), not built; a subclade WITH a query is still built
        eagerly."""
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True,
                          megatree_lazy=True)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=2)
        registered = []
        monkeypatch.setattr(w, "_register_lazy_subclade",
                            lambda **k: registered.append(k))

        raw = w.orthophyl_dir / "downloads" / "Andreesenella" / "assemblies_all.TMP"
        raw.mkdir(parents=True, exist_ok=True)
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")
        gtk = w.orthophyl_dir / "downloads" / "Andreesenella_1" / "genomes_to_keep"
        gtk.mkdir(parents=True, exist_ok=True)
        (gtk / "a.fna").write_text(">c\nAC\n")

        # fake_partition assigns the query to subclade "Andreesenella_1" only,
        # so "Andreesenella_2" has no query and must be lazily registered.
        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[_query(tmp_path)],
                        taxonomy="d__Bacteria;g__Andreesenella")

        assert len(rec["build"]) == 1
        assert rec["build"][0]["entry"]["name"] == "Andreesenella_1"
        assert len(registered) == 1
        assert registered[0]["entry"]["name"] == "Andreesenella_2"

    def test_megatree_without_lazy_builds_every_subclade(
            self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        """Without --megatree-lazy (default), every subclade is built
        eagerly regardless of whether it holds a query."""
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=2)
        registered = []
        monkeypatch.setattr(w, "_register_lazy_subclade",
                            lambda **k: registered.append(k))

        raw = w.orthophyl_dir / "downloads" / "Andreesenella" / "assemblies_all.TMP"
        raw.mkdir(parents=True, exist_ok=True)
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")
        for sc in ("Andreesenella_1", "Andreesenella_2"):
            gtk = w.orthophyl_dir / "downloads" / sc / "genomes_to_keep"
            gtk.mkdir(parents=True, exist_ok=True)
            (gtk / f"{sc}_a.fna").write_text(">c\nAC\n")

        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[_query(tmp_path)],
                        taxonomy="d__Bacteria;g__Andreesenella")

        assert len(rec["build"]) == 2
        assert registered == []


class TestMegatreeHmmReuse:
    """--megatree-hmm-reuse: backbone built FIRST (from reps picked off each
    subclade's QC'd member pool), then every subclade is built passing
    --hmm-assign-dir pointed at the backbone's hmms_final/, instead of each
    subclade running its own independent OrthoFinder orthogroup inference.
    Default (flag off) order/behavior is covered by TestRunMegatree above and
    must stay byte-identical; these tests only exercise the opt-in path.
    """

    def _wire(self, w, wrapper_module, monkeypatch, n_subclades=2,
             hmms_final=True):
        """Stub every heavy dependency; hmms_final controls whether the fake
        backbone run's _run_orthophyl stub actually creates hmms_final/*.hmm
        (so the fallback-when-missing path can be tested too)."""
        rec = {"ceiling": [], "partition": [], "qc": [], "build": [],
               "subsample": [], "orthophyl": [], "graft_cmd": [], "db": [],
               "register": []}

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

        def fake_qc(subclade_dir, raw_members, taxon_label, query_assemblies=None):
            rec["qc"].append(taxon_label)
            gtk = subclade_dir / "genomes_to_keep"
            gtk.mkdir(parents=True, exist_ok=True)
            (gtk / f"{taxon_label}_a.fna").write_text(">c\nAC\n")
            return gtk
        monkeypatch.setattr(w, "_qc_subclade", fake_qc)

        monkeypatch.setattr(w, "_build_subclade",
                            lambda **k: rec["build"].append(k))
        monkeypatch.setattr(w, "_register_lazy_subclade",
                            lambda **k: rec["register"].append(k))

        def fake_subsample(taxon, raw_dir, target, must_keep_stems=None):
            rec["subsample"].append({"taxon": taxon, "target": target,
                                     "seeds": must_keep_stems})
            sel = w.orthophyl_dir / "subsample" / taxon / "selected"
            sel.mkdir(parents=True, exist_ok=True)
            (sel / f"{taxon}_rep.fna").write_text(">c\nAC\n")
            return sel
        monkeypatch.setattr(w, "_subsample_genomes", fake_subsample)

        def fake_run_orthophyl(input_dir, output_dir, taxon_name, assemblies,
                               hmm_assign_dir=None):
            rec["orthophyl"].append({"taxon_name": taxon_name,
                                     "hmm_assign_dir": hmm_assign_dir})
            if taxon_name.endswith("_backbone") and hmms_final:
                hmm_dir = output_dir / "OG_alignmentsToHMM" / "hmms_final"
                hmm_dir.mkdir(parents=True, exist_ok=True)
                (hmm_dir / "OG0000001.hmm").write_text("HMM\n")
        monkeypatch.setattr(w, "_run_orthophyl", fake_run_orthophyl)
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

    def test_backbone_built_before_subclades(
            self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True,
                          megatree_hmm_reuse=True)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=2)

        raw = tmp_path / "raw"
        raw.mkdir()
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")

        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[_query(tmp_path)],
                        taxonomy="d__Bacteria;g__Andreesenella")

        # QC ran for both subclades BEFORE any _build_subclade call -- the
        # backbone's OrthoPhyl run must be the first entry in rec["orthophyl"].
        assert len(rec["qc"]) == 2
        assert rec["orthophyl"][0]["taxon_name"] == "Andreesenella_backbone"
        # Every subclade build received the backbone's hmms_final/ dir.
        assert len(rec["build"]) == 2
        for b in rec["build"]:
            assert b["hmm_assign_dir"] is not None
            assert b["hmm_assign_dir"].name == "hmms_final"
            # genomes_to_keep was passed through (no redundant QC in _build_subclade).
            assert b["genomes_to_keep"] is not None

    def test_falls_back_to_independent_when_no_hmms_final(
            self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        """If the backbone run never produces hmms_final/ (e.g. too few
        pooled reps to cross --ani-shortlist), subclades still build, just
        without hmm_assign_dir -- not a hard failure."""
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True,
                          megatree_hmm_reuse=True)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=2,
                         hmms_final=False)

        raw = tmp_path / "raw"
        raw.mkdir()
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")

        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[_query(tmp_path)],
                        taxonomy="d__Bacteria;g__Andreesenella")

        assert len(rec["build"]) == 2
        for b in rec["build"]:
            assert b["hmm_assign_dir"] is None

    def test_megatree_lazy_still_registers_with_hmm_reuse(
            self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        """--megatree-lazy + --megatree-hmm-reuse together: a subclade with
        no query is still only registered, contributing no backbone reps and
        never QC'd."""
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True,
                          megatree_hmm_reuse=True, megatree_lazy=True)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=2)

        raw = tmp_path / "raw"
        raw.mkdir()
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")

        # fake_partition assigns the query to subclade "_1" only, so "_2" has
        # no query.
        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[_query(tmp_path)],
                        taxonomy="d__Bacteria;g__Andreesenella")

        assert len(rec["register"]) == 1
        assert rec["register"][0]["entry"]["name"] == "Andreesenella_2"
        assert len(rec["build"]) == 1
        assert rec["build"][0]["entry"]["name"] == "Andreesenella_1"
        # Only the built subclade was QC'd -- the lazy one never ran QC.
        assert rec["qc"] == ["Andreesenella_1"]

    def test_default_order_independent_when_flag_off(
            self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        """Regression guard: with --megatree-hmm-reuse NOT set, no QC call
        happens before the FIRST subclade build (today's independent-build
        order), confirming the reorder is fully opt-in."""
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, megatree=True)
        rec = self._wire(w, wrapper_module, monkeypatch, n_subclades=2)

        raw = tmp_path / "raw"
        raw.mkdir()
        for i in range(8):
            (raw / f"g{i}.fna").write_text(">c\nAC\n")

        w._run_megatree(taxon_name="Andreesenella", raw_dir=raw,
                        query_assemblies=[_query(tmp_path)],
                        taxonomy="d__Bacteria;g__Andreesenella")

        # Default path never calls _qc_subclade directly (it's inside the
        # real _build_subclade, which is stubbed here), and the backbone
        # OrthoPhyl run is LAST, not first.
        assert rec["qc"] == []
        assert rec["orthophyl"][-1]["taxon_name"] == "Andreesenella_backbone"
        for b in rec["build"]:
            assert b.get("hmm_assign_dir") is None


class TestSubcladeBuildPhase:
    """Phase 3c: build a lazily-registered subclade on demand, then ReLeaf
    the waiting queries onto it."""

    def test_process_subclade_build_uses_subclade_taxonomy_not_query(
            self, Wrapper, tmp_path, monkeypatch):
        """Hazard 6: rebuilding a lazily-registered subclade must use the
        subclade's OWN recorded taxonomy (subclade_taxonomy), not whatever
        taxonomy the routing query happened to carry."""
        w = _make_wrapper(Wrapper, tmp_path)
        raw_dir = tmp_path / "raw"
        raw_dir.mkdir()

        build_calls = []
        monkeypatch.setattr(w, "_build_subclade",
                            lambda **k: build_calls.append(k))
        releaf_calls = []
        monkeypatch.setattr(w, "_run_releaf",
                            lambda **k: releaf_calls.append(k))

        query = _query(tmp_path, stem="GCF_query")
        query.update({
            "subclade_name": "Andreesenella_2",
            "parent_taxon": "Andreesenella",
            "subclade_id": 2,
            "database_dir": str(tmp_path / "databases" / "Andreesenella_2_db"),
            "subclade_taxonomy": "d__Bacteria;g__Andreesenella",
            "members_file": str(tmp_path / "members.txt"),
            "sketch_file": str(tmp_path / "sketch.msh"),
            "source_genome_dir": str(raw_dir),
            "taxonomy": "d__Bacteria;g__Andreesenella;s__some_query_species",
        })

        w._process_subclade_build("Andreesenella_2", [query])

        assert len(build_calls) == 1
        # Uses subclade_taxonomy, NOT the query's own (more specific) taxonomy.
        assert build_calls[0]["taxonomy"] == "d__Bacteria;g__Andreesenella"
        assert build_calls[0]["force"] is True
        assert build_calls[0]["query_assemblies"] == []
        assert len(releaf_calls) == 1

    def test_process_subclade_build_requires_source_genome_dir(
            self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path)
        query = _query(tmp_path)
        query.update({"subclade_name": "X", "database_dir": str(tmp_path / "X_db")})
        with pytest.raises(RuntimeError, match="source_genome_dir"):
            w._process_subclade_build("X", [query])


class TestRegisterLazySubclade:
    def test_writes_register_checkpoint_not_database_checkpoint(
            self, Wrapper, wrapper_module, tmp_path, monkeypatch):
        """Hazard 1: registration must use its OWN register_<name> checkpoint,
        distinct from _build_subclade's database_<name> key, so a later
        on-demand build is not skipped by a checkpoint set at registration
        time."""
        w = _make_wrapper(Wrapper, tmp_path)

        class FakeCompleted:
            returncode = 0
        monkeypatch.setattr(wrapper_module.subprocess, "run",
                            lambda *a, **k: FakeCompleted())

        entry = {"name": "Andreesenella_2", "subclade_id": 2, "n_genomes": 90,
                 "sketch_file": str(tmp_path / "s.msh"),
                 "members_file": str(tmp_path / "m.txt")}
        w._register_lazy_subclade(
            taxon_name="Andreesenella", entry=entry, raw_dir=tmp_path / "raw",
            taxonomy="d__Bacteria;g__Andreesenella")

        assert w._check_checkpoint("register_Andreesenella_2")
        assert not w._check_checkpoint("database_Andreesenella_2")
