"""Tests for large-taxon wiring in orthophyl_pipeline_wrapper.py.

The DEFAULT large-taxon behavior is diverse subsampling to one tree:
  * raw set under the ceiling -> single build (is_subclade=False), no subsample
  * raw set over the ceiling  -> diverse-subsample to --subsample-size, single build
  * create mode over ceiling  -> subsample, single taxon-flavored DB
  * gather split issues --download-only then --qc-only
The router emits ReLeaf / OrthoPhyl / OrthoPhyl_subclade_build decisions (the
third for lazily registered megatree subclades a query has just matched); the
opt-in --megatree build path itself is covered by
tests/unit/test_wrapper_megatree.py.
All heavy steps (mash, subsample, gather, OrthoPhyl, DB creator) are mocked.
"""

from pathlib import Path

import pytest


@pytest.fixture
def Wrapper(wrapper_module):
    return wrapper_module.PipelineWrapper


@pytest.fixture
def recording_run(monkeypatch, wrapper_module):
    calls = []

    class FakeCompleted:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, *args, **kwargs):
        calls.append([str(x) for x in cmd])
        return FakeCompleted()

    monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
    return calls


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


class TestProcessOrthophylTaxon:
    def test_under_ceiling_single_build(self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=10)

        # _download_raw materializes a small raw set (3 genomes).
        raw_dir = w.orthophyl_dir / "downloads" / "Andreesenella" / "assemblies_all.TMP"
        def fake_dl_raw(taxon, output_dir, query_assemblies=None):
            d = output_dir / "assemblies_all.TMP"
            d.mkdir(parents=True, exist_ok=True)
            for i in range(3):
                (d / f"g{i}.fna").write_text(">c\nAC\n")
            return d
        monkeypatch.setattr(w, "_download_raw", fake_dl_raw)

        part_calls = []
        monkeypatch.setattr(w, "_partition_genomes",
                            lambda *a, **k: part_calls.append(1))

        built = []
        monkeypatch.setattr(w, "_build_subclade",
                            lambda **k: built.append(k))

        w._process_orthophyl_taxon("Andreesenella", [_query(tmp_path)])
        # No partitioning under the ceiling.
        assert part_calls == []
        assert len(built) == 1
        assert built[0]["is_subclade"] is False

    def test_over_ceiling_subsamples_then_single_build(
            self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, subsample_size=4)

        def fake_dl_raw(taxon, output_dir, query_assemblies=None):
            d = output_dir / "assemblies_all.TMP"
            d.mkdir(parents=True, exist_ok=True)
            for i in range(8):  # 8 > ceiling (5) -> subsample
                (d / f"g{i}.fna").write_text(">c\nAC\n")
            return d
        monkeypatch.setattr(w, "_download_raw", fake_dl_raw)

        # Subsample returns a dir of 4 selected genomes; query stem seeds it.
        sub_calls = []
        def fake_subsample(taxon, raw_dir, target, must_keep_stems=None):
            sub_calls.append({"taxon": taxon, "target": target,
                              "seeds": must_keep_stems})
            sel = w.orthophyl_dir / "subsample" / taxon / "selected"
            sel.mkdir(parents=True, exist_ok=True)
            for i in range(4):
                (sel / f"g{i}.fna").write_text(">c\nAC\n")
            return sel
        monkeypatch.setattr(w, "_subsample_genomes", fake_subsample)
        # Partitioning must NOT be used on the default path.
        monkeypatch.setattr(w, "_partition_genomes",
                            lambda *a, **k: pytest.fail("default path must subsample, not partition"))

        built = []
        monkeypatch.setattr(w, "_build_subclade", lambda **k: built.append(k))

        w._process_orthophyl_taxon("Andreesenella", [_query(tmp_path)])

        # Subsampled once, to --subsample-size, seeded by the query stem.
        assert len(sub_calls) == 1
        assert sub_calls[0]["target"] == 4
        assert "GCF_query" in sub_calls[0]["seeds"]
        # One tree over the subsampled set (not a subclade).
        assert len(built) == 1
        assert built[0]["is_subclade"] is False
        assert built[0]["entry"]["name"] == "Andreesenella"

    def test_subsample_size_threads_through(self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, max_tree_genomes=5, subsample_size=42)

        def fake_dl_raw(taxon, output_dir, query_assemblies=None):
            d = output_dir / "assemblies_all.TMP"
            d.mkdir(parents=True, exist_ok=True)
            for i in range(8):
                (d / f"g{i}.fna").write_text(">c\nAC\n")
            return d
        monkeypatch.setattr(w, "_download_raw", fake_dl_raw)

        seen = {}
        def fake_subsample(taxon, raw_dir, target, must_keep_stems=None):
            seen["target"] = target
            sel = w.orthophyl_dir / "subsample" / taxon / "selected"
            sel.mkdir(parents=True, exist_ok=True)
            (sel / "g0.fna").write_text(">c\nAC\n")
            return sel
        monkeypatch.setattr(w, "_subsample_genomes", fake_subsample)
        monkeypatch.setattr(w, "_build_subclade", lambda **k: None)

        w._process_orthophyl_taxon("Andreesenella", [_query(tmp_path)])
        assert seen["target"] == 42


class TestGatherSplitArgv:
    def test_download_raw_issues_download_only(self, Wrapper, tmp_path, recording_run, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path)
        out = tmp_path / "dl"
        # download-only post-check needs raw genomes present.
        raw = out / "assemblies_all.TMP"
        raw.mkdir(parents=True)
        (raw / "g0.fna").write_text(">c\nAC\n")
        w._download_raw("Andreesenella", out)
        cmd = [c for c in recording_run if "fake_gather.sh" in " ".join(c)][0]
        assert "--download-only" in cmd
        assert "--qc-only" not in cmd

    def test_qc_subclade_issues_qc_only(self, Wrapper, tmp_path, recording_run, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path)
        subclade_dir = tmp_path / "sc"
        # qc post-check needs genomes_to_keep + stats.
        gtk = subclade_dir / "genomes_to_keep"
        gtk.mkdir(parents=True)
        (gtk / "g0.fna").write_text(">c\nAC\n")
        (subclade_dir / "assemblies_all.stats.txt").write_text("h\ng0\tx\n")

        raw = tmp_path / "raw"
        raw.mkdir()
        members = [raw / f"g{i}.fna" for i in range(3)]
        for m in members:
            m.write_text(">c\nAC\n")

        w._qc_subclade(subclade_dir, members, taxon_label="Andreesenella_1")
        cmd = [c for c in recording_run if "fake_gather.sh" in " ".join(c)][0]
        assert "--qc-only" in cmd
        assert "--download-only" not in cmd

    def test_must_keep_rides_qc_only_not_download(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, must_keep="GCF_x")
        # download-only call
        out = tmp_path / "dl"
        raw = out / "assemblies_all.TMP"
        raw.mkdir(parents=True)
        (raw / "g0.fna").write_text(">c\nAC\n")
        w._download_raw("Andreesenella", out)
        dl_cmd = [c for c in recording_run if "--download-only" in c][0]
        assert "--must-keep" not in dl_cmd

        recording_run.clear()
        # qc-only call
        subclade_dir = tmp_path / "sc"
        gtk = subclade_dir / "genomes_to_keep"
        gtk.mkdir(parents=True)
        (gtk / "g0.fna").write_text(">c\nAC\n")
        (subclade_dir / "assemblies_all.stats.txt").write_text("h\ng0\tx\n")
        rawm = tmp_path / "raw2"
        rawm.mkdir()
        members = [rawm / "g0.fna"]
        members[0].write_text(">c\nAC\n")
        w._qc_subclade(subclade_dir, members, taxon_label="Andreesenella_1")
        qc_cmd = [c for c in recording_run if "--qc-only" in c][0]
        assert "--must-keep" in qc_cmd


class TestCreateModeSubsamples:
    def test_over_ceiling_create_subsamples_single_db(self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, input_file=None, taxon="Andreesenella",
                          max_tree_genomes=5, subsample_size=4)

        # Stub the gatherer import path used inside _run_taxon_create_mode.
        class FakeGatherer:
            taxid = "123"
            taxon_rank = "genus"
            def __init__(self, **k): pass
            def get_taxonomy_string(self): return "d__Bacteria;g__Andreesenella"
        import types
        fake_mod = types.ModuleType("taxon_assembly_gatherer")
        fake_mod.TaxonAssemblyGatherer = FakeGatherer
        monkeypatch.setitem(__import__("sys").modules, "taxon_assembly_gatherer", fake_mod)

        def fake_dl_raw(taxon, output_dir, query_assemblies=None):
            d = output_dir / "assemblies_all.TMP"
            d.mkdir(parents=True, exist_ok=True)
            for i in range(8):  # 8 > ceiling (5) -> subsample
                (d / f"g{i}.fna").write_text(">c\nAC\n")
            return d
        monkeypatch.setattr(w, "_download_raw", fake_dl_raw)

        sub_calls = []
        def fake_subsample(taxon, raw_dir, target, must_keep_stems=None):
            sub_calls.append(target)
            sel = w.orthophyl_dir / "subsample" / taxon / "selected"
            sel.mkdir(parents=True, exist_ok=True)
            for i in range(4):
                (sel / f"g{i}.fna").write_text(">c\nAC\n")
            return sel
        monkeypatch.setattr(w, "_subsample_genomes", fake_subsample)
        monkeypatch.setattr(w, "_partition_genomes",
                            lambda *a, **k: pytest.fail("create mode must subsample, not partition"))

        # Single-tree QC + build + DB path.
        def fake_qc(subclade_dir, raw, taxon_label, query_assemblies=None):
            gtk = subclade_dir / "genomes_to_keep"
            gtk.mkdir(parents=True, exist_ok=True)
            for i in range(4):
                (gtk / f"g{i}.fna").write_text(">c\nAC\n")
            return gtk
        monkeypatch.setattr(w, "_qc_subclade", fake_qc)
        monkeypatch.setattr(w, "_run_orthophyl", lambda **k: None)
        db_calls = []
        monkeypatch.setattr(w, "_create_taxon_database", lambda **k: db_calls.append(k))
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        rc = w._run_taxon_create_mode()
        assert rc == 0
        assert sub_calls == [4]                     # subsampled to --subsample-size
        assert len(db_calls) == 1                    # one taxon-flavored DB
        assert db_calls[0]["taxon_name"] == "Andreesenella"


class TestTotalGenomeGuardrail:
    """The O(n^2) MASH matrix guardrail for the opt-in partition/megatree path.

    The DEFAULT large-taxon path subsamples and never builds the matrix, so it no
    longer routes through this check; the guardrail is exercised directly here and
    will gate the megatree path once that is wired (commit 3)."""

    def test_helper_message_reports_matrix_size(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, max_total_genomes=1000)
        with pytest.raises(RuntimeError) as exc:
            w._enforce_total_genome_ceiling("BigGenus", 50000)
        msg = str(exc.value)
        assert "50000" in msg
        assert "BigGenus" in msg
        # ~10 GB condensed array at 50k (50000*49999/2 * 8 bytes ~= 10 GB)
        assert "10.0 GB" in msg

    def test_at_ceiling_is_allowed(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, max_total_genomes=10)
        # Exactly at the ceiling must not raise.
        w._enforce_total_genome_ceiling("Foo", 10)


def _write_decision(routing_dir, idx, decision):
    import json
    routing_dir.mkdir(parents=True, exist_ok=True)
    (routing_dir / f"routing_decision_{idx}.json").write_text(json.dumps(decision))


class TestParseRoutingResults:
    """The router emits ReLeaf / OrthoPhyl / OrthoPhyl_subclade_build decisions;
    _parse_routing_results splits them three ways (the third for lazily
    registered megatree subclades a query has just matched -- see
    tests/unit/test_wrapper_megatree.py for the --megatree build path itself)."""

    def test_mixed_decisions_split_two_ways(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path)
        _write_decision(w.routing_dir, 0, {
            "pipeline": "ReLeaf", "assembly_id": "r1",
            "assembly": str(tmp_path / "r1.fna"),
            "matched_database": "Foo", "database_dir": str(tmp_path / "Foo_db"),
        })
        _write_decision(w.routing_dir, 1, {
            "pipeline": "OrthoPhyl", "assembly_id": "o1",
            "assembly": str(tmp_path / "o1.fna"),
            "download_value": "Bacillus", "download_rank": "g",
            "query_taxonomy": "d__Bacteria;g__Bacillus",
            "download_taxonomy": "d__Bacteria;g__Bacillus",
        })

        results = w._parse_routing_results()
        assert len(results["releaf_batch"]) == 1
        assert list(results["orthophyl_batch"].keys()) == ["Bacillus"]
        assert results["subclade_build_batch"] == {}

    def test_subclade_build_decision_grouped_by_subclade(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path)
        _write_decision(w.routing_dir, 0, {
            "pipeline": "OrthoPhyl_subclade_build", "assembly_id": "q1",
            "assembly": str(tmp_path / "q1.fna"),
            "subclade_name": "Bacillus_2", "parent_taxon": "Bacillus",
            "subclade_id": 2, "database_dir": str(tmp_path / "Bacillus_2_db"),
            "subclade_taxonomy": "d__Bacteria;g__Bacillus",
            "members_file": str(tmp_path / "Bacillus_2.members.txt"),
            "sketch_file": str(tmp_path / "Bacillus_2.msh"),
            "source_genome_dir": str(tmp_path / "raw"),
            "query_taxonomy": "d__Bacteria;g__Bacillus",
        })

        results = w._parse_routing_results()
        assert results["releaf_batch"] == []
        assert results["orthophyl_batch"] == {}
        assert list(results["subclade_build_batch"].keys()) == ["Bacillus_2"]
        entry = results["subclade_build_batch"]["Bacillus_2"][0]
        assert entry["assembly_id"] == "q1"
        assert entry["parent_taxon"] == "Bacillus"
        assert entry["subclade_taxonomy"] == "d__Bacteria;g__Bacillus"
