"""Tests for orthophyl_pipeline_wrapper.v2.py batch-mode orchestration (plan section 4).

All external process execution is mocked -- no real ReLeaf/OrthoPhyl/gather/NCBI is ever
invoked. Several tests here assert the *intended* behavior and therefore FAIL against the
current code, documenting known bugs (see WRAPPER_UNIT_TEST_PLAN.md, "Pre-existing
issues"):

  * test_releaf_actually_invoked_when_not_dry_run        -> bug B1
  * test_version_creation_runs_when_db_found             -> bugs B2/B3

These are marked xfail(strict=True) so the suite stays green today but will flip to a
hard failure (XPASS) the moment the bugs are fixed -- prompting removal of the marker.
"""

import json
from pathlib import Path

import pytest


@pytest.fixture
def Wrapper(wrapper_module):
    return wrapper_module.PipelineWrapper


@pytest.fixture
def recording_run(monkeypatch, wrapper_module):
    """Patch subprocess.run inside the wrapper module; record argv, return success."""
    calls = []

    class FakeCompleted:
        def __init__(self):
            self.returncode = 0
            self.stdout = ""
            self.stderr = ""

    def fake_run(cmd, *args, **kwargs):
        calls.append({"cmd": list(cmd), "args": args, "kwargs": kwargs})
        return FakeCompleted()

    monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
    return calls


def _make_wrapper(Wrapper, tmp_path, database_dir, **overrides):
    kwargs = dict(
        input_file=tmp_path / "assemblies.tsv",
        database_dir=database_dir,
        output_dir=tmp_path / "out",
        threads=4,
    )
    kwargs.update(overrides)
    return Wrapper(**kwargs)


class TestConstruction:
    def test_output_subdirs_computed(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        assert w.routing_dir == w.output_dir / "00_routing"
        assert w.releaf_dir == w.output_dir / "01_releaf_only"
        assert w.orthophyl_dir == w.output_dir / "02_orthophyl_novel"
        assert w.results_dir == w.output_dir / "03_results"
        assert w.checkpoint_dir == w.output_dir / "checkpoints"

    def test_script_paths_relative_to_wrapper(self, Wrapper, tmp_path, repo_root):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        assert w.script_dir == repo_root
        assert w.assembly_router.name == "assembly_router_multi.cmd_out3.py"
        assert w.orthophyl_script == repo_root / "OrthoPhyl.sh"

    def test_taxon_mode_flag(self, Wrapper, tmp_path):
        w = Wrapper(
            input_file=None,
            database_dir=tmp_path / "db",
            output_dir=tmp_path / "out",
            taxon="Methylorubrum",
        )
        assert w.taxon_mode is True


class TestValidateDependencies:
    def test_reports_missing_scripts(self, Wrapper, tmp_path, monkeypatch):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        # Point a required script at a nonexistent path.
        w.orthophyl_script = tmp_path / "does_not_exist_OrthoPhyl.sh"
        with pytest.raises(FileNotFoundError, match="OrthoPhyl"):
            w._validate_dependencies()

    def test_missing_gather_downgrades_to_none(self, Wrapper, tmp_path):
        w = _make_wrapper(
            Wrapper, tmp_path, tmp_path / "db",
            gather_script=tmp_path / "missing_gather.sh",
        )
        # gather script missing is a warning, not fatal; it gets reset to None.
        w._validate_dependencies()
        assert w.gather_script is None


class TestParseRoutingResults:
    def _write_decision(self, routing_dir, decision):
        routing_dir.mkdir(parents=True, exist_ok=True)
        f = routing_dir / f"routing_decision_{decision['assembly_id']}.json"
        f.write_text(json.dumps(decision))

    def test_parse_releaf_and_orthophyl(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        self._write_decision(w.routing_dir, {
            "pipeline": "ReLeaf",
            "assembly_id": "m1",
            "assembly": "/data/m1.fna",
            "matched_database": "Escherichia",
            "database_dir": "/db/Escherichia_db",
            "tree_method": "iqtree",
            "tree_data": "CDS",
        })
        self._write_decision(w.routing_dir, {
            "pipeline": "OrthoPhyl",
            "assembly_id": "n1",
            "assembly": "/data/n1.fna",
            "download_value": "Gaiella",
            "download_rank": "genus",
            "query_taxonomy": "d__Bacteria;...",
            "download_taxonomy": "d__Bacteria;...;g__Gaiella",
        })
        results = w._parse_routing_results()
        assert len(results["releaf_batch"]) == 1
        assert results["releaf_batch"][0]["database"] == "Escherichia"
        assert list(results["orthophyl_batch"].keys()) == ["Gaiella"]

    def test_releaf_defaults_applied(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        self._write_decision(w.routing_dir, {
            "pipeline": "ReLeaf",
            "assembly_id": "m1",
            "assembly": "/data/m1.fna",
            "matched_database": "Escherichia",
            "database_dir": "/db/Escherichia_db",
            # tree_method / tree_data intentionally omitted
        })
        results = w._parse_routing_results()
        assert results["releaf_batch"][0]["tree_method"] == "iqtree"
        assert results["releaf_batch"][0]["tree_data"] == "CDS"


class TestReleafPhase:
    def _prepare(self, Wrapper, tmp_path, make_db_dir):
        from conftest import ESCHERICHIA_TAX
        parent = tmp_path / "db"
        make_db_dir(
            clade_name="Escherichia", clade_taxonomy=ESCHERICHIA_TAX,
            clade_rank="g", clade_rank_name="genus", parent=parent,
        )
        asm = tmp_path / "m1.fna"
        asm.write_text(">c\nACGT\n")
        w = _make_wrapper(Wrapper, tmp_path, parent)
        # Directory structure the phase relies on.
        for d in (w.releaf_dir, w.logs_dir, w.checkpoint_dir):
            d.mkdir(parents=True, exist_ok=True)
        batch = [{
            "assembly_id": "m1",
            "assembly_path": str(asm),
            "database": "Escherichia",
            "database_dir": str(parent / "Escherichia_db"),
            "tree_method": "iqtree",
            "tree_data": "CDS",
        }]
        return w, batch

    def test_copies_assemblies_into_input_genomes(
        self, Wrapper, tmp_path, make_db_dir, recording_run
    ):
        w, batch = self._prepare(Wrapper, tmp_path, make_db_dir)
        w._phase_releaf(batch)
        assert (w.releaf_dir / "Escherichia" / "input_genomes" / "m1.fna").exists()

    def test_dry_run_does_not_invoke_releaf(
        self, Wrapper, tmp_path, make_db_dir, recording_run
    ):
        w, batch = self._prepare(Wrapper, tmp_path, make_db_dir)
        w.dry_run = True
        w._phase_releaf(batch)
        releaf_calls = [c for c in recording_run if "ReLeaf.sh" in " ".join(c["cmd"])]
        assert releaf_calls == []

    def test_releaf_actually_invoked_when_not_dry_run(
        self, Wrapper, tmp_path, make_db_dir, recording_run
    ):
        # Regression test for bug B1 (fixed): the misindented 'return' in _run_releaf
        # used to short-circuit execution so ReLeaf.sh was never invoked.
        w, batch = self._prepare(Wrapper, tmp_path, make_db_dir)
        w.dry_run = False
        w._phase_releaf(batch)
        releaf_calls = [c for c in recording_run if "ReLeaf.sh" in " ".join(c["cmd"])]
        assert releaf_calls, "expected ReLeaf.sh to be invoked via subprocess.run"

    def test_releaf_command_shape(
        self, Wrapper, tmp_path, make_db_dir, recording_run
    ):
        w, batch = self._prepare(Wrapper, tmp_path, make_db_dir)
        w.dry_run = False
        w._phase_releaf(batch)
        releaf_calls = [c for c in recording_run if "ReLeaf.sh" in " ".join(c["cmd"])]
        assert releaf_calls, "expected a ReLeaf.sh invocation"
        cmd = releaf_calls[0]["cmd"]
        assert "--store" in cmd
        assert "--input_genomes" in cmd
        assert "--tree_method" in cmd and "iqtree" in cmd
        assert "--TREE_DATA" in cmd and "CDS" in cmd

    def test_version_creation_runs_when_db_found(
        self, Wrapper, tmp_path, make_db_dir, recording_run
    ):
        # Regression test for bugs B2/B3 (fixed): unconditional 'return's in
        # _create_releaf_version made version creation dead code. With a matching
        # *_db present and the versioner script existing, the versioner must run.
        w, batch = self._prepare(Wrapper, tmp_path, make_db_dir)
        w.dry_run = False
        w._phase_releaf(batch)
        version_calls = [
            c for c in recording_run
            if "add_releaf_version.py" in " ".join(c["cmd"])
        ]
        assert version_calls, "expected add_releaf_version.py to be invoked"
        cmd = version_calls[0]["cmd"]
        assert "--database-dir" in cmd
        assert "--releaf-output" in cmd


class TestCheckpoints:
    def test_write_then_check_roundtrip(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        assert w._check_checkpoint("phase_x") is False
        w._write_checkpoint("phase_x")
        assert w._check_checkpoint("phase_x") is True

    def test_save_final_status_writes_json(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        w.output_dir.mkdir(parents=True, exist_ok=True)
        w._save_final_status()
        status = json.loads((w.output_dir / "pipeline_status.json").read_text())
        assert "start_time" in status and "end_time" in status
