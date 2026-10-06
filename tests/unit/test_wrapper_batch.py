"""Tests for orthophyl_pipeline_wrapper.py batch-mode orchestration (plan section 4).

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
        assert w.assembly_router.name == "assembly_router.py"
        assert w.orthophyl_script == repo_root / "OrthoPhyl.sh"

    def test_taxon_mode_flag(self, Wrapper, tmp_path):
        w = Wrapper(
            input_file=None,
            database_dir=tmp_path / "db",
            output_dir=tmp_path / "out",
            taxon="Methylorubrum",
        )
        assert w.taxon_mode is True


class TestGatherScriptDefault:
    """gather_script defaults to the bundled utils/gather_filter_asms.sh when
    --gather-script is omitted, instead of requiring explicit opt-in every
    invocation."""

    def test_defaults_to_bundled_script(self, Wrapper, tmp_path, repo_root):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        assert w.gather_script == repo_root / "utils" / "gather_filter_asms.sh"

    def test_explicit_gather_script_overrides_default(self, Wrapper, tmp_path):
        custom = tmp_path / "my_custom_gather.sh"
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db", gather_script=custom)
        assert w.gather_script == custom

    def test_bundled_default_passes_validate_dependencies(self, Wrapper, tmp_path):
        # The real bundled script must itself be exist+executable, or every
        # taxon-create/local-mode run would hard-fail out of the box.
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        w._validate_dependencies()  # must not raise
        assert w.gather_script is not None


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
        # Batch mode never requires a gather script -- missing is a warning,
        # not fatal; it gets reset to None.
        w._validate_dependencies()
        assert w.gather_script is None


class TestGatherScriptViabilityFailsFast:
    """Taxon create mode and local-genome-ingest mode (without --skip-qc) both
    hard-require a working gather script. _validate_dependencies must catch a
    missing/non-executable script and raise BEFORE any taxonomy-resolution
    work (NCBI taxdump download) happens -- it is called from
    _phase_initialization, which runs before mode dispatch in run()."""

    def test_taxon_create_mode_missing_script_raises(self, Wrapper, tmp_path):
        w = Wrapper(
            input_file=None, database_dir=tmp_path / "db",
            output_dir=tmp_path / "out", taxon="Methylorubrum", threads=4,
            gather_script=tmp_path / "does_not_exist.sh",
        )
        with pytest.raises(FileNotFoundError, match="Gather script"):
            w._validate_dependencies()

    def test_taxon_create_mode_non_executable_script_raises(self, Wrapper, tmp_path):
        bad = tmp_path / "not_executable.sh"
        bad.write_text("#!/bin/bash\necho hi\n")
        bad.chmod(0o644)  # no +x
        w = Wrapper(
            input_file=None, database_dir=tmp_path / "db",
            output_dir=tmp_path / "out", taxon="Methylorubrum", threads=4,
            gather_script=bad,
        )
        with pytest.raises(FileNotFoundError, match="not executable"):
            w._validate_dependencies()

    def test_taxon_create_mode_no_script_at_all_raises(self, Wrapper, tmp_path):
        # gather_script=None at construction now falls back to the bundled
        # utils/gather_filter_asms.sh (see TestGatherScriptDefault below), so
        # there's no longer a public-API way to end up with gather_script
        # actually None. Exercise that defensive branch directly instead.
        w = Wrapper(
            input_file=None, database_dir=tmp_path / "db",
            output_dir=tmp_path / "out", taxon="Methylorubrum", threads=4,
        )
        w.gather_script = None
        with pytest.raises(FileNotFoundError, match="requires a gather script"):
            w._validate_dependencies()

    def test_taxon_update_mode_missing_script_is_soft(self, Wrapper, tmp_path):
        """Update mode (an existing DB matches) tolerates a missing gather
        script -- it only skips QC on the new assemblies, same as before."""
        db_dir = tmp_path / "db"
        db_dir.mkdir()
        existing_db = db_dir / "Methylorubrum_genus_db"
        existing_db.mkdir()
        (existing_db / "database_config.json").write_text(json.dumps({
            "clade_name": "Methylorubrum", "source_taxon_name": "Methylorubrum",
            "assembly_accessions": ["GCF_000001.1"],
        }))
        w = Wrapper(
            input_file=None, database_dir=db_dir,
            output_dir=tmp_path / "out", taxon="Methylorubrum",
            update_existing=True, threads=4,
            gather_script=tmp_path / "does_not_exist.sh",
        )
        w._validate_dependencies()  # must not raise
        assert w.gather_script is None

    def test_local_mode_default_qc_missing_script_raises(self, Wrapper, tmp_path):
        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()
        w = Wrapper(
            input_file=None, database_dir=tmp_path / "db",
            output_dir=tmp_path / "out",
            genome_dir=genome_dir, clade_name="Blorptaxon", threads=4,
            gather_script=tmp_path / "does_not_exist.sh",
        )
        with pytest.raises(FileNotFoundError, match="Gather script"):
            w._validate_dependencies()

    def test_local_mode_skip_qc_missing_script_is_soft(self, Wrapper, tmp_path):
        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()
        w = Wrapper(
            input_file=None, database_dir=tmp_path / "db",
            output_dir=tmp_path / "out",
            genome_dir=genome_dir, clade_name="Blorptaxon", threads=4,
            skip_qc=True,
            gather_script=tmp_path / "does_not_exist.sh",
        )
        w._validate_dependencies()  # must not raise
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

    def _write_decision_raw(self, routing_dir, filename, decision):
        """Like _write_decision, but the caller picks the ON-DISK filename
        directly instead of deriving it from decision['assembly_id'] -- needed
        to simulate a hand-edited/older-router-produced JSON whose assembly_id
        field itself carries an unsanitized value."""
        routing_dir.mkdir(parents=True, exist_ok=True)
        (routing_dir / filename).write_text(json.dumps(decision))

    def test_defense_in_depth_sanitizes_download_value(self, Wrapper, tmp_path):
        """Even if a routing_decision_*.json on disk was NOT sanitized (e.g.
        hand-edited, or produced by an older assembly_router.py), taxon_name
        derived from it must still be safe -- this taxon_name feeds real Path
        joins like self.orthophyl_dir / "downloads" / taxon_name."""
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        self._write_decision_raw(w.routing_dir, "routing_decision_n1.json", {
            "pipeline": "OrthoPhyl",
            "assembly_id": "n1",
            "assembly": "/data/n1.fna",
            "download_value": "../../../../tmp/evil_taxon",
            "download_rank": "genus",
            "query_taxonomy": "d__Bacteria;...",
            "download_taxonomy": "d__Bacteria;...;g__evil_taxon",
        })
        results = w._parse_routing_results()
        taxon_names = list(results["orthophyl_batch"].keys())
        assert len(taxon_names) == 1
        assert "/" not in taxon_names[0]
        assert taxon_names[0] not in (".", "..")

    def test_defense_in_depth_sanitizes_assembly_id(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        self._write_decision_raw(w.routing_dir, "routing_decision_evil.json", {
            "pipeline": "ReLeaf",
            "assembly_id": "../../../../tmp/evil_id",
            "assembly": "/data/m1.fna",
            "matched_database": "Escherichia",
            "database_dir": "/db/Escherichia_db",
        })
        results = w._parse_routing_results()
        assembly_id = results["releaf_batch"][0]["assembly_id"]
        assert "/" not in assembly_id
        assert assembly_id not in (".", "..")

    def test_defense_in_depth_sanitizes_subclade_build_names(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        self._write_decision_raw(w.routing_dir, "routing_decision_evil2.json", {
            "pipeline": "OrthoPhyl_subclade_build",
            "assembly_id": "../../../../tmp/evil_id2",
            "assembly": "/data/q1.fna",
            "subclade_name": "../../../../tmp/evil_subclade",
            "parent_taxon": "../../../../tmp/evil_parent",
            "database_dir": "/db/Foo_1_db",
        })
        results = w._parse_routing_results()
        sc_names = list(results["subclade_build_batch"].keys())
        assert len(sc_names) == 1
        assert "/" not in sc_names[0]
        entry = results["subclade_build_batch"][sc_names[0]][0]
        assert "/" not in entry["assembly_id"]
        assert "/" not in entry["parent_taxon"]


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
        # Fixed flags: -s (not --store), -g (not --input_genomes), -p (not --tree_method), -o (not --TREE_DATA)
        assert "-s" in cmd
        assert "-g" in cmd
        assert "-p" in cmd and "iqtree" in cmd
        assert "-o" in cmd and "CDS" in cmd

    def test_version_creation_runs_when_db_found(
        self, Wrapper, tmp_path, make_db_dir, recording_run, monkeypatch
    ):
        # Regression test for bugs B2/B3 (fixed): unconditional 'return's in
        # _create_releaf_version made version creation dead code. With a matching
        # *_db present and the versioner script existing, the versioner must run.
        w, batch = self._prepare(Wrapper, tmp_path, make_db_dir)
        w.dry_run = False
        
        # Mock subprocess.run to also create expected ReLeaf output files
        original_run = recording_run
        def fake_run_with_outputs(cmd, *args, **kwargs):
            result = {"cmd": list(cmd), "args": args, "kwargs": kwargs}
            original_run.append(result)
            
            # If this is a ReLeaf.sh call, create the expected output files
            if "ReLeaf.sh" in " ".join(cmd):
                # Find the -s/--storage_dir argument
                for i, arg in enumerate(cmd):
                    if arg == "-s" and i + 1 < len(cmd):
                        storage_dir = Path(cmd[i + 1])
                        releaf_dir = storage_dir / "ReLeaf_dir"
                        releaf_dir.mkdir(parents=True, exist_ok=True)
                        (releaf_dir / "new_prot_alignments.trm.nm").write_text("fake")
                        (releaf_dir / "new_CDS_alignments.trm.nm").write_text("fake")
                        (releaf_dir / "new_trees").mkdir(exist_ok=True)
                        break
            
            class FakeCompleted:
                returncode = 0
                stdout = ""
                stderr = ""
            return FakeCompleted()
        
        import subprocess
        monkeypatch.setattr(subprocess, "run", fake_run_with_outputs)
        
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


class TestRunOrthophylCommand:
    """_run_orthophyl must always pass -n low enough to force OrthoPhyl.sh's
    MASH-shortlist/HMM-building branch (ANI_ORTHOFINDER_TO_ALL_SEQS) -- otherwise
    a database built from a small genome set never gets HMMs and can never be
    ReLeaf'd onto later. See ani_shortlist's docstring in __init__."""

    def _make_genome_dir(self, tmp_path, n, ext=".fna"):
        d = tmp_path / "genomes_in"
        d.mkdir(parents=True, exist_ok=True)
        for i in range(n):
            (d / f"g{i}{ext}").write_text(">c\nACGT\n")
        return d

    def test_small_input_capped_by_count_minus_one(
            self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db", ani_shortlist=20)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        input_dir = self._make_genome_dir(tmp_path, 6)

        w._run_orthophyl(
            input_dir=input_dir, output_dir=tmp_path / "out_run",
            taxon_name="Foo", assemblies=[])

        cmd = recording_run[0]["cmd"]
        assert "-n" in cmd
        assert cmd[cmd.index("-n") + 1] == "5"  # min(20, 6-1)

    def test_large_input_capped_by_ani_shortlist(
            self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db", ani_shortlist=20)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        input_dir = self._make_genome_dir(tmp_path, 50)

        w._run_orthophyl(
            input_dir=input_dir, output_dir=tmp_path / "out_run",
            taxon_name="Foo", assemblies=[])

        cmd = recording_run[0]["cmd"]
        assert "-n" in cmd
        assert cmd[cmd.index("-n") + 1] == "20"  # min(20, 50-1)

    def test_custom_ani_shortlist_threaded_through(
            self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db", ani_shortlist=5)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        input_dir = self._make_genome_dir(tmp_path, 50)

        w._run_orthophyl(
            input_dir=input_dir, output_dir=tmp_path / "out_run",
            taxon_name="Foo", assemblies=[])

        cmd = recording_run[0]["cmd"]
        assert cmd[cmd.index("-n") + 1] == "5"  # min(5, 50-1)


class TestCreateDatabaseEntryLogFile:
    """_create_database_entry builds log_file = self.logs_dir /
    f"database_{taxon_name}.log" directly from taxon_name (which can carry a
    batch TSV's free-text taxonomy value verbatim) in BOTH the subclade/backbone
    branch and the classic --update branch. A '/' in taxon_name used to crash
    with FileNotFoundError (open() doesn't create missing parent dirs) rather
    than just silently mis-routing like the db_dir lookups fixed elsewhere."""

    def test_update_branch_log_file_stays_inside_logs_dir(
            self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        w.database_dir.mkdir(parents=True, exist_ok=True)

        w._create_database_entry(
            taxon_name="Foo/Bar",
            orthophyl_output=tmp_path / "op_out",
            taxonomy="g__FooBar",
        )

        logs = list(w.logs_dir.glob("database_*.log"))
        assert len(logs) == 1
        assert "/" not in logs[0].name.removeprefix("database_").removesuffix(".log")
        assert logs[0].parent == w.logs_dir

    def test_subclade_branch_log_file_stays_inside_logs_dir(
            self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, tmp_path / "db")
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        w.database_dir.mkdir(parents=True, exist_ok=True)

        w._create_database_entry(
            taxon_name="Foo/Bar_1",
            orthophyl_output=tmp_path / "op_out",
            taxonomy="g__FooBar",
            subclade_meta={"is_backbone": False, "parent_taxon": "Foo/Bar"},
        )

        logs = list(w.logs_dir.glob("database_*.log"))
        assert len(logs) == 1
        assert "/" not in logs[0].name.removeprefix("database_").removesuffix(".log")
        assert logs[0].parent == w.logs_dir
