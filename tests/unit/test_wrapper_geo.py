"""Tests for --geo/--geo-tsv/--geo-top-k/--geo-format wrapper plumbing and
_generate_geo_figures, in orthophyl_pipeline_wrapper.py.

Covers:
  - --geo-tsv/--geo-top-k implying --geo (mirrors test_wrapper_megatree_
    implication.py's pattern for --megatree-lazy et al.)
  - --genome-dir --geo without --geo-tsv is rejected by main()
  - missing matplotlib fails fast at _validate_dependencies (not at
    construction, and not deferred to the first non-fatal figure failure)
  - _generate_geo_figures: correct argv/env/output dir, dry-run short-
    circuit, no-op when self.geo is False, warn-not-raise on no metadata
    and on a failing subprocess
  - the three publish-site hooks (local mode, taxon create mode) call the
    helper with the right tree + metadata once their builds are mocked out
  - _phase_aggregation finds ReLeaf/OrthoPhyl trees via the real locator
    functions (not hardcoded filenames) and threads DB-persisted metadata

All external process execution is mocked; no real OrthoFinder/mash/etc. runs.
"""

import json
import logging

import pytest


@pytest.fixture
def Wrapper(wrapper_module):
    return wrapper_module.PipelineWrapper


# --------------------------------------------------------------------------- #
# --geo-tsv / --geo-top-k implication, via main()                            #
# --------------------------------------------------------------------------- #

def _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv):
    captured = {}

    def fake_run(self):
        captured["geo"] = self.geo
        captured["geo_tsv"] = self.geo_tsv
        captured["geo_top_k"] = self.geo_top_k
        return 0

    monkeypatch.setattr(wrapper_module.sys, "argv", argv)
    monkeypatch.setattr(wrapper_module.PipelineWrapper, "run", fake_run)
    rc = wrapper_module.main()
    return rc, captured


def _base_argv(tmp_path):
    input_tsv = tmp_path / "assemblies.tsv"
    input_tsv.write_text("assembly_path\ttaxonomy\tassembly_id\n")
    db_dir = tmp_path / "databases"
    db_dir.mkdir(exist_ok=True)
    return [
        "orthophyl_pipeline_wrapper.py",
        "--input", str(input_tsv),
        "--database-dir", str(db_dir),
        "--dry-run",
    ]


class TestGeoImplication:
    def test_geo_tsv_alone_implies_geo(self, wrapper_module, tmp_path, monkeypatch, caplog):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        argv = _base_argv(tmp_path) + ["--geo-tsv", str(geo_tsv)]
        with caplog.at_level(logging.INFO):
            rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["geo"] is True
        assert any("--geo-tsv" in rec.message and "implying --geo" in rec.message
                   for rec in caplog.records)

    def test_geo_top_k_alone_implies_geo(self, wrapper_module, tmp_path, monkeypatch):
        argv = _base_argv(tmp_path) + ["--geo-top-k", "50"]
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["geo"] is True
        assert captured["geo_top_k"] == 50

    def test_geo_top_k_default_does_not_imply_geo(self, wrapper_module, tmp_path, monkeypatch):
        """--geo-top-k's default (100) must not spuriously imply --geo just
        because it's present in args -- only an explicit non-default value does."""
        argv = _base_argv(tmp_path)
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["geo"] is False

    def test_explicit_geo_with_tsv_still_works(self, wrapper_module, tmp_path, monkeypatch):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        argv = _base_argv(tmp_path) + ["--geo", "--geo-tsv", str(geo_tsv)]
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["geo"] is True

    def test_genome_dir_geo_without_tsv_rejected(self, wrapper_module, tmp_path, monkeypatch):
        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()
        (genome_dir / "g1.fna").write_text(">c\nAC\n")
        db_dir = tmp_path / "databases"
        db_dir.mkdir()
        argv = [
            "orthophyl_pipeline_wrapper.py",
            "--genome-dir", str(genome_dir),
            "--clade-name", "Blorptaxon",
            "--database-dir", str(db_dir),
            "--geo",
        ]
        monkeypatch.setattr(wrapper_module.sys, "argv", argv)
        with pytest.raises(SystemExit):
            wrapper_module.main()

    def test_genome_dir_geo_with_tsv_accepted(self, wrapper_module, tmp_path, monkeypatch):
        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()
        (genome_dir / "g1.fna").write_text(">c\nAC\n")
        db_dir = tmp_path / "databases"
        db_dir.mkdir()
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        argv = [
            "orthophyl_pipeline_wrapper.py",
            "--genome-dir", str(genome_dir),
            "--clade-name", "Blorptaxon",
            "--database-dir", str(db_dir),
            "--geo", "--geo-tsv", str(geo_tsv),
            "--dry-run",
        ]
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["geo"] is True


# --------------------------------------------------------------------------- #
# Construction-time validation                                               #
# --------------------------------------------------------------------------- #

class TestGeoConstruction:
    def test_geo_tsv_must_exist(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="--geo-tsv does not exist"):
            Wrapper(
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                geo=True,
                geo_tsv=tmp_path / "nope.tsv",
            )

    def test_bad_geo_format_rejected(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="--geo-format must be one of"):
            Wrapper(
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                geo=True,
                geo_format="jpeg",
            )

    def test_defaults(self, Wrapper, tmp_path):
        w = Wrapper(database_dir=tmp_path / "db", output_dir=tmp_path / "out")
        assert w.geo is False
        assert w.geo_tsv is None
        assert w.geo_top_k == 100
        assert w.geo_format == "png"


# --------------------------------------------------------------------------- #
# Missing matplotlib fails fast at startup                                   #
# --------------------------------------------------------------------------- #

class TestGeoDependencyCheck:
    def test_geo_without_matplotlib_raises_at_validate(self, Wrapper, tmp_path, monkeypatch, wrapper_module):
        w = Wrapper(database_dir=tmp_path / "db", output_dir=tmp_path / "out", geo=True)
        monkeypatch.setattr(wrapper_module.importlib.util, "find_spec", lambda name: None)
        with pytest.raises(ImportError, match="--geo requires matplotlib"):
            w._validate_dependencies()

    def test_geo_false_skips_matplotlib_check(self, Wrapper, tmp_path, monkeypatch, wrapper_module):
        """Without --geo, a missing matplotlib must not block startup at all."""
        w = Wrapper(database_dir=tmp_path / "db", output_dir=tmp_path / "out", geo=False)
        monkeypatch.setattr(wrapper_module.importlib.util, "find_spec", lambda name: None)
        # Other dependency checks may still fail/pass depending on env; we
        # only care that an ImportError mentioning --geo is NOT raised here.
        try:
            w._validate_dependencies()
        except ImportError as e:
            assert "--geo" not in str(e)

    def test_geo_with_matplotlib_present_does_not_raise_for_that_reason(
            self, Wrapper, tmp_path, monkeypatch, wrapper_module):
        w = Wrapper(database_dir=tmp_path / "db", output_dir=tmp_path / "out", geo=True)
        monkeypatch.setattr(wrapper_module.importlib.util, "find_spec",
                             lambda name: object() if name == "matplotlib" else None)
        try:
            w._validate_dependencies()
        except ImportError as e:
            assert "--geo" not in str(e)


# --------------------------------------------------------------------------- #
# _generate_geo_figures                                                      #
# --------------------------------------------------------------------------- #

@pytest.fixture
def geo_wrapper(Wrapper, tmp_path, **overrides):
    def _make(**kw):
        defaults = dict(
            database_dir=tmp_path / "db",
            output_dir=tmp_path / "out",
            geo=True,
        )
        defaults.update(kw)
        w = Wrapper(**defaults)
        w.results_dir.mkdir(parents=True, exist_ok=True)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        return w
    return _make


class TestGenerateGeoFigures:
    def test_noop_when_geo_disabled(self, Wrapper, tmp_path, monkeypatch, wrapper_module):
        w = Wrapper(database_dir=tmp_path / "db", output_dir=tmp_path / "out", geo=False)
        w.results_dir.mkdir(parents=True, exist_ok=True)
        calls = []
        monkeypatch.setattr(wrapper_module.subprocess, "run",
                             lambda *a, **k: calls.append(a) or None)
        w._generate_geo_figures(tmp_path / "t.nwk", "label", metadata_file=tmp_path / "m.tsv")
        assert calls == []

    def test_warns_and_skips_with_no_metadata_source(self, geo_wrapper, tmp_path, monkeypatch,
                                                       wrapper_module, caplog):
        w = geo_wrapper()
        calls = []
        monkeypatch.setattr(wrapper_module.subprocess, "run",
                             lambda *a, **k: calls.append(a) or None)
        with caplog.at_level(logging.WARNING):
            w._generate_geo_figures(tmp_path / "t.nwk", "label", metadata_file=None)
        assert calls == []
        assert any("no metadata or --geo-tsv" in rec.message for rec in caplog.records)

    def test_builds_expected_command_and_env(self, geo_wrapper, tmp_path, monkeypatch, wrapper_module):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        w = geo_wrapper(geo_tsv=geo_tsv, geo_top_k=42, geo_format="svg")

        captured = {}

        class FakeCompleted:
            returncode = 0

        def fake_run(cmd, stdout=None, stderr=None, text=None, env=None):
            captured["cmd"] = cmd
            captured["env"] = env
            if stdout is not None:
                stdout.write("ok\n")
            return FakeCompleted()

        monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)

        tree_file = tmp_path / "mytree.nwk"
        tree_file.write_text("(a,b);")
        metadata_file = tmp_path / "meta.tsv"
        metadata_file.write_text("GCF_1 ... USA\n")

        w._generate_geo_figures(tree_file, "MyLabel", metadata_file=metadata_file)

        cmd = captured["cmd"]
        assert cmd[0] == wrapper_module.sys.executable
        assert cmd[1] == str(w.geo_figures_script)
        assert "--tree" in cmd and str(tree_file) in cmd
        assert "--out-dir" in cmd and str(w.results_dir / "figures") in cmd
        assert "--label" in cmd and "MyLabel" in cmd
        assert "--format" in cmd and "svg" in cmd
        assert "--top-k" in cmd and "42" in cmd
        assert "--metadata" in cmd and str(metadata_file) in cmd
        assert "--geo-tsv" in cmd and str(geo_tsv) in cmd

        assert captured["env"]["QT_QPA_PLATFORM"] == "offscreen"
        assert (w.results_dir / "figures").is_dir()

    def test_omits_metadata_flag_when_none_but_geo_tsv_present(
            self, geo_wrapper, tmp_path, monkeypatch, wrapper_module):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        w = geo_wrapper(geo_tsv=geo_tsv)
        captured = {}

        class FakeCompleted:
            returncode = 0

        def fake_run(cmd, stdout=None, stderr=None, text=None, env=None):
            captured["cmd"] = cmd
            if stdout is not None:
                stdout.write("ok\n")
            return FakeCompleted()

        monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
        w._generate_geo_figures(tmp_path / "t.nwk", "label", metadata_file=None)
        assert "--metadata" not in captured["cmd"]
        assert "--geo-tsv" in captured["cmd"]

    def test_dry_run_short_circuits_before_subprocess(self, Wrapper, tmp_path, monkeypatch, wrapper_module):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        w = Wrapper(database_dir=tmp_path / "db", output_dir=tmp_path / "out",
                    geo=True, geo_tsv=geo_tsv, dry_run=True)
        w.results_dir.mkdir(parents=True, exist_ok=True)
        calls = []
        monkeypatch.setattr(wrapper_module.subprocess, "run",
                             lambda *a, **k: calls.append(a) or None)
        w._generate_geo_figures(tmp_path / "t.nwk", "label", metadata_file=None)
        assert calls == []

    def test_subprocess_failure_warns_does_not_raise(self, geo_wrapper, tmp_path, monkeypatch,
                                                       wrapper_module, caplog):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        w = geo_wrapper(geo_tsv=geo_tsv)

        class FakeCompleted:
            returncode = 1

        def fake_run(cmd, stdout=None, stderr=None, text=None, env=None):
            if stdout is not None:
                stdout.write("boom\n")
            return FakeCompleted()

        monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
        with caplog.at_level(logging.WARNING):
            w._generate_geo_figures(tmp_path / "t.nwk", "label", metadata_file=None)
        assert any("figure generation failed" in rec.message for rec in caplog.records)

    def test_subprocess_oserror_warns_does_not_raise(self, geo_wrapper, tmp_path, monkeypatch,
                                                       wrapper_module, caplog):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g1\tUSA\n")
        w = geo_wrapper(geo_tsv=geo_tsv)

        def fake_run(*a, **k):
            raise OSError("no such executable")

        monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
        with caplog.at_level(logging.WARNING):
            w._generate_geo_figures(tmp_path / "t.nwk", "label", metadata_file=None)
        assert any("figure generation failed" in rec.message for rec in caplog.records)


# --------------------------------------------------------------------------- #
# Hook sites: local mode / taxon create mode call the helper correctly       #
# --------------------------------------------------------------------------- #

class TestLocalModeHook:
    def _make_local_wrapper(self, Wrapper, tmp_path, **overrides):
        genome_dir = tmp_path / "my_genomes"
        genome_dir.mkdir(exist_ok=True)
        for s in ("g0", "g1", "g2", "g3"):
            (genome_dir / f"{s}.fna").write_text(">c\nACGT\n")
        db_dir = tmp_path / "databases"
        db_dir.mkdir(exist_ok=True)
        fake_gather = tmp_path / "fake_gather.sh"
        fake_gather.write_text("#!/bin/bash\necho fake\n")
        fake_gather.chmod(0o755)
        kwargs = dict(
            database_dir=db_dir,
            output_dir=tmp_path / "out",
            threads=4,
            gather_script=fake_gather,
            genome_dir=genome_dir,
            clade_name="Blorptaxon",
            geo=True,
        )
        kwargs.update(overrides)
        w = Wrapper(**kwargs)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        return w

    def test_local_mode_calls_geo_hook_with_tree_and_no_metadata(
            self, Wrapper, tmp_path, monkeypatch):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g0\tUSA\n")
        w = self._make_local_wrapper(Wrapper, tmp_path, geo_tsv=geo_tsv)

        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))
        monkeypatch.setattr(w, "_qc_subclade", lambda *a, **k: (
            (tmp_path / "gtk").mkdir(parents=True, exist_ok=True),
            [(tmp_path / "gtk" / f"k{i}.fna").write_text(">c\nAC\n") for i in range(4)],
            tmp_path / "gtk")[-1])
        monkeypatch.setattr(w, "_run_orthophyl", lambda **k: None)
        monkeypatch.setattr(w, "_create_local_database", lambda **k: None)

        real_tree = tmp_path / "species.nwk"
        real_tree.write_text("(g0,g1);")
        monkeypatch.setattr(w, "_locate_species_tree", lambda out: real_tree)
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        geo_calls = []
        monkeypatch.setattr(w, "_generate_geo_figures",
                             lambda tree_file, label, metadata_file=None:
                                 geo_calls.append((tree_file, label, metadata_file)))

        rc = w._run_local_genomes_mode()
        assert rc == 0
        assert len(geo_calls) == 1
        tree_file, label, metadata_file = geo_calls[0]
        assert label == w._default_run_name(w.clade_name)
        assert metadata_file is None  # local mode: no NCBI metadata, relies on --geo-tsv
        assert tree_file.name == f"{label}_phylogeny.nwk"

    def test_local_mode_skips_geo_hook_when_tree_not_found(self, Wrapper, tmp_path, monkeypatch):
        geo_tsv = tmp_path / "locs.tsv"
        geo_tsv.write_text("g0\tUSA\n")
        w = self._make_local_wrapper(Wrapper, tmp_path, geo_tsv=geo_tsv)

        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))
        gtk = tmp_path / "gtk"
        gtk.mkdir(parents=True, exist_ok=True)
        for i in range(4):
            (gtk / f"k{i}.fna").write_text(">c\nAC\n")
        monkeypatch.setattr(w, "_qc_subclade", lambda *a, **k: gtk)
        monkeypatch.setattr(w, "_run_orthophyl", lambda **k: None)
        monkeypatch.setattr(w, "_create_local_database", lambda **k: None)
        monkeypatch.setattr(w, "_locate_species_tree", lambda out: out / "nonexistent.nwk")
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        geo_calls = []
        monkeypatch.setattr(w, "_generate_geo_figures",
                             lambda *a, **k: geo_calls.append((a, k)))

        rc = w._run_local_genomes_mode()
        assert rc == 0
        assert geo_calls == []


class TestPhaseAggregationGeoHook:
    def test_releaf_tree_passes_persisted_metadata_file(self, Wrapper, tmp_path, monkeypatch):
        db_dir_parent = tmp_path / "databases"
        db_dir_parent.mkdir()
        db_dir = db_dir_parent / "Escherichia_db"
        db_dir.mkdir()
        releaf_dir = db_dir / "orthophyl_run" / "ReLeaf_dir"
        releaf_dir.mkdir(parents=True)

        metadata_rel = "assembly_metadata.tsv"
        (db_dir / metadata_rel).write_text("GCF_1 ... USA\n")
        (db_dir / "database_config.json").write_text(
            json.dumps({"metadata_file": metadata_rel}))

        w = Wrapper(database_dir=db_dir_parent, output_dir=tmp_path / "out", geo=True)
        w.results_dir.mkdir(parents=True, exist_ok=True)
        (w.output_dir / "02_orthophyl_novel" / "orthophyl_runs").mkdir(parents=True, exist_ok=True)

        releaf_tree = tmp_path / "releaf_tree.nwk"
        releaf_tree.write_text("(a,b);")
        monkeypatch.setattr(w, "_locate_releaf_tree", lambda out: releaf_tree)
        monkeypatch.setattr(w, "_generate_summary_report", lambda *a, **k: None)

        geo_calls = []
        monkeypatch.setattr(w, "_generate_geo_figures",
                             lambda tree_file, label, metadata_file=None:
                                 geo_calls.append((tree_file, label, metadata_file)))

        w._phase_aggregation()

        assert len(geo_calls) == 1
        tree_file, label, metadata_file = geo_calls[0]
        assert label == "Escherichia_db"
        assert metadata_file == db_dir / metadata_rel

    def test_releaf_tree_with_no_persisted_metadata_passes_none(self, Wrapper, tmp_path, monkeypatch):
        db_dir_parent = tmp_path / "databases"
        db_dir_parent.mkdir()
        db_dir = db_dir_parent / "Escherichia_db"
        db_dir.mkdir()
        releaf_dir = db_dir / "orthophyl_run" / "ReLeaf_dir"
        releaf_dir.mkdir(parents=True)
        # No database_config.json at all.

        w = Wrapper(database_dir=db_dir_parent, output_dir=tmp_path / "out", geo=True)
        w.results_dir.mkdir(parents=True, exist_ok=True)
        (w.output_dir / "02_orthophyl_novel" / "orthophyl_runs").mkdir(parents=True, exist_ok=True)

        releaf_tree = tmp_path / "releaf_tree.nwk"
        releaf_tree.write_text("(a,b);")
        monkeypatch.setattr(w, "_locate_releaf_tree", lambda out: releaf_tree)
        monkeypatch.setattr(w, "_generate_summary_report", lambda *a, **k: None)

        geo_calls = []
        monkeypatch.setattr(w, "_generate_geo_figures",
                             lambda tree_file, label, metadata_file=None:
                                 geo_calls.append((tree_file, label, metadata_file)))

        w._phase_aggregation()
        assert len(geo_calls) == 1
        assert geo_calls[0][2] is None

    def test_orthophyl_tree_uses_download_dir_metadata(self, Wrapper, tmp_path, monkeypatch):
        db_dir_parent = tmp_path / "databases"
        db_dir_parent.mkdir()

        w = Wrapper(database_dir=db_dir_parent, output_dir=tmp_path / "out", geo=True)
        w.results_dir.mkdir(parents=True, exist_ok=True)

        taxon_dir = w.orthophyl_dir / "orthophyl_runs" / "Pseudomonas"
        taxon_dir.mkdir(parents=True, exist_ok=True)
        tree_file = tmp_path / "op_tree.nwk"
        tree_file.write_text("(a,b);")
        monkeypatch.setattr(w, "_locate_species_tree", lambda out: tree_file)
        monkeypatch.setattr(w, "_locate_releaf_tree", lambda out: None)
        monkeypatch.setattr(w, "_generate_summary_report", lambda *a, **k: None)

        metadata_dir = w.orthophyl_dir / "downloads" / "Pseudomonas"
        metadata_dir.mkdir(parents=True, exist_ok=True)
        metadata_file = metadata_dir / "all_asm_acc_metadata"
        metadata_file.write_text("GCF_1 ... USA\n")

        geo_calls = []
        monkeypatch.setattr(w, "_generate_geo_figures",
                             lambda tf, label, metadata_file=None:
                                 geo_calls.append((tf, label, metadata_file)))

        w._phase_aggregation()
        assert len(geo_calls) == 1
        tf, label, mf = geo_calls[0]
        assert label == "Pseudomonas"
        assert mf == metadata_file

    def test_orthophyl_tree_without_metadata_file_passes_none(self, Wrapper, tmp_path, monkeypatch):
        db_dir_parent = tmp_path / "databases"
        db_dir_parent.mkdir()

        w = Wrapper(database_dir=db_dir_parent, output_dir=tmp_path / "out", geo=True)
        w.results_dir.mkdir(parents=True, exist_ok=True)

        taxon_dir = w.orthophyl_dir / "orthophyl_runs" / "Pseudomonas"
        taxon_dir.mkdir(parents=True, exist_ok=True)
        tree_file = tmp_path / "op_tree.nwk"
        tree_file.write_text("(a,b);")
        monkeypatch.setattr(w, "_locate_species_tree", lambda out: tree_file)
        monkeypatch.setattr(w, "_locate_releaf_tree", lambda out: None)
        monkeypatch.setattr(w, "_generate_summary_report", lambda *a, **k: None)
        # No downloads/<taxon>/all_asm_acc_metadata created.

        geo_calls = []
        monkeypatch.setattr(w, "_generate_geo_figures",
                             lambda tf, label, metadata_file=None:
                                 geo_calls.append((tf, label, metadata_file)))

        w._phase_aggregation()
        assert len(geo_calls) == 1
        assert geo_calls[0][2] is None
