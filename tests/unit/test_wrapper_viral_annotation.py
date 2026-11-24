"""Tests for the wrapper's viral-annotation threading (--annotator viral,
--vadr-models, --vigor-db, --vigor-db-dir).

Covers:
  - __init__ validation (viral requires --genome-dir, rejects --megatree,
    requires --vadr-models/--vigor-db-dir)
  - --genome-dir + viral auto-forces skip_qc
  - _run_orthophyl passthrough of the viral flags
  - _run_releaf reads the TARGET DATABASE's own annotator (not self.annotator)
    and raises a clear error if this run has no vadr/vigor paths to use
  - OP_database_tool.py --annotator passthrough from _create_database_entry /
    _create_local_database

All external process execution is mocked -- no real OrthoPhyl.sh/ReLeaf.sh/
OP_database_tool.py is invoked.
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
        # ReLeaf.sh calls are checked for expected output files afterward --
        # fabricate them the same way test_wrapper_batch.py does.
        if "ReLeaf.sh" in " ".join(str(x) for x in cmd):
            for i, arg in enumerate(cmd):
                if arg == "-s" and i + 1 < len(cmd):
                    releaf_dir = Path(str(cmd[i + 1])) / "ReLeaf_dir"
                    releaf_dir.mkdir(parents=True, exist_ok=True)
                    (releaf_dir / "new_prot_alignments.trm.nm").write_text("fake")
                    (releaf_dir / "new_CDS_alignments.trm.nm").write_text("fake")
                    (releaf_dir / "new_trees").mkdir(exist_ok=True)
                    break
        return FakeCompleted()

    monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
    return calls


def _make_genome_dir(tmp_path, stems=("g0", "g1", "g2", "g3")):
    d = tmp_path / "my_genomes"
    d.mkdir(exist_ok=True)
    for s in stems:
        (d / f"{s}.fna").write_text(">c\nACGT\n")
    return d


def _make_local_wrapper(Wrapper, tmp_path, **overrides):
    db_dir = tmp_path / "databases"
    db_dir.mkdir(exist_ok=True)
    fake_gather = tmp_path / "fake_gather.sh"
    fake_gather.write_text("#!/bin/bash\necho fake\n")
    fake_gather.chmod(0o755)
    vigor_dir = tmp_path / "vigor_dbs"
    vigor_dir.mkdir(exist_ok=True)
    kwargs = dict(
        database_dir=db_dir,
        output_dir=tmp_path / "out",
        threads=4,
        gather_script=fake_gather,
        genome_dir=_make_genome_dir(tmp_path),
        clade_name="Noroviralia",
        annotator="viral",
        vadr_model_libs="/opt/vadr/calici-models",
        vigor_db_dir=vigor_dir,
    )
    kwargs.update(overrides)
    w = Wrapper(**kwargs)
    w.logs_dir.mkdir(parents=True, exist_ok=True)
    w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
    return w


class TestConstructionValidation:
    def test_defaults_to_prodigal(self, Wrapper, tmp_path):
        w = Wrapper(taxon="Methylorubrum", database_dir=tmp_path / "db",
                    output_dir=tmp_path / "out", threads=4)
        assert w.annotator == "prodigal"
        assert w.vadr_model_libs is None
        assert w.vigor_db == "auto"
        assert w.vigor_db_dir is None

    def test_invalid_annotator_rejected(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="annotator"):
            Wrapper(taxon="Methylorubrum", database_dir=tmp_path / "db",
                    output_dir=tmp_path / "out", threads=4, annotator="bogus")

    def test_viral_requires_vadr_models(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="--vadr-models"):
            _make_local_wrapper(Wrapper, tmp_path, vadr_model_libs=None)

    def test_viral_requires_vigor_db_dir(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="--vigor-db-dir"):
            _make_local_wrapper(Wrapper, tmp_path, vigor_db_dir=None)

    def test_viral_rejects_taxon_mode(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="genome-dir"):
            Wrapper(taxon="SomeVirus", database_dir=tmp_path / "db",
                    output_dir=tmp_path / "out", threads=4, annotator="viral",
                    vadr_model_libs="d", vigor_db_dir=tmp_path)

    def test_viral_rejects_input_batch_mode(self, Wrapper, tmp_path):
        fake_input = tmp_path / "input.tsv"
        fake_input.write_text("a\tfoo\n")
        with pytest.raises(ValueError, match="genome-dir"):
            Wrapper(input_file=fake_input, database_dir=tmp_path / "db",
                    output_dir=tmp_path / "out", threads=4, annotator="viral",
                    vadr_model_libs="d", vigor_db_dir=tmp_path)

    def test_viral_rejects_megatree(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="megatree"):
            _make_local_wrapper(Wrapper, tmp_path, megatree=True)

    def test_viral_auto_forces_skip_qc(self, Wrapper, tmp_path, caplog):
        w = _make_local_wrapper(Wrapper, tmp_path)
        assert w.skip_qc is True

    def test_explicit_skip_qc_stays_true(self, Wrapper, tmp_path):
        w = _make_local_wrapper(Wrapper, tmp_path, skip_qc=True)
        assert w.skip_qc is True

    def test_prodigal_genome_dir_does_not_force_skip_qc(self, Wrapper, tmp_path):
        w = _make_local_wrapper(Wrapper, tmp_path, annotator="prodigal",
                                 vadr_model_libs=None, vigor_db_dir=None)
        assert w.skip_qc is False


class TestRunOrthophylPassthrough:
    def test_viral_flags_appended(self, Wrapper, tmp_path, recording_run):
        w = _make_local_wrapper(Wrapper, tmp_path)
        input_dir = tmp_path / "op_input"
        input_dir.mkdir()
        (input_dir / "a.fna").write_text(">a\nACGT\n")
        output_dir = tmp_path / "op_out"
        w._run_orthophyl(input_dir, output_dir, "Noroviralia", assemblies=[])
        assert len(recording_run) == 1
        cmd = recording_run[0]
        assert "--annotator" in cmd and cmd[cmd.index("--annotator") + 1] == "viral"
        assert "--vadr-models" in cmd
        assert cmd[cmd.index("--vadr-models") + 1] == "/opt/vadr/calici-models"
        assert "--vigor-db-dir" in cmd

    def test_prodigal_omits_viral_flags(self, Wrapper, tmp_path, recording_run):
        w = Wrapper(taxon="Methylorubrum", database_dir=tmp_path / "db",
                    output_dir=tmp_path / "out", threads=4)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        input_dir = tmp_path / "op_input"
        input_dir.mkdir()
        (input_dir / "a.fna").write_text(">a\nACGT\n")
        output_dir = tmp_path / "op_out"
        w._run_orthophyl(input_dir, output_dir, "Methylorubrum", assemblies=[])
        cmd = recording_run[0]
        assert "--annotator" not in cmd


class TestRunReleafUsesTargetDatabaseAnnotator:
    def _common(self, Wrapper, tmp_path, **overrides):
        w = Wrapper(taxon="Methylorubrum", database_dir=tmp_path / "db",
                    output_dir=tmp_path / "out", threads=4, **overrides)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        db_dir = tmp_path / "db" / "SomeVirus_db"
        (db_dir / "orthophyl_run").mkdir(parents=True, exist_ok=True)
        input_genomes = tmp_path / "releaf_input"
        input_genomes.mkdir()
        return w, db_dir, input_genomes

    def test_viral_db_with_cli_paths_appends_flags(self, Wrapper, tmp_path, recording_run):
        w, db_dir, input_genomes = self._common(
            Wrapper, tmp_path,
            vadr_model_libs="/opt/vadr/calici-models",
            vigor_db_dir=tmp_path / "vigor_dbs")
        (tmp_path / "vigor_dbs").mkdir(exist_ok=True)
        w._run_releaf(
            database_dir=db_dir, input_genomes=input_genomes,
            output_dir=tmp_path / "releaf_out", tree_method="iqtree",
            tree_data="CDS", database_name="SomeVirus", annotator="viral")
        cmd = recording_run[0]
        assert "--annotator" in cmd and cmd[cmd.index("--annotator") + 1] == "viral"

    def test_viral_db_without_cli_paths_raises(self, Wrapper, tmp_path, recording_run):
        # self.annotator stays 'prodigal' (this run's own default) but the
        # matched DB is viral and no --vadr-models/--vigor-db-dir were given.
        w, db_dir, input_genomes = self._common(Wrapper, tmp_path)
        with pytest.raises(RuntimeError, match="vadr-models|vigor-db-dir"):
            w._run_releaf(
                database_dir=db_dir, input_genomes=input_genomes,
                output_dir=tmp_path / "releaf_out", tree_method="iqtree",
                tree_data="CDS", database_name="SomeVirus", annotator="viral")

    def test_prodigal_db_omits_flags_even_with_viral_cli_paths(
            self, Wrapper, tmp_path, recording_run):
        w, db_dir, input_genomes = self._common(
            Wrapper, tmp_path,
            vadr_model_libs="/opt/vadr/calici-models",
            vigor_db_dir=tmp_path / "vigor_dbs")
        (tmp_path / "vigor_dbs").mkdir(exist_ok=True)
        w._run_releaf(
            database_dir=db_dir, input_genomes=input_genomes,
            output_dir=tmp_path / "releaf_out", tree_method="iqtree",
            tree_data="CDS", database_name="SomeVirus", annotator="prodigal")
        cmd = recording_run[0]
        assert "--annotator" not in cmd

    def test_default_annotator_param_is_prodigal(self, Wrapper, tmp_path, recording_run):
        w, db_dir, input_genomes = self._common(Wrapper, tmp_path)
        w._run_releaf(
            database_dir=db_dir, input_genomes=input_genomes,
            output_dir=tmp_path / "releaf_out", tree_method="iqtree",
            tree_data="CDS", database_name="SomeVirus")
        cmd = recording_run[0]
        assert "--annotator" not in cmd


class TestCreateDatabaseEntryAnnotatorPassthrough:
    def test_local_database_passes_annotator(self, Wrapper, tmp_path, recording_run):
        w = _make_local_wrapper(Wrapper, tmp_path)
        orthophyl_output = tmp_path / "op_run"
        orthophyl_output.mkdir()
        w._create_database_entry(
            taxon_name="Noroviralia",
            orthophyl_output=orthophyl_output,
            taxonomy="g__Noroviralia",
            taxonomy_source="user_supplied",
            qc_applied=False,
            annotator="viral",
        )
        cmd = recording_run[0]
        assert "--annotator" in cmd and cmd[cmd.index("--annotator") + 1] == "viral"

    def test_default_annotator_param_is_prodigal(self, Wrapper, tmp_path, recording_run):
        w = Wrapper(taxon="Methylorubrum", database_dir=tmp_path / "db",
                    output_dir=tmp_path / "out", threads=4)
        w.logs_dir.mkdir(parents=True, exist_ok=True)
        w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        w.database_dir.mkdir(parents=True, exist_ok=True)
        orthophyl_output = tmp_path / "op_run"
        orthophyl_output.mkdir()
        w._create_database_entry(
            taxon_name="Methylorubrum",
            orthophyl_output=orthophyl_output,
            taxonomy="g__Methylorubrum",
        )
        cmd = recording_run[0]
        assert "--annotator" in cmd and cmd[cmd.index("--annotator") + 1] == "prodigal"
