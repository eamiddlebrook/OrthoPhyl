"""Tests for the wrapper's tree-building speed controls (--tree-method,
--iqtree-model, --iqtree-mset, --iqtree-fast, --iqtree-options).

All external process execution is mocked -- no real OrthoPhyl.sh is invoked.
"""

import subprocess

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


def _make_wrapper(Wrapper, tmp_path, **overrides):
    kwargs = dict(
        taxon="Methylorubrum",
        database_dir=tmp_path / "db",
        output_dir=tmp_path / "out",
        threads=4,
    )
    kwargs.update(overrides)
    w = Wrapper(**kwargs)
    for d in (w.logs_dir, w.checkpoint_dir):
        d.mkdir(parents=True, exist_ok=True)
    return w


def _run(w, tmp_path, taxon_name="Methylorubrum"):
    input_dir = tmp_path / "input"
    input_dir.mkdir(exist_ok=True)
    (input_dir / "a.fna").write_text(">a\nACGT\n")
    output_dir = tmp_path / "op_out"
    w._run_orthophyl(input_dir, output_dir, taxon_name, assemblies=[])


class TestConstruction:
    def test_defaults(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path)
        assert w.tree_method == "iqtree"
        assert w.iqtree_model is None
        assert w.iqtree_mset is None
        assert w.iqtree_fast is False
        assert w.iqtree_options is None

    def test_invalid_tree_method_rejected(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="tree-method"):
            _make_wrapper(Wrapper, tmp_path, tree_method="raxml")

    def test_iqtree_flags_with_fasttree_rejected(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="fasttree"):
            _make_wrapper(Wrapper, tmp_path, tree_method="fasttree", iqtree_model="GTR+F+G4")

    def test_model_and_mset_mutually_exclusive(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="mutually exclusive"):
            _make_wrapper(Wrapper, tmp_path, iqtree_model="GTR+F+G4", iqtree_mset="GTR")


class TestRunOrthophylCommandShape:
    def test_default_command_is_unchanged(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path)
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        assert "-p" in cmd
        assert cmd[cmd.index("-p") + 1] == "iqtree"
        assert "-c" not in cmd

    def test_tree_method_fasttree(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, tree_method="fasttree")
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        assert cmd[cmd.index("-p") + 1] == "fasttree"
        assert "-c" not in cmd

    def test_tree_method_both(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, tree_method="both")
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        assert cmd[cmd.index("-p") + 1] == "fasttree iqtree"

    def test_iqtree_model_adds_control_file(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, iqtree_model="GTR+F+G4")
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        assert "-c" in cmd
        control_path = cmd[cmd.index("-c") + 1]
        content = open(control_path).read()
        assert "IQtree_CDS_speciestree_options" in content
        assert "-m GTR+F+G4" in content
        assert "$IQtree_CDS_speciestree_options" in content

    def test_iqtree_mset_adds_control_file(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, iqtree_mset="GTR,HKY")
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        control_path = cmd[cmd.index("-c") + 1]
        content = open(control_path).read()
        assert "--mset GTR,HKY" in content

    def test_iqtree_fast_with_clean_options_adds_control_file(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, iqtree_fast=True,
                           iqtree_options="--symtest --lmap 2000 -t PARS --ninit 2")
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        control_path = cmd[cmd.index("-c") + 1]
        content = open(control_path).read()
        assert "--fast" in content
        assert "$IQtree_CDS_speciestree_options" not in content

    def test_iqtree_options_replaces_base(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, iqtree_options="-m GTR+F+G4")
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        control_path = cmd[cmd.index("-c") + 1]
        content = open(control_path).read()
        assert "$IQtree_CDS_speciestree_options" not in content
        assert "-m GTR+F+G4" in content

    def test_combined_flags_appended_in_order(self, Wrapper, tmp_path, recording_run):
        w = _make_wrapper(Wrapper, tmp_path, iqtree_options="-t PARS", iqtree_fast=True)
        _run(w, tmp_path)
        cmd = recording_run[0]["cmd"]
        control_path = cmd[cmd.index("-c") + 1]
        content = open(control_path).read()
        assert "-t PARS --fast" in content

    def test_iqtree_fast_alone_rejected_due_to_default_ufboot(self, Wrapper, tmp_path):
        # The default IQtree_CDS_speciestree_options always includes -B 1000
        # (UltraFast Bootstrap), which IQ-TREE refuses to combine with --fast.
        with pytest.raises(ValueError, match="UltraFast Bootstrap"):
            _make_wrapper(Wrapper, tmp_path, iqtree_fast=True)

    def test_iqtree_fast_with_bootstrap_options_rejected(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="UltraFast Bootstrap"):
            _make_wrapper(Wrapper, tmp_path, iqtree_fast=True, iqtree_options="-B 1000")


class TestControlFileSourcing:
    """The control file is meant to be sourced by bash after
    control_file.defaults, so $IQtree_CDS_speciestree_options must expand to
    the pre-existing default rather than being treated literally, and
    adversarial content must not escape the double-quoted assignment.

    These call _tree_control_file directly (no subprocess mocking) since
    recording_run patches the process-wide subprocess module -- going
    through a mocked _run_orthophyl would also intercept this test's own
    real bash invocations used to verify sourcing behavior.
    """

    def test_sourced_value_appends_to_existing_default(self, Wrapper, tmp_path):
        w = _make_wrapper(Wrapper, tmp_path, iqtree_model="GTR+F+G4")
        control_path = w._tree_control_file("Methylorubrum")
        out = subprocess.run(
            ["bash", "-c",
             f'IQtree_CDS_speciestree_options="BASE"; source {control_path}; '
             f'echo "$IQtree_CDS_speciestree_options"'],
            capture_output=True, text=True, check=True,
        )
        assert out.stdout.strip() == "BASE -m GTR+F+G4"

    def test_malicious_iqtree_options_stays_inert(self, Wrapper, tmp_path):
        marker = tmp_path / "pwned"
        w = _make_wrapper(Wrapper, tmp_path, iqtree_options=f'"; touch {marker}; echo "')
        control_path = w._tree_control_file("Methylorubrum")
        subprocess.run(
            ["bash", "-c",
             f'IQtree_CDS_speciestree_options="BASE"; source {control_path}'],
            capture_output=True, text=True, check=True,
        )
        assert not marker.exists()

    def test_malicious_iqtree_model_stays_inert(self, Wrapper, tmp_path):
        marker = tmp_path / "pwned2"
        w = _make_wrapper(Wrapper, tmp_path, iqtree_model=f'$(touch {marker})')
        control_path = w._tree_control_file("Methylorubrum")
        subprocess.run(
            ["bash", "-c",
             f'IQtree_CDS_speciestree_options="BASE"; source {control_path}'],
            capture_output=True, text=True, check=True,
        )
        assert not marker.exists()


class TestMainArgParseValidation:
    def _parse(self, wrapper_module, monkeypatch, tmp_path, argv):
        monkeypatch.setattr(wrapper_module.sys, "argv", ["orthophyl_pipeline_wrapper.py"] + argv)
        with pytest.raises(SystemExit) as exc_info:
            wrapper_module.main()
        return exc_info.value

    def test_fasttree_with_iqtree_model_errors(self, wrapper_module, monkeypatch, tmp_path, capsys):
        monkeypatch.setattr(
            wrapper_module.sys, "argv",
            ["orthophyl_pipeline_wrapper.py", "--taxon", "X",
             "--database-dir", str(tmp_path / "db"),
             "--tree-method", "fasttree", "--iqtree-model", "GTR+F+G4"],
        )
        with pytest.raises(SystemExit):
            wrapper_module.main()
        assert "iqtree-model" in capsys.readouterr().err

    def test_model_and_mset_together_errors(self, wrapper_module, monkeypatch, tmp_path, capsys):
        monkeypatch.setattr(
            wrapper_module.sys, "argv",
            ["orthophyl_pipeline_wrapper.py", "--taxon", "X",
             "--database-dir", str(tmp_path / "db"),
             "--iqtree-model", "GTR+F+G4", "--iqtree-mset", "GTR"],
        )
        with pytest.raises(SystemExit):
            wrapper_module.main()
        assert "mutually exclusive" in capsys.readouterr().err

    def test_iqtree_fast_alone_errors(self, wrapper_module, monkeypatch, tmp_path, capsys):
        monkeypatch.setattr(
            wrapper_module.sys, "argv",
            ["orthophyl_pipeline_wrapper.py", "--taxon", "X",
             "--database-dir", str(tmp_path / "db"), "--iqtree-fast"],
        )
        with pytest.raises(SystemExit):
            wrapper_module.main()
        assert "UltraFast Bootstrap" in capsys.readouterr().err

    def test_iqtree_fast_with_clean_options_passes_validation(
        self, wrapper_module, monkeypatch, tmp_path, capsys
    ):
        monkeypatch.setattr(
            wrapper_module.sys, "argv",
            ["orthophyl_pipeline_wrapper.py", "--taxon", "X",
             "--database-dir", str(tmp_path / "db"), "--dry-run",
             "--iqtree-fast", "--iqtree-options", "--symtest --lmap 2000 -t PARS --ninit 2"],
        )
        # Should get past arg validation (may exit later for other reasons,
        # e.g. missing NCBI datasets tooling in the test env, but not with
        # the UltraFast Bootstrap error).
        try:
            wrapper_module.main()
        except SystemExit:
            pass
        assert "UltraFast Bootstrap" not in capsys.readouterr().err
