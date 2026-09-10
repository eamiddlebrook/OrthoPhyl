"""Tests for gather_filter_asms.sh command construction in _download_genomes.

These assert the wrapper builds the correct argv for the genome-retention
controls added alongside CheckM2 QC of query genomes:

  * --must-keep <value>          (comma-sep list or a file path, passed raw)
  * --query-genomes <paths>      (query genomes QC'd alongside downloads)
  * --keep-failing-query         (downgrade query-QC failure to a warning)

All subprocess execution is mocked -- the real gather script is never run.
"""

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


def _make_wrapper(Wrapper, tmp_path, **overrides):
    """Wrapper with a fake gather script, ready for _download_genomes."""
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
    )
    kwargs.update(overrides)
    w = Wrapper(**kwargs)
    # _download_genomes writes a per-taxon log into logs_dir (normally created
    # during _phase_initialization, which we skip here).
    w.logs_dir.mkdir(parents=True, exist_ok=True)
    return w


def _prep_genomes_to_keep(output_dir):
    """Create a genomes_to_keep/ + stats file so post-download checks pass."""
    gtk = output_dir / "genomes_to_keep"
    gtk.mkdir(parents=True, exist_ok=True)
    (gtk / "GCF_000001.1.fna").write_text(">x\nATCG\n")
    stats = output_dir / "assemblies_all.stats.txt"
    stats.write_text("header\nGCF_000001.1\tcheckm2\n")


def _get_gather_cmd(recording_run):
    for c in recording_run:
        if "fake_gather.sh" in " ".join(str(x) for x in c["cmd"]):
            return [str(x) for x in c["cmd"]]
    raise AssertionError("gather script was not invoked")


def test_must_keep_passed_through(Wrapper, tmp_path, recording_run):
    w = _make_wrapper(Wrapper, tmp_path, must_keep="GCF_000001.1,GCF_000002.1")
    out = tmp_path / "dl"
    _prep_genomes_to_keep(out)
    w._download_genomes("SomeTaxon", out)
    cmd = _get_gather_cmd(recording_run)
    assert "--must-keep" in cmd
    assert cmd[cmd.index("--must-keep") + 1] == "GCF_000001.1,GCF_000002.1"


def test_query_genomes_joined_paths(Wrapper, tmp_path, recording_run):
    w = _make_wrapper(Wrapper, tmp_path)
    out = tmp_path / "dl"
    _prep_genomes_to_keep(out)
    q1 = tmp_path / "q1.fasta"
    q2 = tmp_path / "q2.fna"
    q1.write_text(">a\nAC\n")
    q2.write_text(">b\nGT\n")
    assemblies = [
        {"assembly_id": "q1", "assembly_path": str(q1)},
        {"assembly_id": "q2", "assembly_path": str(q2)},
    ]
    w._download_genomes("SomeTaxon", out, query_assemblies=assemblies)
    cmd = _get_gather_cmd(recording_run)
    assert "--query-genomes" in cmd
    joined = cmd[cmd.index("--query-genomes") + 1]
    assert joined == f"{q1},{q2}"


def test_keep_failing_query_flag(Wrapper, tmp_path, recording_run):
    w = _make_wrapper(Wrapper, tmp_path, keep_failing_query=True)
    out = tmp_path / "dl"
    _prep_genomes_to_keep(out)
    w._download_genomes("SomeTaxon", out)
    cmd = _get_gather_cmd(recording_run)
    assert "--keep-failing-query" in cmd


def test_no_retention_flags_by_default(Wrapper, tmp_path, recording_run):
    w = _make_wrapper(Wrapper, tmp_path)
    out = tmp_path / "dl"
    _prep_genomes_to_keep(out)
    w._download_genomes("SomeTaxon", out)
    cmd = _get_gather_cmd(recording_run)
    assert "--must-keep" not in cmd
    assert "--query-genomes" not in cmd
    assert "--keep-failing-query" not in cmd
