"""Tests for implicit --megatree: passing a megatree-only boolean flag alone
(--megatree-lazy, --megatree-hmm-reuse, --megatree-hmm-reuse-skip-leftover)
must imply --megatree, since every internal check gates on self.megatree and
these flags were otherwise silent no-ops without it.

--placement is deliberately excluded (it routes future queries against ANY
existing megatree-shaped DB, not behavior gated on self.megatree for this
run) and is tested here as a negative case.

All external process execution is mocked -- main() is driven via sys.argv
with PipelineWrapper.run stubbed, mirroring the pattern in
test_wrapper_local_genomes.py::TestLocalGenomesMegatree.
"""

import pytest


@pytest.fixture
def Wrapper(wrapper_module):
    return wrapper_module.PipelineWrapper


def _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv):
    """Run main() with PipelineWrapper.run stubbed, capturing the
    constructed instance (specifically its self.megatree) before run() would
    have executed anything real."""
    captured = {}

    def fake_run(self):
        captured["megatree"] = self.megatree
        return 0

    monkeypatch.setattr(wrapper_module.sys, "argv", argv)
    monkeypatch.setattr(wrapper_module.PipelineWrapper, "run", fake_run)
    rc = wrapper_module.main()
    return rc, captured


class TestMegatreeImplication:
    def _base_argv(self, tmp_path):
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

    def test_megatree_lazy_alone_implies_megatree(
            self, wrapper_module, tmp_path, monkeypatch, caplog):
        argv = self._base_argv(tmp_path) + ["--megatree-lazy"]
        import logging
        with caplog.at_level(logging.INFO):
            rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["megatree"] is True
        assert any("--megatree-lazy" in rec.message and "implying --megatree"
                   in rec.message for rec in caplog.records)

    def test_megatree_hmm_reuse_alone_implies_megatree(
            self, wrapper_module, tmp_path, monkeypatch):
        argv = self._base_argv(tmp_path) + ["--megatree-hmm-reuse"]
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["megatree"] is True

    def test_megatree_hmm_reuse_skip_leftover_alone_implies_megatree(
            self, wrapper_module, tmp_path, monkeypatch):
        argv = self._base_argv(tmp_path) + ["--megatree-hmm-reuse-skip-leftover"]
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["megatree"] is True

    def test_explicit_megatree_with_lazy_still_works(
            self, wrapper_module, tmp_path, monkeypatch):
        """Passing --megatree explicitly alongside a megatree-only flag must
        still work -- the implication is a no-op when --megatree is already set."""
        argv = self._base_argv(tmp_path) + ["--megatree", "--megatree-lazy"]
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["megatree"] is True

    def test_placement_alone_does_not_imply_megatree(
            self, wrapper_module, tmp_path, monkeypatch):
        """--placement is used by routing regardless of --megatree (it
        disambiguates ties against existing megatree DBs from prior runs),
        so it must NOT force --megatree on for this run."""
        argv = self._base_argv(tmp_path) + ["--placement", "backbone"]
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["megatree"] is False

    def test_no_megatree_flags_at_all_leaves_megatree_false(
            self, wrapper_module, tmp_path, monkeypatch):
        argv = self._base_argv(tmp_path)
        rc, captured = _run_main_capturing_wrapper(wrapper_module, monkeypatch, argv)
        assert rc == 0
        assert captured["megatree"] is False

    def test_genome_dir_skip_qc_megatree_lazy_still_errors(
            self, wrapper_module, tmp_path, monkeypatch):
        """The existing --skip-qc + --megatree conflict check must still
        fire when --megatree was only implied (not passed explicitly) --
        the restriction is real regardless of how --megatree became true."""
        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()
        db_dir = tmp_path / "databases"
        db_dir.mkdir()
        argv = [
            "orthophyl_pipeline_wrapper.py",
            "--genome-dir", str(genome_dir),
            "--clade-name", "Blorptaxon",
            "--database-dir", str(db_dir),
            "--megatree-lazy",
            "--skip-qc",
        ]
        monkeypatch.setattr(wrapper_module.sys, "argv", argv)
        with pytest.raises(SystemExit):
            wrapper_module.main()
