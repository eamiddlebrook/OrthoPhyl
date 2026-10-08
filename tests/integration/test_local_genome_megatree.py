"""Integration test for orthophyl_pipeline_wrapper.py --genome-dir --megatree,
on real chloroplast data.

Reproduces (and pins as a regression test) a manual verification run during
development of this feature: --genome-dir local-ingest mode has no query
assemblies (the whole --genome-dir IS the input), so an oversized local
genome set takes the same full-coverage partition/backbone/graft path
--taxon create mode already uses in that scenario, instead of being
rejected outright. Confirms the resulting subclade/backbone databases
correctly carry taxonomy_source="user_supplied"/qc_applied=true, and that
the merged megatree + conflict report are published.
"""

import json
import subprocess
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.slow
class TestLocalGenomeMegatreeEndToEnd:
    """Drives orthophyl_pipeline_wrapper.py as a subprocess (not mocked) --
    catches bugs mocked-subprocess unit tests structurally cannot."""

    @pytest.fixture(scope="class")
    def database_dir(self, tmp_path_factory):
        return tmp_path_factory.mktemp("local_megatree_db")

    @pytest.fixture(scope="class")
    def output_dir(self, tmp_path_factory):
        return tmp_path_factory.mktemp("local_megatree_out")

    @pytest.fixture(scope="class")
    def run_result(
        self, project_root, wrapper_script, gather_script, database_dir,
        output_dir, local_megatree_genomes,
    ):
        """Run the wrapper once per class: 20 genomes, --max-tree-genomes 8
        forces the megatree branch (20 > 8), --ani-shortlist 3 forces
        OrthoPhyl.sh's MASH-shortlist/HMM-building branch for both the
        subclade and backbone builds, --backbone-reps 5 keeps IQ-TREE's
        bootstrap step above its own >=4-sequence floor."""
        cmd = [
            "python3", str(wrapper_script),
            "--genome-dir", str(local_megatree_genomes),
            "--clade-name", "LocalMegatreeTestClade",
            "--database-dir", str(database_dir),
            "--output-dir", str(output_dir),
            "--gather-script", str(gather_script),
            "--use-bbmap",
            "--megatree",
            "--max-tree-genomes", "8",
            "--subclade-size", "8",
            "--backbone-reps", "5",
            "--ani-shortlist", "3",
            "--threads", "4",
        ]
        result = subprocess.run(
            cmd, cwd=project_root, capture_output=True, text=True, timeout=1700)
        return result

    def test_run_completes_successfully(self, run_result):
        assert run_result.returncode == 0, (
            f"--genome-dir --megatree run failed:\n"
            f"STDOUT:\n{run_result.stdout}\n\nSTDERR:\n{run_result.stderr}")

    def test_merged_megatree_published(self, run_result, output_dir):
        tree_dir = output_dir / "03_results" / "trees" / "orthophyl"
        merged_tree = tree_dir / "LocalMegatreeTestClade_megatree.nwk"
        conflict_report = tree_dir / "LocalMegatreeTestClade_megatree_conflicts.json"
        assert merged_tree.exists(), (
            f"Merged megatree not found at {merged_tree}.\n"
            f"STDOUT:\n{run_result.stdout}")
        assert merged_tree.stat().st_size > 0
        assert conflict_report.exists(), (
            f"Conflict report not found at {conflict_report}")
        # Must be valid JSON (a list, possibly empty).
        json.loads(conflict_report.read_text())

    def test_backbone_db_has_user_supplied_provenance(self, run_result, database_dir):
        config_path = database_dir / "LocalMegatreeTestClade_db" / "database_config.json"
        assert config_path.exists(), (
            f"Backbone database not created at {config_path}.\n"
            f"STDOUT:\n{run_result.stdout}")
        config = json.loads(config_path.read_text())
        assert config.get("taxonomy_source") == "user_supplied", (
            f"Expected taxonomy_source='user_supplied' on the backbone DB, "
            f"got: {config.get('taxonomy_source')!r}")
        assert config.get("qc_applied") is True
        assert config.get("is_backbone") is True
        assert config.get("parent_taxon") == "LocalMegatreeTestClade"

    def test_subclade_dbs_have_user_supplied_provenance(self, run_result, database_dir):
        subclade_dbs = sorted(database_dir.glob("LocalMegatreeTestClade_*_db"))
        assert subclade_dbs, (
            f"No subclade databases found under {database_dir}.\n"
            f"STDOUT:\n{run_result.stdout}")
        for db_dir in subclade_dbs:
            config = json.loads((db_dir / "database_config.json").read_text())
            assert config.get("taxonomy_source") == "user_supplied", (
                f"{db_dir.name}: expected taxonomy_source='user_supplied', "
                f"got {config.get('taxonomy_source')!r}")
            assert config.get("qc_applied") is True, (
                f"{db_dir.name}: expected qc_applied=True")
            assert config.get("is_subclade") is True
            assert config.get("parent_taxon") == "LocalMegatreeTestClade"
            assert config.get("built") is True, (
                f"{db_dir.name}: subclade DB should be fully built (no query -> "
                f"every subclade built eagerly, not lazily registered)")
