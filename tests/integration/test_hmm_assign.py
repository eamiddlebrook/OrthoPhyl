"""Integration tests for OrthoPhyl.sh --hmm-assign-dir, on real chloroplast data.

Reproduces (and pins as regression tests) manual verification scenarios run
during development of dev_plans/near_term_improvements.md items 3 and 4:
assigning genes into a precomputed external HMM set instead of running
OrthoFinder from scratch, and the now-default behavior of routing genes that
match none of those HMMs through a real OrthoFinder clustering pass of their
own (--hmm-assign-leftover-orthofinder, default on; --skip-hmm-assign-leftover
restores the old log-and-drop behavior).

Drives OrthoPhyl.sh directly (bash level), not through
orthophyl_pipeline_wrapper.py -- --hmm-assign-dir has no wrapper-level flag
(deliberate, see dev_plans/near_term_improvements.md item 4's "Non-goals").
"""

import subprocess
from pathlib import Path

import pytest


def _run_orthophyl(orthophyl_script, project_root, control_file, args, timeout=900):
    cmd = [str(orthophyl_script), *args]
    return subprocess.run(
        cmd, cwd=project_root, capture_output=True, text=True, timeout=timeout)


@pytest.mark.integration
@pytest.mark.slow
class TestHmmAssignDirBasic:
    """--hmm-assign-dir with --skip-hmm-assign-leftover: genes matching no
    external HMM are logged and dropped (the pre-leftover-routing behavior),
    pinning the plain external-assignment-only path."""

    @pytest.fixture(scope="class")
    def backbone_dir(self, tmp_path_factory, project_root, orthophyl_script,
                      control_file, hmm_assign_backbone_genomes):
        """Build a real 8-genome backbone, forced through the MASH-shortlist/
        HMM-building branch via -n 3 (8 genomes > 3), exactly as manually
        verified during development."""
        output_dir = tmp_path_factory.mktemp("hmm_assign_backbone_out")
        result = _run_orthophyl(
            orthophyl_script, project_root, control_file,
            [
                "-g", str(hmm_assign_backbone_genomes),
                "-s", str(output_dir),
                "-t", "4",
                "-c", str(control_file),
                "-n", "3",
            ],
        )
        assert result.returncode == 0, (
            f"Backbone build failed:\nSTDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}")
        return output_dir

    @pytest.fixture(scope="class")
    def backbone_hmm_dir(self, backbone_dir):
        # OG_alignmentsToHMM lives directly under $store (the -s output
        # dir), NOT under phylo_current ($wd) -- confirmed by inspecting a
        # real build's output tree.
        hmm_dir = backbone_dir / "OG_alignmentsToHMM" / "hmms_final"
        assert hmm_dir.exists() and list(hmm_dir.glob("*.hmm")), (
            f"Backbone did not produce hmms_final/ at {hmm_dir}")
        return hmm_dir

    @pytest.fixture(scope="class")
    def subclade_result(self, tmp_path_factory, project_root, orthophyl_script,
                         control_file, hmm_assign_subclade_genomes,
                         backbone_hmm_dir):
        output_dir = tmp_path_factory.mktemp("hmm_assign_basic_subclade_out")
        result = _run_orthophyl(
            orthophyl_script, project_root, control_file,
            [
                "-g", str(hmm_assign_subclade_genomes),
                "-s", str(output_dir),
                "-t", "4",
                "-c", str(control_file),
                "--hmm-assign-dir", str(backbone_hmm_dir),
                "--skip-hmm-assign-leftover",
            ],
        )
        return result, output_dir

    def test_subclade_build_succeeds(self, subclade_result):
        result, _ = subclade_result
        assert result.returncode == 0, (
            f"Subclade build failed:\nSTDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}")

    def test_no_orthofinder_run_in_subclade(self, subclade_result):
        """--hmm-assign-dir must skip OrthoFinder/MASH-shortlisting
        entirely for the externally-assigned genes."""
        _, output_dir = subclade_result
        results_dirs = list(output_dir.rglob("Results_*"))
        assert not results_dirs, (
            f"Found unexpected OrthoFinder Results_* dir(s) under a pure "
            f"--hmm-assign-dir (no leftover) run: {results_dirs}")

    def test_assigned_ogs_are_subset_of_backbone_hmms(
            self, subclade_result, backbone_hmm_dir):
        _, output_dir = subclade_result
        alignments_dir = output_dir / "phylo_current" / "AlignmentsProts"
        assert alignments_dir.exists(), f"{alignments_dir} not found"
        backbone_og_ids = {f.stem for f in backbone_hmm_dir.glob("*.hmm")}
        assigned_og_ids = {f.stem for f in alignments_dir.glob("*.faa")}
        assert assigned_og_ids, "No OGs assigned in subclade AlignmentsProts/"
        foreign = assigned_og_ids - backbone_og_ids
        assert not foreign, (
            f"Found OG IDs in subclade AlignmentsProts/ not present in the "
            f"backbone's HMM set (unexpected with --skip-hmm-assign-leftover): "
            f"{foreign}")
        lft_ids = {og for og in assigned_og_ids if og.startswith("OG0_LFT_")}
        assert not lft_ids, (
            f"Found OG0_LFT_* leftover-routed OGs despite "
            f"--skip-hmm-assign-leftover: {lft_ids}")

    def test_unmatched_ogs_logged(self, subclade_result):
        _, output_dir = subclade_result
        unmatched_log = output_dir / "phylo_current" / "hmm_assign_unmatched_OGs.txt"
        assert unmatched_log.exists(), (
            f"hmm_assign_unmatched_OGs.txt not found at {unmatched_log}")

    def test_species_tree_produced(self, subclade_result):
        _, output_dir = subclade_result
        tree_dir = output_dir / "phylo_current" / "SpeciesTree"
        assert tree_dir.exists() and list(tree_dir.glob("*.tree")), (
            f"No species tree produced under {tree_dir}")


@pytest.mark.integration
@pytest.mark.slow
class TestHmmAssignLeftoverDefault:
    """--hmm-assign-dir with a deliberately PARTIAL external HMM set and no
    extra flag: leftover-OrthoFinder routing is on by default, so genes
    matching none of the (partial) external HMMs get pooled and clustered
    into additional OG0_LFT_* orthogroups."""

    @pytest.fixture(scope="class")
    def backbone_dir(self, tmp_path_factory, project_root, orthophyl_script,
                      control_file, hmm_assign_backbone_genomes):
        output_dir = tmp_path_factory.mktemp("hmm_assign_leftover_backbone_out")
        result = _run_orthophyl(
            orthophyl_script, project_root, control_file,
            [
                "-g", str(hmm_assign_backbone_genomes),
                "-s", str(output_dir),
                "-t", "4",
                "-c", str(control_file),
                "-n", "3",
            ],
        )
        assert result.returncode == 0, (
            f"Backbone build failed:\nSTDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}")
        return output_dir

    @pytest.fixture(scope="class")
    def partial_hmm_dir(self, backbone_dir, tmp_path_factory):
        """Take only the first 40 of the backbone's HMMs -- a deliberately
        incomplete external set, mirroring the manual verification -- so a
        meaningful number of genes have no external HMM to match."""
        full_hmm_dir = backbone_dir / "OG_alignmentsToHMM" / "hmms_final"
        all_hmms = sorted(full_hmm_dir.glob("*.hmm"))
        assert len(all_hmms) > 40, (
            f"Expected > 40 backbone HMMs to take a meaningful partial slice, "
            f"found {len(all_hmms)}")
        partial_dir = tmp_path_factory.mktemp("hmm_assign_partial_hmms")
        for hmm_file in all_hmms[:40]:
            (partial_dir / hmm_file.name).write_bytes(hmm_file.read_bytes())
        return partial_dir

    @pytest.fixture(scope="class")
    def subclade_result(self, tmp_path_factory, project_root, orthophyl_script,
                         control_file, hmm_assign_subclade_genomes,
                         partial_hmm_dir):
        output_dir = tmp_path_factory.mktemp("hmm_assign_leftover_subclade_out")
        result = _run_orthophyl(
            orthophyl_script, project_root, control_file,
            [
                "-g", str(hmm_assign_subclade_genomes),
                "-s", str(output_dir),
                "-t", "4",
                "-c", str(control_file),
                "--hmm-assign-dir", str(partial_hmm_dir),
            ],
        )
        return result, output_dir

    def test_subclade_build_succeeds(self, subclade_result):
        result, _ = subclade_result
        assert result.returncode == 0, (
            f"Subclade build failed:\nSTDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}")

    def test_leftover_prots_materialized(self, subclade_result):
        _, output_dir = subclade_result
        leftover_dir = output_dir / "phylo_current" / "hmm_assign_leftover_prots"
        assert leftover_dir.exists(), (
            f"hmm_assign_leftover_prots/ not found at {leftover_dir} -- "
            f"expected with a partial external HMM set and default "
            f"leftover-routing")
        assert list(leftover_dir.glob("*.faa")), (
            "hmm_assign_leftover_prots/ exists but has no per-genome .faa files")

    def test_real_orthofinder_ran_on_leftovers(self, subclade_result):
        _, output_dir = subclade_result
        leftover_dir = output_dir / "phylo_current" / "hmm_assign_leftover_prots"
        results_dirs = list(leftover_dir.rglob("Results_*"))
        assert results_dirs, (
            f"No OrthoFinder Results_* dir found under {leftover_dir} -- "
            f"leftover genes should be routed through a real ORTHO_RUN")

    def test_both_namespaces_present_no_collision(self, subclade_result):
        _, output_dir = subclade_result
        alignments_dir = output_dir / "phylo_current" / "AlignmentsProts"
        assert alignments_dir.exists(), f"{alignments_dir} not found"
        all_ids = {f.stem for f in alignments_dir.glob("*.faa")}
        lft_ids = {og for og in all_ids if og.startswith("OG0_LFT_")}
        external_ids = all_ids - lft_ids
        assert external_ids, "No externally-assigned OGs found"
        assert lft_ids, (
            "No OG0_LFT_* leftover-routed OGs found -- expected with a "
            "partial external HMM set and default leftover-routing")
        # By construction (disjoint prefixes), there can be no collision --
        # this assertion documents the invariant rather than detects one.
        assert external_ids.isdisjoint(lft_ids)

    def test_species_tree_produced(self, subclade_result):
        _, output_dir = subclade_result
        tree_dir = output_dir / "phylo_current" / "SpeciesTree"
        assert tree_dir.exists() and list(tree_dir.glob("*.tree")), (
            f"No species tree produced under {tree_dir}")


@pytest.mark.integration
@pytest.mark.slow
class TestSkipHmmAssignLeftover:
    """Identical setup to TestHmmAssignLeftoverDefault's partial-HMM scenario,
    but with --skip-hmm-assign-leftover: restores the old log-and-drop
    behavior even with a partial external HMM set."""

    @pytest.fixture(scope="class")
    def backbone_dir(self, tmp_path_factory, project_root, orthophyl_script,
                      control_file, hmm_assign_backbone_genomes):
        output_dir = tmp_path_factory.mktemp("hmm_assign_skip_backbone_out")
        result = _run_orthophyl(
            orthophyl_script, project_root, control_file,
            [
                "-g", str(hmm_assign_backbone_genomes),
                "-s", str(output_dir),
                "-t", "4",
                "-c", str(control_file),
                "-n", "3",
            ],
        )
        assert result.returncode == 0, (
            f"Backbone build failed:\nSTDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}")
        return output_dir

    @pytest.fixture(scope="class")
    def partial_hmm_dir(self, backbone_dir, tmp_path_factory):
        full_hmm_dir = backbone_dir / "OG_alignmentsToHMM" / "hmms_final"
        all_hmms = sorted(full_hmm_dir.glob("*.hmm"))
        assert len(all_hmms) > 40, (
            f"Expected > 40 backbone HMMs to take a meaningful partial slice, "
            f"found {len(all_hmms)}")
        partial_dir = tmp_path_factory.mktemp("hmm_assign_skip_partial_hmms")
        for hmm_file in all_hmms[:40]:
            (partial_dir / hmm_file.name).write_bytes(hmm_file.read_bytes())
        return partial_dir

    @pytest.fixture(scope="class")
    def subclade_result(self, tmp_path_factory, project_root, orthophyl_script,
                         control_file, hmm_assign_subclade_genomes,
                         partial_hmm_dir):
        output_dir = tmp_path_factory.mktemp("hmm_assign_skip_subclade_out")
        result = _run_orthophyl(
            orthophyl_script, project_root, control_file,
            [
                "-g", str(hmm_assign_subclade_genomes),
                "-s", str(output_dir),
                "-t", "4",
                "-c", str(control_file),
                "--hmm-assign-dir", str(partial_hmm_dir),
                "--skip-hmm-assign-leftover",
            ],
        )
        return result, output_dir

    def test_subclade_build_succeeds(self, subclade_result):
        result, _ = subclade_result
        assert result.returncode == 0, (
            f"Subclade build failed:\nSTDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}")

    def test_no_leftover_prots_directory(self, subclade_result):
        _, output_dir = subclade_result
        leftover_dir = output_dir / "phylo_current" / "hmm_assign_leftover_prots"
        assert not leftover_dir.exists(), (
            f"hmm_assign_leftover_prots/ should not exist with "
            f"--skip-hmm-assign-leftover, found at {leftover_dir}")

    def test_no_leftover_ogs_in_alignments(self, subclade_result):
        _, output_dir = subclade_result
        alignments_dir = output_dir / "phylo_current" / "AlignmentsProts"
        assert alignments_dir.exists(), f"{alignments_dir} not found"
        lft_ids = {f.stem for f in alignments_dir.glob("*.faa")
                   if f.stem.startswith("OG0_LFT_")}
        assert not lft_ids, (
            f"Found OG0_LFT_* entries despite --skip-hmm-assign-leftover: "
            f"{lft_ids}")

    def test_unmatched_ogs_still_logged(self, subclade_result):
        _, output_dir = subclade_result
        unmatched_log = output_dir / "phylo_current" / "hmm_assign_unmatched_OGs.txt"
        assert unmatched_log.exists(), (
            f"hmm_assign_unmatched_OGs.txt not found at {unmatched_log} -- "
            f"the OG-level audit log should still be written")
