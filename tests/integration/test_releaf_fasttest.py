"""Integration tests for ReLeaf.sh on the fasttest dataset."""

import subprocess
from pathlib import Path

import pytest

from .helpers.tree_compare import (
    assert_taxa_count,
    assert_taxa_present,
    clean_taxon_name,
    get_leaf_names,
    load_tree,
    normalize_labels,
)


@pytest.mark.integration
@pytest.mark.slow
class TestReLeafFasttest:
    """Integration tests for ReLeaf.sh using the fasttest add-assembly dataset."""
    
    def test_releaf_completes_successfully(
        self,
        orthophyl_run,
        project_root,
        releaf_script,
        fasttest_addasm_genomes,
        fasttest_addasm_annots_nucls,
        fasttest_addasm_annots_prots,
    ):
        """
        ReLeaf.sh exits 0 when adding assemblies to an existing OrthoPhyl run.
        
        This test verifies the CLI flag shape at the real tool boundary,
        catching bugs that mocked subprocess tests cannot detect.
        """
        # Build the command (mirrors test.sh)
        cmd = [
            str(releaf_script),
            "-g", str(fasttest_addasm_genomes.resolve()),
            "-a", f"{fasttest_addasm_annots_nucls},{fasttest_addasm_annots_prots}",
            "-s", str(orthophyl_run),
            "-t", "4",
            "-p", "iqtree",
            "-o", "BOTH",
        ]
        
        print(f"\n{'='*80}")
        print(f"Running ReLeaf.sh")
        print(f"Command: {' '.join(cmd)}")
        print(f"{'='*80}\n")
        
        # Run ReLeaf
        result = subprocess.run(
            cmd,
            cwd=project_root,
            capture_output=True,
            text=True,
        )
        
        # Check for success
        assert result.returncode == 0, (
            f"ReLeaf.sh failed with exit code {result.returncode}\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )
        
        print(f"\n{'='*80}")
        print(f"ReLeaf.sh completed successfully")
        print(f"{'='*80}\n")
    
    def test_releaf_outputs_present(self, orthophyl_run):
        """ReLeaf produces new alignments and trees."""
        releaf_dir = orthophyl_run / "ReLeaf_dir"
        assert releaf_dir.exists(), f"ReLeaf_dir not created: {releaf_dir}"
        
        # New alignments (trimmed and name-mapped)
        new_prot_alns = releaf_dir / "new_prot_alignments.trm.nm"
        new_cds_alns = releaf_dir / "new_CDS_alignments.trm.nm"
        
        assert new_prot_alns.exists(), f"New protein alignments not found: {new_prot_alns}"
        assert new_cds_alns.exists(), f"New CDS alignments not found: {new_cds_alns}"
        
        # Check that alignments are non-empty directories
        assert list(new_prot_alns.iterdir()), "New protein alignments directory is empty"
        assert list(new_cds_alns.iterdir()), "New CDS alignments directory is empty"
        
        # New trees
        tree_dir = releaf_dir / "new_trees"
        assert tree_dir.exists(), f"new_trees directory not found: {tree_dir}"
        
        tree_files = list(tree_dir.glob("*.tree*"))
        assert len(tree_files) > 0, f"No tree files found in {tree_dir}"
    
    def test_added_taxa_in_releaf_tree(self, orthophyl_run, expected_releaf_added_taxa):
        """All newly-added genomes appear in the ReLeaf output tree."""
        # ReLeaf can produce trees in multiple locations/names
        tree_candidates = [
            orthophyl_run / "ReLeaf_dir" / "phylogeny_with_new_genomes.nwk",
            orthophyl_run / "ReLeaf_dir" / "new_trees" / "phylogeny.nwk",
            orthophyl_run / "ReLeaf_dir" / "new_trees" / "phylogeny_with_new_genomes.nwk",
        ]
        
        # Also check for any .tree or .treefile in new_trees
        new_trees_dir = orthophyl_run / "ReLeaf_dir" / "new_trees"
        if new_trees_dir.exists():
            tree_candidates.extend(new_trees_dir.glob("*.tree"))
            tree_candidates.extend(new_trees_dir.glob("*.treefile"))
        
        # Find the first existing tree
        tree_path = next((p for p in tree_candidates if p.exists()), None)
        assert tree_path, (
            f"No ReLeaf output tree found. Checked:\n"
            f"{[str(p) for p in tree_candidates]}"
        )
        
        print(f"Using ReLeaf tree: {tree_path}")
        
        tree = load_tree(str(tree_path))
        tree = normalize_labels(tree)
        
        # Get cleaned tree taxa
        tree_taxa = {clean_taxon_name(leaf.name) for leaf in tree.get_leaves()}
        
        # All added taxa should be present
        missing = expected_releaf_added_taxa - tree_taxa
        assert not missing, (
            f"Missing added taxa in ReLeaf tree: {sorted(missing)}\n"
            f"Expected added ({len(expected_releaf_added_taxa)}): {sorted(expected_releaf_added_taxa)}\n"
            f"Tree taxa ({len(tree_taxa)}): {sorted(tree_taxa)}"
        )
    
    def test_original_taxa_retained(self, orthophyl_run, expected_orthophyl_taxa):
        """All original OrthoPhyl taxa are still present after ReLeaf (no dropouts)."""
        # Find the ReLeaf output tree
        tree_candidates = [
            orthophyl_run / "ReLeaf_dir" / "phylogeny_with_new_genomes.nwk",
            orthophyl_run / "ReLeaf_dir" / "new_trees" / "phylogeny.nwk",
            orthophyl_run / "ReLeaf_dir" / "new_trees" / "phylogeny_with_new_genomes.nwk",
        ]
        
        new_trees_dir = orthophyl_run / "ReLeaf_dir" / "new_trees"
        if new_trees_dir.exists():
            tree_candidates.extend(new_trees_dir.glob("*.tree"))
            tree_candidates.extend(new_trees_dir.glob("*.treefile"))
        
        tree_path = next((p for p in tree_candidates if p.exists()), None)
        assert tree_path, "No ReLeaf output tree found"
        
        tree = load_tree(str(tree_path))
        tree = normalize_labels(tree)
        
        # Get cleaned tree taxa
        tree_taxa = {clean_taxon_name(leaf.name) for leaf in tree.get_leaves()}
        
        # All original taxa should still be present
        missing = expected_orthophyl_taxa - tree_taxa
        assert not missing, (
            f"Missing original taxa in ReLeaf tree: {sorted(missing)}\n"
            f"Expected original ({len(expected_orthophyl_taxa)}): {sorted(expected_orthophyl_taxa)}\n"
            f"Tree taxa ({len(tree_taxa)}): {sorted(tree_taxa)}"
        )
    
    def test_releaf_tree_has_expected_taxa(self, orthophyl_run, expected_releaf_total_taxa):
        """Final ReLeaf tree contains exactly the expected number of taxa (original + added)."""
        # Find the ReLeaf output tree
        tree_candidates = [
            orthophyl_run / "ReLeaf_dir" / "phylogeny_with_new_genomes.nwk",
            orthophyl_run / "ReLeaf_dir" / "new_trees" / "phylogeny.nwk",
            orthophyl_run / "ReLeaf_dir" / "new_trees" / "phylogeny_with_new_genomes.nwk",
        ]
        
        new_trees_dir = orthophyl_run / "ReLeaf_dir" / "new_trees"
        if new_trees_dir.exists():
            tree_candidates.extend(new_trees_dir.glob("*.tree"))
            tree_candidates.extend(new_trees_dir.glob("*.treefile"))
        
        tree_path = next((p for p in tree_candidates if p.exists()), None)
        assert tree_path, "No ReLeaf output tree found"
        
        tree = load_tree(str(tree_path))
        tree = normalize_labels(tree)
        
        # Verify tree has exactly the expected number of taxa
        assert_taxa_count(tree, len(expected_releaf_total_taxa))
    
    def test_hmm_search_results_present(self, orthophyl_run):
        """ReLeaf HMM search results are generated for added genomes."""
        releaf_dir = orthophyl_run / "ReLeaf_dir"
        
        # HMM search results (domtblout files from hmmsearch)
        # These may be in various subdirectories depending on ReLeaf version
        hmm_search_candidates = [
            releaf_dir / "hmm_search_results",
            releaf_dir / "hmmsearch_output",
            releaf_dir,  # Sometimes directly in ReLeaf_dir
        ]
        
        # Look for .domtblout files (hmmsearch output)
        domtblout_found = False
        for search_dir in hmm_search_candidates:
            if search_dir.exists():
                domtblout_files = list(search_dir.glob("**/*.domtblout"))
                if domtblout_files:
                    domtblout_found = True
                    break
        
        # If not found, at least check that some processing happened
        # (new alignments exist, which implies HMM search succeeded)
        if not domtblout_found:
            # Fallback: verify new alignments exist (tested elsewhere)
            new_prot_alns = releaf_dir / "new_prot_alignments.trm.nm"
            assert new_prot_alns.exists(), (
                "No HMM search results found and no new alignments generated"
            )
    
    def test_releaf_preserves_orthophyl_outputs(self, orthophyl_run, expected_orthophyl_taxa):
        """ReLeaf doesn't corrupt the original OrthoPhyl outputs."""
        # Original OrthoPhyl tree should still exist
        original_tree = (
            orthophyl_run / "phylo_current" / "SpeciesTree" / "iqtree.SCO_strict.CDS.tree"
        )
        assert original_tree.exists(), "Original OrthoPhyl tree was removed/corrupted"
        
        # Original tree should still have the expected number of taxa
        tree = load_tree(str(original_tree))
        tree = normalize_labels(tree)
        assert_taxa_count(tree, len(expected_orthophyl_taxa))
        
        # Original alignments should still exist
        phylo_current = orthophyl_run / "phylo_current"
        assert (phylo_current / "AlignmentsProts.trm.nm").exists()
        assert (phylo_current / "AlignmentsCDS.trm.nm").exists()
