"""Integration tests for OrthoPhyl.sh on the fasttest dataset."""

from pathlib import Path

import pytest

from .helpers.tree_compare import (
    assert_taxa_count,
    assert_taxa_present,
    clean_taxon_name,
    load_tree,
    normalize_labels,
    same_topology,
)


@pytest.mark.integration
@pytest.mark.slow
class TestOrthoPhylFasttest:
    """Integration tests for OrthoPhyl.sh using the fasttest (orchid chloroplast) dataset."""
    
    def test_orthophyl_completes_successfully(self, orthophyl_run):
        """OrthoPhyl.sh exits 0 and produces expected directory structure."""
        assert orthophyl_run.exists(), f"Output directory not created: {orthophyl_run}"
        assert (orthophyl_run / "genome_list").exists(), "genome_list not found"
        assert (orthophyl_run / "phylo_current").exists(), "phylo_current directory not found"
    
    def test_species_trees_generated(self, orthophyl_run):
        """Final species trees exist for both SCO_strict and SCO_3."""
        tree_dir = orthophyl_run / "phylo_current" / "SpeciesTree"
        assert tree_dir.exists(), f"SpeciesTree directory not found: {tree_dir}"
        
        # Top-level tree copies (these are the main output trees)
        strict_tree = tree_dir / "iqtree.SCO_strict.CDS.tree"
        sco3_tree = tree_dir / "iqtree.SCO_3.CDS.tree"
        
        assert strict_tree.exists(), f"SCO_strict tree not found: {strict_tree}"
        assert sco3_tree.exists(), f"SCO_3 tree not found: {sco3_tree}"
        
        # Detailed iqtree outputs
        iqtree_dir = tree_dir / "iqtree"
        assert iqtree_dir.exists(), f"iqtree subdirectory not found: {iqtree_dir}"
        
        treefile = iqtree_dir / "iqtree.SCO_strict.CDS.treefile"
        iqtree_report = iqtree_dir / "iqtree.SCO_strict.CDS.iqtree"
        
        assert treefile.exists(), f"IQ-TREE treefile not found: {treefile}"
        assert iqtree_report.exists(), f"IQ-TREE report not found: {iqtree_report}"
    
    def test_hmms_and_alignments_present(self, orthophyl_run):
        """HMMs and trimmed alignments are generated."""
        # HMMs can be in multiple locations depending on OrthoPhyl version
        hmm_candidates = [
            orthophyl_run / "OG_alignmentsToHMM" / "hmms_final",
            orthophyl_run / "phylo_current" / "OG_alignmentsToHMM" / "hmms_final",
        ]
        
        hmm_found = False
        for hmm_dir in hmm_candidates:
            if hmm_dir.exists() and list(hmm_dir.glob("*.hmm")):
                hmm_found = True
                break
        
        assert hmm_found, (
            f"No HMMs found in any of the expected locations:\n"
            f"{[str(d) for d in hmm_candidates]}"
        )
        
        # Trimmed alignments (used by ReLeaf)
        phylo_current = orthophyl_run / "phylo_current"
        
        prot_alns = phylo_current / "AlignmentsProts.trm.nm"
        cds_alns = phylo_current / "AlignmentsCDS.trm.nm"
        
        assert prot_alns.exists(), f"Protein alignments not found: {prot_alns}"
        assert cds_alns.exists(), f"CDS alignments not found: {cds_alns}"
        
        # Check that alignments are non-empty directories
        assert list(prot_alns.iterdir()), "Protein alignments directory is empty"
        assert list(cds_alns.iterdir()), "CDS alignments directory is empty"
    
    def test_all_input_taxa_in_tree(self, orthophyl_run):
        """All input taxa (genomes + pre-annotated) appear in the final species tree."""
        tree_path = orthophyl_run / "phylo_current" / "SpeciesTree" / "iqtree.SCO_strict.CDS.tree"
        tree = load_tree(str(tree_path))
        tree = normalize_labels(tree)
        
        # Read actual input taxa from genome_list and pre_annotated_list
        genome_list_path = orthophyl_run / "genome_list"
        pre_annotated_path = orthophyl_run / "pre_annotated_list"
        
        expected_taxa = set()
        
        # Read genome_list (skip comments and blank lines)
        if genome_list_path.exists():
            with open(genome_list_path) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        # Remove .fasta extension if present
                        taxon = line.replace(".fasta", "").replace(".fa", "")
                        expected_taxa.add(taxon)
        
        # Read pre_annotated_list (skip comments and blank lines)
        if pre_annotated_path.exists():
            with open(pre_annotated_path) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        # Remove .faa extension if present
                        taxon = line.replace(".faa", "").replace(".fa", "")
                        expected_taxa.add(taxon)
        
        # With -n 5, fasttest (9 inputs) goes through MASH shortlist
        # All 9 inputs should make it to the final tree
        assert len(expected_taxa) == 9, f"Expected 9 input taxa, found {len(expected_taxa)}"
        
        # Check that all inputs made it to the tree
        tree_taxa = {leaf.name for leaf in tree.get_leaves()}
        present_count = sum(1 for taxon in expected_taxa if taxon in tree_taxa)
        
        assert present_count == 9, (
            f"Expected all 9 input taxa in tree, found {present_count}\n"
            f"Input taxa: {expected_taxa}\n"
            f"Tree taxa: {tree_taxa}"
        )
        
        # Verify tree has exactly 9 taxa
        assert_taxa_count(tree, 9)
    
    def test_source_inputs_match_tree(self, orthophyl_run, expected_orthophyl_taxa):
        """
        All source input taxa (from input directories) appear in the final tree.
        
        This test validates against the SOURCE input directories rather than
        OrthoPhyl's internal bookkeeping files. It catches bugs where inputs
        are dropped during the copy/clean phase.
        """
        tree_path = orthophyl_run / "phylo_current" / "SpeciesTree" / "iqtree.SCO_strict.CDS.tree"
        tree = load_tree(str(tree_path))
        tree = normalize_labels(tree)
        
        # Get cleaned tree taxa (apply same normalization as fixture)
        tree_taxa = {clean_taxon_name(leaf.name) for leaf in tree.get_leaves()}
        
        # All expected taxa from source directories should be present
        missing = expected_orthophyl_taxa - tree_taxa
        assert not missing, (
            f"Missing taxa in tree: {sorted(missing)}\n"
            f"Expected taxa from source dirs ({len(expected_orthophyl_taxa)}): {sorted(expected_orthophyl_taxa)}\n"
            f"Tree taxa ({len(tree_taxa)}): {sorted(tree_taxa)}"
        )
        
        # Verify tree has exactly the expected number of taxa
        assert_taxa_count(tree, len(expected_orthophyl_taxa))
    
    @pytest.mark.skip(reason="Reference tree needs to be generated from a successful -n 5 run")
    def test_topology_matches_reference(self, orthophyl_run, reference_trees):
        """
        Produced tree has exact topology (RF=0) vs reference tree.
        
        This is the key test — it verifies not just that a tree was produced,
        but that the *correct* tree was produced. If the pipeline silently
        uses wrong parameters, drops taxa, or has a tool integration bug,
        the topology will diverge.
        
        NOTE: Skipped until we generate a reference tree from a successful
        -n 5 run (which triggers MASH shortlist and produces 8-taxon trees).
        """
        produced_tree_path = (
            orthophyl_run / "phylo_current" / "SpeciesTree" / "iqtree.SCO_strict.CDS.tree"
        )
        reference_tree_path = reference_trees / "fasttest_n5.SCO_strict.CDS.tree"
        
        # Load and normalize both trees
        produced = load_tree(str(produced_tree_path))
        reference = load_tree(str(reference_tree_path))
        
        produced = normalize_labels(produced)
        reference = normalize_labels(reference)
        
        # Check topology match (RF distance = 0)
        assert same_topology(produced, reference, tolerance=0), (
            f"Tree topology differs from reference (RF > 0)\n"
            f"Produced: {produced_tree_path}\n"
            f"Reference: {reference_tree_path}\n"
            f"This indicates a potential pipeline bug or parameter change."
        )
    
    def test_genome_list_matches_input(self, orthophyl_run, fasttest_genomes):
        """The genome_list file contains all input genomes."""
        genome_list = orthophyl_run / "genome_list"
        assert genome_list.exists(), "genome_list not found"
        
        # Read genome list (skip comments and blank lines), clean names
        with open(genome_list) as f:
            listed_genomes = {
                clean_taxon_name(line.strip()) for line in f
                if line.strip() and not line.startswith("#")
            }
        
        # Get input genome basenames (without extensions), clean names
        input_genomes = {
            clean_taxon_name(path.stem) for path in fasttest_genomes.glob("*.fasta")
        }
        
        # genome_list should contain all input genomes
        # (it may have more if genomes were copied/renamed during processing)
        for genome in input_genomes:
            assert any(genome in listed for listed in listed_genomes), (
                f"Input genome {genome} not found in genome_list.\n"
                f"Input genomes: {sorted(input_genomes)}\n"
                f"Listed genomes: {sorted(listed_genomes)}"
            )
    
    def test_orthofinder_completed(self, orthophyl_run):
        """OrthoFinder ran successfully and produced orthogroups."""
        # OrthoFinder results can be in multiple locations
        of_candidates = [
            orthophyl_run / "annots_prots.fixed" / "OrthoFinder" / "Results_ortho",
            orthophyl_run / "phylo_current" / "Results_ortho",
        ]
        
        of_found = False
        for of_dir in of_candidates:
            if of_dir.exists():
                # Check for key OrthoFinder outputs
                orthogroups = of_dir / "Orthogroups" / "Orthogroups.tsv"
                if orthogroups.exists():
                    of_found = True
                    break
        
        assert of_found, (
            f"OrthoFinder results not found in any expected location:\n"
            f"{[str(d) for d in of_candidates]}"
        )
    
    def test_alignment_quality_files_present(self, orthophyl_run):
        """Alignment quality assessment files are generated."""
        phylo_current = orthophyl_run / "phylo_current"
        
        # Check for SCO (single-copy orthologs) lists
        sco_strict = phylo_current / "SCO_strict"
        sco_3 = phylo_current / "SCO_3"
        
        assert sco_strict.exists(), f"SCO_strict list not found: {sco_strict}"
        assert sco_3.exists(), f"SCO_3 list not found: {sco_3}"
        
        # Check that they're non-empty
        with open(sco_strict) as f:
            strict_ogs = [line.strip() for line in f if line.strip()]
        with open(sco_3) as f:
            sco3_ogs = [line.strip() for line in f if line.strip()]
        
        assert len(strict_ogs) > 0, "SCO_strict list is empty"
        assert len(sco3_ogs) > 0, "SCO_3 list is empty"
        assert len(sco3_ogs) >= len(strict_ogs), (
            "SCO_3 should have at least as many orthogroups as SCO_strict"
        )
