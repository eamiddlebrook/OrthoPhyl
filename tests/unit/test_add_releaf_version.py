"""Tests for add_releaf_version.py's relative-symlink fix.

Builds a minimal base-version DB + ReLeaf output, runs
create_releaf_version, and checks that every composite-structure symlink
is relative and survives moving the whole database tree -- simulating a
different Singularity bind-mount root.
"""

import json
import os
import shutil
from pathlib import Path

import pytest


def _make_base_version(db_dir, n_genomes=3):
    """A db_dir with a 'current' version pointing at a built OrthoPhyl run."""
    version_dir = db_dir / "v1_orthophyl_initial"
    run = version_dir / "orthophyl_run"

    hmm_dir = run / "OG_alignmentsToHMM" / "hmms_final"
    hmm_dir.mkdir(parents=True)
    (hmm_dir / "OG0000001.hmm").write_text("HMMER3/f placeholder\n")

    trim_dir = run / "phylo_current" / "trimmed_columns"
    trim_dir.mkdir(parents=True)
    (trim_dir / "OG0000001.trm").write_text(">a\nMAA\n")

    iqtree_dir = run / "phylo_current" / "SpeciesTree" / "iqtree"
    iqtree_dir.mkdir(parents=True)
    (iqtree_dir / "report.txt").write_text("iqtree report\n")

    (run / "phylo_current" / "SCO_strict.txt").write_text("OG0000001\n")

    version_dir.mkdir(exist_ok=True)
    (version_dir / "version_info.json").write_text(
        json.dumps({"total_genomes": n_genomes})
    )
    (version_dir / "database_config.json").write_text(
        json.dumps({"n_genomes": n_genomes})
    )
    (version_dir / "genome_list.txt").write_text(
        "\n".join(f"genome_{i}" for i in range(n_genomes)) + "\n"
    )

    current_link = db_dir / "current"
    current_link.symlink_to(version_dir.name)
    return version_dir, run


def _make_releaf_output(releaf_dir, n_new=1):
    prot_aln = releaf_dir / "new_prot_alignments.trm.nm"
    prot_aln.mkdir(parents=True)
    cds_aln = releaf_dir / "new_CDS_alignments.trm.nm"
    cds_aln.mkdir(parents=True)

    tree_dir = releaf_dir / "new_trees"
    tree_dir.mkdir(parents=True)
    (tree_dir / "SCO_strict.CDS.iqtree.treefile.addasm").write_text("(a,b,c);\n")

    annots_dir = releaf_dir / "annots_prots"
    annots_dir.mkdir(parents=True)
    for i in range(n_new):
        (annots_dir / f"new_genome_{i}.faa").write_text(">s\nM\n")

    return releaf_dir


@pytest.fixture
def base_and_releaf(tmp_path):
    db_dir = tmp_path / "Escherichia_db"
    version_dir, run = _make_base_version(db_dir)
    releaf_dir = _make_releaf_output(tmp_path / "releaf_out")
    return db_dir, releaf_dir


class TestCompositeLinksAreRelative:
    def test_hmm_and_trim_links_relative(
        self, add_releaf_version_module, base_and_releaf
    ):
        db_dir, releaf_dir = base_and_releaf
        new_version_dir = add_releaf_version_module.create_releaf_version(
            db_dir=db_dir, releaf_output=releaf_dir, version_name="v2_releaf_test"
        )
        composite = new_version_dir / "orthophyl_run_composite"

        hmm_link = composite / "OG_alignmentsToHMM" / "hmms_final"
        assert not os.path.isabs(os.readlink(hmm_link))
        assert hmm_link.resolve().is_dir()

        trim_link = composite / "phylo_current" / "trimmed_columns"
        assert not os.path.isabs(os.readlink(trim_link))
        assert trim_link.resolve().is_dir()

    def test_alignment_links_relative(
        self, add_releaf_version_module, base_and_releaf
    ):
        db_dir, releaf_dir = base_and_releaf
        new_version_dir = add_releaf_version_module.create_releaf_version(
            db_dir=db_dir, releaf_output=releaf_dir, version_name="v2_releaf_test"
        )
        composite = new_version_dir / "orthophyl_run_composite"

        prot_link = composite / "phylo_current" / "AlignmentsProts.trm"
        assert not os.path.isabs(os.readlink(prot_link))
        assert prot_link.resolve().is_dir()

        cds_link = composite / "phylo_current" / "AlignmentsCDS.trm"
        assert not os.path.isabs(os.readlink(cds_link))
        assert cds_link.resolve().is_dir()

    def test_tree_and_phylogeny_links_relative(
        self, add_releaf_version_module, base_and_releaf
    ):
        db_dir, releaf_dir = base_and_releaf
        new_version_dir = add_releaf_version_module.create_releaf_version(
            db_dir=db_dir, releaf_output=releaf_dir, version_name="v2_releaf_test"
        )
        composite = new_version_dir / "orthophyl_run_composite"
        tree_link = composite / "FINAL_SPECIES_TREES" / "SCO_strict.CDS.iqtree.treefile"
        assert not os.path.isabs(os.readlink(tree_link))
        assert tree_link.resolve().exists()

        phylogeny_link = new_version_dir / "phylogeny.nwk"
        assert not os.path.isabs(os.readlink(phylogeny_link))
        assert phylogeny_link.resolve().exists()

        releaf_output_link = new_version_dir / "releaf_output"
        assert not os.path.isabs(os.readlink(releaf_output_link))
        assert releaf_output_link.resolve().is_dir()

    def test_orthophyl_run_link_is_already_relative_bare_name(
        self, add_releaf_version_module, base_and_releaf
    ):
        """orthophyl_run -> orthophyl_run_composite was already a bare
        relative name before this fix and should stay untouched."""
        db_dir, releaf_dir = base_and_releaf
        new_version_dir = add_releaf_version_module.create_releaf_version(
            db_dir=db_dir, releaf_output=releaf_dir, version_name="v2_releaf_test"
        )
        link = new_version_dir / "orthophyl_run"
        assert os.readlink(link) == "orthophyl_run_composite"


class TestCompositeLinksSurviveTreeMove:
    def test_links_resolve_after_moving_whole_database(
        self, add_releaf_version_module, tmp_path
    ):
        """Simulates a different Singularity bind root: build the composite
        version under tmp/a/, move the whole tree to tmp/b/, and confirm
        every relative link still resolves."""
        a = tmp_path / "a"
        a.mkdir()
        db_dir = a / "Escherichia_db"
        _make_base_version(db_dir)
        releaf_dir = _make_releaf_output(a / "releaf_out")

        new_version_dir = add_releaf_version_module.create_releaf_version(
            db_dir=db_dir, releaf_output=releaf_dir, version_name="v2_releaf_test"
        )
        rel_new_version_dir = new_version_dir.relative_to(tmp_path)

        b = tmp_path / "b"
        shutil.move(str(a), str(b))

        moved_version_dir = tmp_path / "b" / rel_new_version_dir.relative_to("a")
        composite = moved_version_dir / "orthophyl_run_composite"

        assert (composite / "OG_alignmentsToHMM" / "hmms_final" / "OG0000001.hmm").exists()
        assert (composite / "phylo_current" / "trimmed_columns" / "OG0000001.trm").exists()
        assert (composite / "phylo_current" / "AlignmentsProts.trm").resolve().is_dir()
        assert (composite / "phylo_current" / "AlignmentsCDS.trm").resolve().is_dir()
        assert (moved_version_dir / "phylogeny.nwk").resolve().exists()
