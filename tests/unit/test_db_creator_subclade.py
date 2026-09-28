"""Tests for subclade support in OP_database_tool.py.

Covers subclade_meta config fields + file copies for a built subclade, the
is_backbone flag for a megatree's sparse overview DB, and the lazy
--register-only path (built=false placeholder, no orthophyl_run symlink,
genome_list.txt sourced from members_file) used for on-demand subclade
builds.
"""

import json
from pathlib import Path

import pytest

from conftest import ESCHERICHIA_TAX


@pytest.fixture
def sketch_and_members(tmp_path):
    """A source sketch + members file to be copied into the DB dir."""
    src = tmp_path / "partition_src"
    src.mkdir()
    sketch = src / "Escherichia_1.msh"
    sketch.write_bytes(b"MASH-placeholder")
    members = src / "Escherichia_1.members.txt"
    members.write_text("GCF_001.fna\nGCF_002.fna\nGCF_003.fna\nGCF_004.fna\n")
    return sketch, members


class TestBuiltSubclade:
    def test_writes_subclade_fields_and_copies(
            self, db_creator_module, orthophyl_run_skeleton, sketch_and_members, tmp_path):
        run = orthophyl_run_skeleton(n_genomes=4)
        sketch, members = sketch_and_members
        out = tmp_path / "databases"
        out.mkdir()

        subclade_meta = {
            "is_subclade": True,
            "parent_taxon": "Escherichia",
            "subclade_id": 1,
            "sketch_file": str(sketch),
            "members_file": str(members),
            "source_genome_dir": str(tmp_path / "raw"),
            "built": True,
        }
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=run, clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_1", output_dir=out,
            subclade_meta=subclade_meta,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["is_subclade"] is True
        assert config["parent_taxon"] == "Escherichia"
        assert config["subclade_id"] == 1
        assert config["built"] is True
        # Files copied into the DB dir.
        assert (db_dir / "subclade_sketch.msh").exists()
        assert (db_dir / "subclade_members.txt").exists()
        # In-DB copies referenced by config.
        assert config["sketch_file"] == str(db_dir / "subclade_sketch.msh")
        assert config["members_file"] == str(db_dir / "subclade_members.txt")

    def test_non_subclade_unchanged(
            self, db_creator_module, orthophyl_run_skeleton, tmp_path):
        run = orthophyl_run_skeleton(n_genomes=3)
        out = tmp_path / "databases"
        out.mkdir()
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=run, clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia", output_dir=out,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config.get("is_subclade") is False
        assert config.get("is_backbone") is False
        assert config.get("built") is True
        assert not (db_dir / "subclade_sketch.msh").exists()


class TestBackbone:
    def test_writes_is_backbone_field(
            self, db_creator_module, orthophyl_run_skeleton, tmp_path):
        run = orthophyl_run_skeleton(n_genomes=5)
        out = tmp_path / "databases"
        out.mkdir()
        subclade_meta = {
            "is_backbone": True,
            "parent_taxon": "Escherichia",
        }
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=run, clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia", output_dir=out,
            subclade_meta=subclade_meta,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["is_backbone"] is True
        assert config["is_subclade"] is False
        assert config["parent_taxon"] == "Escherichia"
        assert config["built"] is True

    def test_backbone_gets_real_tree_and_symlink(
            self, db_creator_module, orthophyl_run_skeleton, tmp_path):
        """A backbone DB is BUILT (unlike a lazily-registered subclade) -- it
        gets the real tree copy and orthophyl_run symlink like any other
        built DB."""
        run = orthophyl_run_skeleton(n_genomes=5)
        out = tmp_path / "databases"
        out.mkdir()
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=run, clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia", output_dir=out,
            subclade_meta={"is_backbone": True, "parent_taxon": "Escherichia"},
        )
        assert (db_dir / "orthophyl_run").exists()
        tree_text = (db_dir / "phylogeny.nwk").read_text()
        assert "Placeholder" not in tree_text


class TestRegisterOnly:
    """Lazy registration: built=false placeholder, no tree/HMMs, no
    orthophyl_run symlink -- used so --megatree can partition a huge taxon
    once and build subclades incrementally as queries route to them."""

    def test_requires_is_subclade(self, db_creator_module, tmp_path):
        """register_only without subclade_meta['is_subclade'] is a caller
        bug at the CLI layer (--register-only requires --is-subclade); the
        function itself doesn't enforce it, but a register-only call with no
        subclade_meta at all should still produce a built=false DB."""
        out = tmp_path / "databases"
        out.mkdir()
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=Path("/nonexistent/unused"),
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_2", output_dir=out,
            register_only=True,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["built"] is False

    def test_register_only_writes_placeholder_no_symlink(
            self, db_creator_module, sketch_and_members, tmp_path):
        sketch, members = sketch_and_members
        out = tmp_path / "databases"
        out.mkdir()
        subclade_meta = {
            "is_subclade": True,
            "parent_taxon": "Escherichia",
            "subclade_id": 1,
            "sketch_file": str(sketch),
            "members_file": str(members),
            "source_genome_dir": str(tmp_path / "raw"),
            "n_genomes": 4,
        }
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=Path("/nonexistent/unused"),
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_1", output_dir=out,
            subclade_meta=subclade_meta, register_only=True,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["built"] is False
        assert config["has_trees"] is False
        assert config["is_subclade"] is True
        assert config["n_genomes"] == 4
        # No real OrthoPhyl run to symlink -- registration writes a
        # placeholder tree and skips the symlink entirely.
        assert not (db_dir / "orthophyl_run").exists()
        tree_text = (db_dir / "phylogeny.nwk").read_text()
        assert "Placeholder" in tree_text or "placeholder" in tree_text.lower()

    def test_register_only_genome_list_from_members_file(
            self, db_creator_module, sketch_and_members, tmp_path):
        sketch, members = sketch_and_members
        out = tmp_path / "databases"
        out.mkdir()
        subclade_meta = {
            "is_subclade": True,
            "parent_taxon": "Escherichia",
            "subclade_id": 1,
            "sketch_file": str(sketch),
            "members_file": str(members),
            "source_genome_dir": str(tmp_path / "raw"),
            "n_genomes": 4,
        }
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=Path("/nonexistent/unused"),
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_1", output_dir=out,
            subclade_meta=subclade_meta, register_only=True,
        )
        genome_list = (db_dir / "genome_list.txt").read_text()
        for genome in ["GCF_001.fna", "GCF_002.fna", "GCF_003.fna", "GCF_004.fna"]:
            assert genome in genome_list

    def test_register_only_records_source_genome_dir(
            self, db_creator_module, sketch_and_members, tmp_path):
        sketch, members = sketch_and_members
        out = tmp_path / "databases"
        out.mkdir()
        raw_dir = tmp_path / "raw"
        subclade_meta = {
            "is_subclade": True,
            "parent_taxon": "Escherichia",
            "subclade_id": 1,
            "sketch_file": str(sketch),
            "members_file": str(members),
            "source_genome_dir": str(raw_dir),
            "n_genomes": 4,
        }
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=Path("/nonexistent/unused"),
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_1", output_dir=out,
            subclade_meta=subclade_meta, register_only=True,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["source_genome_dir"] == str(raw_dir)

    def test_register_then_force_build_promotes_to_built_true(
            self, db_creator_module, orthophyl_run_skeleton, sketch_and_members, tmp_path):
        """Checkpoint-collision regression (hazard 1): registering, then
        later force-building the SAME clade name must promote built=false ->
        built=true, not silently keep the placeholder."""
        sketch, members = sketch_and_members
        out = tmp_path / "databases"
        out.mkdir()
        subclade_meta = {
            "is_subclade": True,
            "parent_taxon": "Escherichia",
            "subclade_id": 1,
            "sketch_file": str(sketch),
            "members_file": str(members),
            "source_genome_dir": str(tmp_path / "raw"),
            "n_genomes": 4,
        }
        db_creator_module.create_database_for_run(
            orthophyl_dir=Path("/nonexistent/unused"),
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_1", output_dir=out,
            subclade_meta=subclade_meta, register_only=True,
        )

        run = orthophyl_run_skeleton(n_genomes=4)
        db_dir = db_creator_module.create_database_for_run(
            orthophyl_dir=run, clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_1", output_dir=out,
            subclade_meta={**subclade_meta, "built": True},
            force=True, register_only=False,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["built"] is True
        assert config["has_trees"] is True
        assert (db_dir / "orthophyl_run").exists()

    def test_register_only_without_force_refuses_existing(
            self, db_creator_module, sketch_and_members, tmp_path):
        sketch, members = sketch_and_members
        out = tmp_path / "databases"
        out.mkdir()
        subclade_meta = {
            "is_subclade": True,
            "parent_taxon": "Escherichia",
            "subclade_id": 1,
            "sketch_file": str(sketch),
            "members_file": str(members),
            "source_genome_dir": str(tmp_path / "raw"),
            "n_genomes": 4,
        }
        db_creator_module.create_database_for_run(
            orthophyl_dir=Path("/nonexistent/unused"),
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_name="Escherichia_1", output_dir=out,
            subclade_meta=subclade_meta, register_only=True,
        )
        with pytest.raises(FileExistsError):
            db_creator_module.create_database_for_run(
                orthophyl_dir=Path("/nonexistent/unused"),
                clade_taxonomy=ESCHERICHIA_TAX,
                clade_name="Escherichia_1", output_dir=out,
                subclade_meta=subclade_meta, register_only=True,
            )
