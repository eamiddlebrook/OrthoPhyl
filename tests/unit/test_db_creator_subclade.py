"""Tests for subclade support in create_hierarchical_database.py.

Covers subclade_meta config fields + file copies for a built subclade. (The
lazy --register-only path was removed; every subclade DB is now built.)
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
        assert config.get("built") is True
        assert not (db_dir / "subclade_sketch.msh").exists()
