"""Tests for OP_database_tool.py --migrate (plan Part 2/3: host<->container
portability repair for pre-existing, pre-nesting/pre-relative database dirs).

Builds a fake *legacy* database tree by hand (not via the make_db_dir
fixture, which already writes the current relative-path convention) --
absolute /databases/... paths baked into the index and configs, a subclade
that needs to move into the nested layout, an absolute orthophyl_run
symlink, and one dangling link whose target was deleted.
"""

import json
import os
from pathlib import Path

import pytest

ESCHERICHIA_TAX = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia"
)


def _write_config(db_dir, **fields):
    config = {
        "created": "2025-01-01T00:00:00",
        "version": "1.0",
        "database_type": "hierarchical_taxonomy",
        "n_genomes": 5,
        "has_hmms": True,
        "has_trees": True,
        "available_tree_methods": ["iqtree"],
        "available_data_types": ["CDS"],
    }
    config.update(fields)
    (db_dir / "database_config.json").write_text(json.dumps(config, indent=2))
    return config


@pytest.fixture
def legacy_tree(tmp_path):
    """Build a fake legacy (pre-852ddd8, pre-nesting) database tree rooted
    at tmp_path/databases, simulating one that was built inside a container
    bind-mounted at /databases and is now being inspected on the host.

    Returns (output_dir, old_prefix) where old_prefix is the fake container
    mount path ("/databases") baked into the legacy absolute fields.
    """
    output_dir = tmp_path / "databases"
    output_dir.mkdir()
    old_prefix = "/databases"

    # A flat, non-megatree DB with no path fields beyond orthophyl_source --
    # should end up untouched except for relativizing orthophyl_source.
    flat_dir = output_dir / "Flat_db"
    flat_dir.mkdir()
    run_src = tmp_path / "run_src_flat"
    run_src.mkdir()
    (flat_dir / "orthophyl_run").symlink_to(run_src, target_is_directory=True)
    _write_config(
        flat_dir,
        clade_name="Flat",
        clade_taxonomy="d__Bacteria;g__Flat",
        clade_rank="g", clade_rank_name="genus",
        orthophyl_source=f"{old_prefix}/Flat_db/orthophyl_run",
    )

    # A subclade DB sitting FLAT (pre-nesting bug) with absolute, mount-
    # prefixed sketch_file/members_file/source_genome_dir, that needs to:
    #  (a) move under output_dir/Escherichia/Escherichia_1_db
    #  (b) have those fields remapped off /databases and relativized
    sub_dir = output_dir / "Escherichia_1_db"
    sub_dir.mkdir()
    sketch_file = sub_dir / "subclade_sketch.msh"
    sketch_file.write_bytes(b"MASH-placeholder")
    members_file = sub_dir / "subclade_members.txt"
    members_file.write_text("g0.fna\ng1.fna\n")
    source_genome_dir = tmp_path / "source_genomes" / "Escherichia_1"
    source_genome_dir.mkdir(parents=True)
    sub_run_src = tmp_path / "run_src_sub"
    sub_run_src.mkdir()
    (sub_dir / "orthophyl_run").symlink_to(sub_run_src, target_is_directory=True)
    _write_config(
        sub_dir,
        clade_name="Escherichia_1",
        clade_taxonomy=ESCHERICHIA_TAX,
        clade_rank="g", clade_rank_name="genus",
        is_subclade=True, is_backbone=False,
        parent_taxon="Escherichia", subclade_id=1, built=True,
        sketch_file=f"{old_prefix}/Escherichia_1_db/subclade_sketch.msh",
        members_file=f"{old_prefix}/Escherichia_1_db/subclade_members.txt",
        source_genome_dir=str(source_genome_dir),
        orthophyl_source=f"{old_prefix}/Escherichia_1_db/orthophyl_run",
    )

    # A DB with a DANGLING orthophyl_run target (its .pipeline_runs dir was
    # deleted) -- must be reported, never deleted.
    dangling_dir = output_dir / "Dangling_db"
    dangling_dir.mkdir()
    deleted_target = tmp_path / "deleted_run_src"  # never created
    (dangling_dir / "orthophyl_run").symlink_to(deleted_target, target_is_directory=True)
    _write_config(
        dangling_dir,
        clade_name="Dangling",
        clade_taxonomy="d__Bacteria;g__Dangling",
        clade_rank="g", clade_rank_name="genus",
        orthophyl_source=f"{old_prefix}/Dangling_db/orthophyl_run",
    )

    # Legacy database_index.json with absolute database_dir entries, as
    # written before commit 852ddd8.
    index = {
        "databases": [
            {"clade_name": "Flat", "database_dir": f"{old_prefix}/Flat_db"},
            {"clade_name": "Escherichia_1", "database_dir": f"{old_prefix}/Escherichia_1_db"},
            {"clade_name": "Dangling", "database_dir": f"{old_prefix}/Dangling_db"},
        ]
    }
    (output_dir / "database_index.json").write_text(json.dumps(index, indent=2))

    # Legacy orthophyl_runs.tsv with an absolute orthophyl_dir column.
    tsv_lines = [
        "clade_name\torthophyl_dir",
        f"Flat\t{old_prefix}/orthophyl_run_dir_flat",
    ]
    (output_dir / "orthophyl_runs.tsv").write_text("\n".join(tsv_lines) + "\n")

    return output_dir, old_prefix


class TestInferStripPrefix:
    def test_infers_from_index(self, db_creator_module, legacy_tree):
        output_dir, old_prefix = legacy_tree
        assert db_creator_module._infer_strip_prefix(output_dir) == old_prefix

    def test_none_when_no_index(self, db_creator_module, tmp_path):
        assert db_creator_module._infer_strip_prefix(tmp_path / "nope") is None


class TestMigrateDryRun:
    def test_dry_run_changes_nothing_on_disk(self, db_creator_module, legacy_tree):
        output_dir, old_prefix = legacy_tree
        before = sorted(str(p) for p in output_dir.rglob("*"))

        report = db_creator_module.migrate_databases(output_dir, dry_run=True)

        after = sorted(str(p) for p in output_dir.rglob("*"))
        assert before == after
        # The subclade move and field rewrites are still *reported*.
        assert report["dry_run"] is True
        assert any(m["clade_name"] == "Escherichia_1" for m in report["moved"])
        assert any(r["clade_name"] == "Escherichia_1" for r in report["rewritten"])
        assert any(d["clade_name"] == "Dangling" for d in report["dangling_links"])
        # database_config.json for the subclade is untouched (still absolute).
        sub_config = json.loads(
            (output_dir / "Escherichia_1_db" / "database_config.json").read_text())
        assert sub_config["sketch_file"].startswith(old_prefix)


class TestMigrateLive:
    def test_subclade_moved_to_nested_layout(self, db_creator_module, legacy_tree):
        output_dir, _ = legacy_tree
        db_creator_module.migrate_databases(output_dir, dry_run=False)

        assert not (output_dir / "Escherichia_1_db").exists()
        new_dir = output_dir / "Escherichia" / "Escherichia_1_db"
        assert new_dir.is_dir()
        config = json.loads((new_dir / "database_config.json").read_text())
        assert config["clade_name"] == "Escherichia_1"

    def test_fields_rewritten_relative(self, db_creator_module, legacy_tree):
        output_dir, _ = legacy_tree
        db_creator_module.migrate_databases(output_dir, dry_run=False)

        new_dir = output_dir / "Escherichia" / "Escherichia_1_db"
        config = json.loads((new_dir / "database_config.json").read_text())
        assert not Path(config["sketch_file"]).is_absolute()
        assert not Path(config["members_file"]).is_absolute()
        assert (new_dir / config["sketch_file"]).exists()
        assert (new_dir / config["members_file"]).exists()
        # source_genome_dir lived outside db_dir -- relativized but not moved.
        assert not Path(config["source_genome_dir"]).is_absolute()
        assert (new_dir / config["source_genome_dir"]).exists()

    def test_flat_db_orthophyl_source_relativized(self, db_creator_module, legacy_tree):
        output_dir, _ = legacy_tree
        db_creator_module.migrate_databases(output_dir, dry_run=False)

        flat_dir = output_dir / "Flat_db"
        config = json.loads((flat_dir / "database_config.json").read_text())
        assert not Path(config["orthophyl_source"]).is_absolute()
        link = flat_dir / "orthophyl_run"
        assert link.is_symlink()
        assert not os.path.isabs(os.readlink(link))
        assert link.resolve().is_dir()

    def test_dangling_link_reported_and_left_in_place(self, db_creator_module, legacy_tree):
        output_dir, _ = legacy_tree
        report = db_creator_module.migrate_databases(output_dir, dry_run=False)

        assert any(d["clade_name"] == "Dangling" for d in report["dangling_links"])
        dangling_dir = output_dir / "Dangling_db"
        link = dangling_dir / "orthophyl_run"
        assert link.is_symlink()
        assert not link.resolve().exists()

    def test_router_loads_all_databases_after_migration(
            self, db_creator_module, router_module, legacy_tree, tmp_path):
        output_dir, _ = legacy_tree
        db_creator_module.migrate_databases(output_dir, dry_run=False)

        router = router_module.MultiDatabaseRouter(
            database_dir=output_dir, output_dir=tmp_path / "routing_out")
        names = {db["clade_name"] for db in router.databases}
        assert names == {"Flat", "Escherichia_1", "Dangling"}

    def test_runs_tsv_relativized(self, db_creator_module, legacy_tree):
        output_dir, _ = legacy_tree
        db_creator_module.migrate_databases(output_dir, dry_run=False)

        tsv = (output_dir / "orthophyl_runs.tsv").read_text()
        data_line = tsv.splitlines()[1]
        orthophyl_dir = data_line.split("\t")[1]
        assert not Path(orthophyl_dir).is_absolute()

    def test_second_run_is_noop(self, db_creator_module, legacy_tree):
        output_dir, _ = legacy_tree
        db_creator_module.migrate_databases(output_dir, dry_run=False)

        second = db_creator_module.migrate_databases(output_dir, dry_run=False)
        assert second["moved"] == []
        for r in second["rewritten"]:
            assert r["changes"] == {}
