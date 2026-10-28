"""Tests for OP_database_tool.py (plan section 3)."""

import json
import os
import shutil
from pathlib import Path

import pytest

from conftest import ESCHERICHIA_TAX


class TestValidateOrthophylRun:
    def test_valid_run(self, db_creator_module, orthophyl_run_skeleton):
        run = orthophyl_run_skeleton(n_genomes=3)
        v = db_creator_module.validate_orthophyl_run(run)
        assert v["valid"] is True
        assert v["has_hmms"] is True
        assert v["has_trees"] is True
        assert v["n_genomes"] == 3

    def test_prefers_iqtree_tree(self, db_creator_module, orthophyl_run_skeleton):
        run = orthophyl_run_skeleton()
        v = db_creator_module.validate_orthophyl_run(run)
        assert "iqtree" in v["tree_file"].name.lower()

    def test_missing_hmms_warns(self, db_creator_module, orthophyl_run_skeleton):
        run = orthophyl_run_skeleton(with_hmms=False)
        v = db_creator_module.validate_orthophyl_run(run)
        assert v["has_hmms"] is False
        # Alignments still present, and genomes counted -> still valid.
        assert v["has_alignments"] is True
        assert v["valid"] is True

    def test_genome_count_fallback_to_faa(self, db_creator_module, tmp_path):
        run = tmp_path / "run"
        hmm = run / "OG_alignmentsToHMM" / "hmms_final"
        hmm.mkdir(parents=True)
        (hmm / "OG.hmm").write_text("x\n")
        prots = run / "annots_prots"
        prots.mkdir(parents=True)
        for i in range(4):
            (prots / f"g{i}.faa").write_text(">s\nM\n")
        v = db_creator_module.validate_orthophyl_run(run)
        assert v["n_genomes"] == 4

    def test_invalid_when_no_genomes(self, db_creator_module, tmp_path):
        run = tmp_path / "run"
        hmm = run / "OG_alignmentsToHMM" / "hmms_final"
        hmm.mkdir(parents=True)
        (hmm / "OG.hmm").write_text("x\n")
        v = db_creator_module.validate_orthophyl_run(run)
        assert v["n_genomes"] == 0
        assert v["valid"] is False

    def test_missing_dir_raises(self, db_creator_module, tmp_path):
        with pytest.raises(FileNotFoundError):
            db_creator_module.validate_orthophyl_run(tmp_path / "nope")


class TestDatabaseExists:
    def test_detects_existing(self, db_creator_module, make_db_dir, tmp_path):
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Escherichia",
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_rank="g",
            clade_rank_name="genus",
            parent=parent,
        )
        assert db_creator_module.database_exists("Escherichia", parent) is not None

    def test_safe_name_sanitization(self, db_creator_module, tmp_path):
        # "Foo Bar/Baz" -> "Foo_Bar_Baz_db"
        parent = tmp_path / "databases"
        parent.mkdir()
        db_dir = parent / "Foo_Bar_Baz_db"
        db_dir.mkdir()
        (db_dir / "database_config.json").write_text("{}")
        assert db_creator_module.database_exists("Foo Bar/Baz", parent) == db_dir

    def test_absent_returns_none(self, db_creator_module, tmp_path):
        parent = tmp_path / "databases"
        parent.mkdir()
        assert db_creator_module.database_exists("Nothing", parent) is None


class TestCreateDatabaseForRun:
    def test_creates_all_artifacts(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton(n_genomes=3)
        out = tmp_path / "databases"
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out
        )
        assert (db_dir / "database_config.json").exists()
        assert (db_dir / "phylogeny.nwk").exists()
        assert (db_dir / "genome_list.txt").exists()
        assert (db_dir / "README.txt").exists()
        assert (db_dir / "orthophyl_run").is_symlink()

    def test_orthophyl_run_symlink_is_relative(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton(n_genomes=3)
        out = tmp_path / "databases"
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out
        )
        link = db_dir / "orthophyl_run"
        raw_target = os.readlink(link)
        assert not os.path.isabs(raw_target)
        assert link.resolve() == run.resolve()

    def test_orthophyl_run_symlink_survives_tree_move(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        """Simulates a different Singularity bind root: build the DB and run
        under one path, move the whole tree elsewhere, and confirm the
        symlink still resolves -- a relative link survives, an absolute one
        (baked to the old mount) would not."""
        a = tmp_path / "a"
        a.mkdir()
        run = orthophyl_run_skeleton(n_genomes=3)
        shutil.move(str(run), str(a / "orthophyl_run"))
        run_in_a = a / "orthophyl_run"
        out = a / "databases"
        db_creator_module.create_database_for_run(
            run_in_a, ESCHERICHIA_TAX, "Escherichia", out
        )

        b = tmp_path / "b"
        shutil.move(str(a), str(b))

        db_dir = b / "databases" / "Escherichia_db"
        link = db_dir / "orthophyl_run"
        assert link.is_symlink()
        assert link.resolve() == (b / "orthophyl_run").resolve()
        assert (link / "genome_list").exists()

    def test_raises_on_existing_without_force(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        db_creator_module.create_database_for_run(run, ESCHERICHIA_TAX, "Escherichia", out)
        with pytest.raises(FileExistsError):
            db_creator_module.create_database_for_run(
                run, ESCHERICHIA_TAX, "Escherichia", out
            )

    def test_force_rebuilds(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        db_creator_module.create_database_for_run(run, ESCHERICHIA_TAX, "Escherichia", out)
        # Should not raise with force=True.
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out, force=True
        )
        assert (db_dir / "database_config.json").exists()

    def test_config_records_rank(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["clade_rank"] == "g"
        assert config["clade_rank_name"] == "genus"

    def test_taxonomy_source_defaults_to_ncbi(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        """Every existing caller omits taxonomy_source/qc_applied; the config
        must still default to 'ncbi'/True so they stay byte-identical."""
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["taxonomy_source"] == "ncbi"
        assert config["qc_applied"] is True

    def test_taxonomy_source_user_supplied_and_qc_not_applied(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out,
            taxonomy_source="user_supplied", qc_applied=False,
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["taxonomy_source"] == "user_supplied"
        assert config["qc_applied"] is False


class TestMetadataFile:
    """--metadata-file (plan section 4): copies an assembly metadata TSV
    (e.g. all_asm_acc_metadata) into the DB for later --geo figure use,
    recorded as a path relative to db_dir per the Singularity relative-path
    convention (see assembly_router.py's sketch/members_file handling)."""

    def test_metadata_file_copied_and_recorded_relative(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        metadata = tmp_path / "all_asm_acc_metadata"
        metadata.write_text("GCF_000001 species1 strain1 bios1 1234 ok ok USA\n")

        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out, metadata_file=str(metadata)
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["metadata_file"] == "assembly_metadata.tsv"
        copied = db_dir / config["metadata_file"]
        assert copied.exists()
        assert copied.read_text() == metadata.read_text()

    def test_metadata_file_resolves_after_move(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        """The recorded path is relative to db_dir, so it must still resolve
        after the whole tree moves -- same guarantee as sketch_file/
        members_file (test_orthophyl_run_symlink_survives_tree_move)."""
        a = tmp_path / "a"
        a.mkdir()
        run = orthophyl_run_skeleton()
        shutil.move(str(run), str(a / "orthophyl_run"))
        run_in_a = a / "orthophyl_run"
        metadata = a / "all_asm_acc_metadata"
        metadata.write_text("GCF_000001 species1 strain1 bios1 1234 ok ok USA\n")
        out = a / "databases"
        db_creator_module.create_database_for_run(
            run_in_a, ESCHERICHIA_TAX, "Escherichia", out, metadata_file=str(metadata)
        )

        b = tmp_path / "b"
        shutil.move(str(a), str(b))

        db_dir = b / "databases" / "Escherichia_db"
        config = json.loads((db_dir / "database_config.json").read_text())
        resolved = db_dir / config["metadata_file"]
        assert resolved.exists()

    def test_no_metadata_file_arg_leaves_field_none(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        """Every pre-existing caller omits --metadata-file; the config must
        still default to metadata_file=None so they stay byte-identical."""
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["metadata_file"] is None

    def test_nonexistent_metadata_file_leaves_field_none(
        self, db_creator_module, orthophyl_run_skeleton, tmp_path
    ):
        run = orthophyl_run_skeleton()
        out = tmp_path / "databases"
        db_dir = db_creator_module.create_database_for_run(
            run, ESCHERICHIA_TAX, "Escherichia", out,
            metadata_file=str(tmp_path / "does_not_exist.tsv"),
        )
        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["metadata_file"] is None
        assert not (db_dir / "assembly_metadata.tsv").exists()


class TestParseInputTable:
    def test_parses_three_columns(self, db_creator_module, tmp_path):
        tsv = tmp_path / "runs.tsv"
        tsv.write_text(f"Escherichia\t/data/ecoli\t{ESCHERICHIA_TAX}\n")
        runs = db_creator_module.parse_input_table(tsv)
        assert len(runs) == 1
        assert runs[0]["clade_name"] == "Escherichia"

    def test_skips_comments_and_short_rows(self, db_creator_module, tmp_path):
        tsv = tmp_path / "runs.tsv"
        tsv.write_text(
            "# comment\n"
            f"Escherichia\t/data/ecoli\t{ESCHERICHIA_TAX}\n"
            "BadRow\tonly_two_fields\n"
        )
        runs = db_creator_module.parse_input_table(tsv)
        assert len(runs) == 1


class TestMasterIndex:
    def test_index_contents(self, db_creator_module, tmp_path):
        out = tmp_path / "databases"
        out.mkdir()
        databases = [
            {
                "clade_name": "Escherichia",
                "clade_taxonomy": ESCHERICHIA_TAX,
                "database_dir": out / "Escherichia_db",
                "n_genomes": 42,
            }
        ]
        db_creator_module.create_master_index(databases, out)
        index = json.loads((out / "database_index.json").read_text())
        assert index["n_databases"] == 1
        assert index["databases"][0]["clade_rank"] == "g"
        assert (out / "database_summary.txt").exists()

    def test_database_dir_is_relative_to_output_dir(self, db_creator_module, tmp_path):
        out = tmp_path / "databases"
        out.mkdir()
        databases = [
            {
                "clade_name": "Escherichia",
                "clade_taxonomy": ESCHERICHIA_TAX,
                "database_dir": out / "Escherichia_db",
                "n_genomes": 42,
            }
        ]
        db_creator_module.create_master_index(databases, out)
        index = json.loads((out / "database_index.json").read_text())
        entry = index["databases"][0]["database_dir"]
        assert not os.path.isabs(entry)
        assert entry == "Escherichia_db"

    def test_taxonomy_source_surfaced_when_present(self, db_creator_module, tmp_path):
        out = tmp_path / "databases"
        out.mkdir()
        databases = [
            {
                "clade_name": "Blorptaxon",
                "clade_taxonomy": "g__Blorptaxon",
                "database_dir": out / "Blorptaxon_db",
                "n_genomes": 12,
                "taxonomy_source": "user_supplied",
            }
        ]
        db_creator_module.create_master_index(databases, out)
        index = json.loads((out / "database_index.json").read_text())
        assert index["databases"][0]["taxonomy_source"] == "user_supplied"

    def test_taxonomy_source_absent_when_not_provided(self, db_creator_module, tmp_path):
        """Mirrors the is_subclade/built conditional-key pattern: omitting the
        key from the caller's dict must omit it from the index entry too, not
        default it to 'ncbi' (that .get() default belongs to the readers)."""
        out = tmp_path / "databases"
        out.mkdir()
        databases = [
            {
                "clade_name": "Escherichia",
                "clade_taxonomy": ESCHERICHIA_TAX,
                "database_dir": out / "Escherichia_db",
                "n_genomes": 42,
            }
        ]
        db_creator_module.create_master_index(databases, out)
        index = json.loads((out / "database_index.json").read_text())
        assert "taxonomy_source" not in index["databases"][0]


class TestRemoveDatabases:
    """remove_databases(names, output_dir, confirm=False) -- dry-run by
    default, deletes + refreshes the master index only with confirm=True."""

    def test_dry_run_does_not_delete(self, db_creator_module, make_db_dir, tmp_path):
        out = tmp_path / "databases"
        db_dir = make_db_dir(
            clade_name="Alpha", clade_taxonomy="d__Bacteria;g__Alpha",
            clade_rank="g", clade_rank_name="genus", n_genomes=5, parent=out,
        )
        result = db_creator_module.remove_databases(["Alpha"], out, confirm=False)
        assert result["dry_run"] is True
        assert db_dir.exists()
        assert [r["clade_name"] for r in result["removed"]] == ["Alpha"]
        assert result["not_found"] == []
        # Dry run must not touch the index at all.
        assert not (out / "database_index.json").exists()

    def test_confirm_removes_directory(self, db_creator_module, make_db_dir, tmp_path):
        out = tmp_path / "databases"
        db_dir = make_db_dir(
            clade_name="Alpha", clade_taxonomy="d__Bacteria;g__Alpha",
            clade_rank="g", clade_rank_name="genus", n_genomes=5, parent=out,
        )
        result = db_creator_module.remove_databases(["Alpha"], out, confirm=True)
        assert result["dry_run"] is False
        assert not db_dir.exists()
        assert [r["clade_name"] for r in result["removed"]] == ["Alpha"]

    def test_not_found_name_warns_and_continues(
            self, db_creator_module, make_db_dir, tmp_path):
        out = tmp_path / "databases"
        alpha_dir = make_db_dir(
            clade_name="Alpha", clade_taxonomy="d__Bacteria;g__Alpha",
            clade_rank="g", clade_rank_name="genus", n_genomes=5, parent=out,
        )
        result = db_creator_module.remove_databases(
            ["Alpha", "NoSuchClade"], out, confirm=True)
        assert not alpha_dir.exists()
        assert [r["clade_name"] for r in result["removed"]] == ["Alpha"]
        assert result["not_found"] == ["NoSuchClade"]

    def test_refreshes_index_after_removal(
            self, db_creator_module, make_db_dir, tmp_path):
        out = tmp_path / "databases"
        make_db_dir(
            clade_name="Alpha", clade_taxonomy="d__Bacteria;g__Alpha",
            clade_rank="g", clade_rank_name="genus", n_genomes=5, parent=out,
        )
        make_db_dir(
            clade_name="Beta", clade_taxonomy="d__Bacteria;g__Beta",
            clade_rank="g", clade_rank_name="genus", n_genomes=7, parent=out,
        )
        db_creator_module.remove_databases(["Alpha"], out, confirm=True)
        index = json.loads((out / "database_index.json").read_text())
        names = [d["clade_name"] for d in index["databases"]]
        assert names == ["Beta"]

    def test_removing_last_database_writes_empty_index(
            self, db_creator_module, make_db_dir, tmp_path):
        out = tmp_path / "databases"
        make_db_dir(
            clade_name="Alpha", clade_taxonomy="d__Bacteria;g__Alpha",
            clade_rank="g", clade_rank_name="genus", n_genomes=5, parent=out,
        )
        db_creator_module.remove_databases(["Alpha"], out, confirm=True)
        index = json.loads((out / "database_index.json").read_text())
        assert index["n_databases"] == 0
        assert index["databases"] == []

    def test_missing_output_dir_is_all_not_found(self, db_creator_module, tmp_path):
        out = tmp_path / "does_not_exist"
        result = db_creator_module.remove_databases(["Alpha"], out, confirm=True)
        assert result["removed"] == []
        assert result["not_found"] == ["Alpha"]
        assert not out.exists()

    def test_dry_run_reports_context_even_with_unreadable_config(
            self, db_creator_module, make_db_dir, tmp_path):
        out = tmp_path / "databases"
        db_dir = make_db_dir(
            clade_name="Alpha", clade_taxonomy="d__Bacteria;g__Alpha",
            clade_rank="g", clade_rank_name="genus", n_genomes=5, parent=out,
        )
        (db_dir / "database_config.json").write_text("{not valid json")
        result = db_creator_module.remove_databases(["Alpha"], out, confirm=False)
        assert [r["clade_name"] for r in result["removed"]] == ["Alpha"]
        assert result["removed"][0]["clade_taxonomy"] is None
        assert result["removed"][0]["n_genomes"] is None
