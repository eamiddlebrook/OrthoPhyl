"""Tests for the multi-database assembly router (plan section 2)."""

import json
from pathlib import Path

import pytest

from conftest import GAIELLALES_TAX, ESCHERICHIA_TAX


@pytest.fixture
def Router(router_module):
    return router_module.MultiDatabaseRouter


@pytest.fixture
def two_db_dir(make_db_dir, tmp_path):
    """A database dir with a family-level and a genus-level DB sharing a lineage."""
    parent = tmp_path / "databases"
    make_db_dir(
        clade_name="Enterobacteriaceae",
        clade_taxonomy=(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae"
        ),
        clade_rank="f",
        clade_rank_name="family",
        n_genomes=200,
        parent=parent,
    )
    make_db_dir(
        clade_name="Escherichia",
        clade_taxonomy=ESCHERICHIA_TAX,
        clade_rank="g",
        clade_rank_name="genus",
        n_genomes=2840,
        tree_methods=["iqtree", "fasttree"],
        data_types=["CDS", "PROT"],
        parent=parent,
    )
    return parent


class TestLoadDatabases:
    def test_load_by_glob_when_no_index(self, Router, two_db_dir, tmp_path):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        names = sorted(db["clade_name"] for db in router.databases)
        assert names == ["Enterobacteriaceae", "Escherichia"]

    def test_empty_dir_raises(self, Router, tmp_path):
        empty = tmp_path / "empty_db_dir"
        empty.mkdir()
        with pytest.raises(ValueError, match="No databases found"):
            Router(database_dir=empty, output_dir=tmp_path / "out")


class TestFindMatching:
    def test_most_specific_first(self, Router, two_db_dir, tmp_path):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        matches = router.find_matching_databases(
            ESCHERICHIA_TAX + ";s__Escherichia_coli"
        )
        assert matches, "expected at least one match"
        # Genus DB (more specific) must sort ahead of the family DB.
        assert matches[0]["database"]["clade_name"] == "Escherichia"

    def test_no_match_for_novel_phylum(self, Router, two_db_dir, tmp_path):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        assert router.find_matching_databases(GAIELLALES_TAX) == []


class TestRouteToReLeaf:
    def _assembly(self, tmp_path):
        asm = tmp_path / "query.fna"
        asm.write_text(">c\nACGT\n")
        return asm

    def test_routes_to_releaf_on_match(self, Router, two_db_dir, tmp_path):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        decision = router.route_assembly(
            self._assembly(tmp_path), ESCHERICHIA_TAX + ";s__", "q1"
        )
        assert decision["pipeline"] == "ReLeaf"
        assert decision["matched_database"] == "Escherichia"
        assert decision["matched_rank"] == "genus"

    def test_prefers_iqtree_cds(self, Router, two_db_dir, tmp_path):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        decision = router.route_assembly(
            self._assembly(tmp_path), ESCHERICHIA_TAX + ";s__", "q1"
        )
        assert decision["tree_method"] == "iqtree"
        assert decision["tree_data"] == "CDS"

    def test_writes_decision_and_summary(self, Router, two_db_dir, tmp_path):
        out = tmp_path / "out"
        router = Router(database_dir=two_db_dir, output_dir=out)
        router.route_assembly(self._assembly(tmp_path), ESCHERICHIA_TAX + ";s__", "q1")
        decision_file = out / "routing_decision_q1.json"
        summary_file = out / "routing_summary_q1.txt"
        assert decision_file.exists() and summary_file.exists()
        # Decision JSON is valid and round-trips.
        json.loads(decision_file.read_text())

    def test_copies_assembly_into_releaf_input(self, Router, two_db_dir, tmp_path):
        out = tmp_path / "out"
        router = Router(database_dir=two_db_dir, output_dir=out)
        router.route_assembly(self._assembly(tmp_path), ESCHERICHIA_TAX + ";s__", "q1")
        assert (out / "releaf_input" / "q1.fna").exists()


class TestRouteToOrthoPhyl:
    def _assembly(self, tmp_path):
        asm = tmp_path / "novel.fna"
        asm.write_text(">c\nACGT\n")
        return asm

    def test_novel_routes_to_orthophyl(self, Router, two_db_dir, tmp_path):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        decision = router.route_assembly(
            self._assembly(tmp_path),
            GAIELLALES_TAX + ";f__Gaiellaceae;g__Gaiella;s__",
            "novel1",
        )
        assert decision["pipeline"] == "OrthoPhyl"

    def test_download_rank_uses_genus_when_species_empty(
        self, Router, two_db_dir, tmp_path
    ):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        decision = router.route_assembly(
            self._assembly(tmp_path),
            GAIELLALES_TAX + ";f__Gaiellaceae;g__Gaiella;s__",
            "novel1",
        )
        assert decision["download_rank"] == "genus"
        assert decision["download_value"] == "Gaiella"

    def test_download_falls_back_to_family(self, Router, two_db_dir, tmp_path):
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        decision = router.route_assembly(
            self._assembly(tmp_path),
            GAIELLALES_TAX + ";f__Gaiellaceae;g__;s__",
            "novel1",
        )
        assert decision["download_rank"] == "family"
        assert decision["download_value"] == "Gaiellaceae"


class TestBatchRoute:
    def test_batch_mixed_counts(self, Router, two_db_dir, tmp_path):
        g1 = tmp_path / "g1.fna"
        g1.write_text(">c\nACGT\n")
        g2 = tmp_path / "g2.fna"
        g2.write_text(">c\nACGT\n")
        batch = tmp_path / "batch.tsv"
        batch.write_text(
            f"{g1}\t{ESCHERICHIA_TAX};s__\tmatched1\n"
            f"{g2}\t{GAIELLALES_TAX};f__Gaiellaceae;g__Gaiella;s__\tnovel1\n"
        )
        out = tmp_path / "out"
        router = Router(database_dir=two_db_dir, output_dir=out)
        decisions = router.batch_route(batch)
        pipelines = sorted(d["pipeline"] for d in decisions)
        assert pipelines == ["OrthoPhyl", "ReLeaf"]
        assert (out / "batch_routing_summary.txt").exists()

    def test_batch_strips_quotes_and_expands_tilde(
        self, Router, two_db_dir, tmp_path, monkeypatch
    ):
        # Place the genome under a fake HOME so ~ expansion resolves to it.
        home = tmp_path / "home"
        home.mkdir()
        monkeypatch.setenv("HOME", str(home))
        g = home / "g.fna"
        g.write_text(">c\nACGT\n")
        batch = tmp_path / "batch.tsv"
        batch.write_text(f'~/g.fna\t"{ESCHERICHIA_TAX};s__"\tq\n')
        router = Router(database_dir=two_db_dir, output_dir=tmp_path / "out")
        decisions = router.batch_route(batch)
        assert len(decisions) == 1
        assert decisions[0]["pipeline"] == "ReLeaf"
