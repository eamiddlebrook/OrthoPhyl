"""Tests for subclade routing in assembly_router.py.

Covers the MASH tie-break among same-parent subclades, the built vs unbuilt
decision (ReLeaf vs OrthoPhyl_subclade_build), and that _save_decision /
_generate_batch_summary handle the new pipeline value without KeyError. mash is
mocked at the _run_mash boundary.
"""

import json
from pathlib import Path

import pytest

# A genus that both subclades share (indistinguishable by taxonomy).
GENUS_TAX = (
    "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;"
    "f__Eubacteriaceae;g__Andreesenella"
)


@pytest.fixture
def Router(router_module):
    return router_module.MultiDatabaseRouter


def _assembly(tmp_path):
    asm = tmp_path / "query.fna"
    asm.write_text(">c\nACGT\n")
    return asm


@pytest.fixture
def two_subclades(make_db_dir, tmp_path):
    """Two built subclades of the same parent genus + a distant genus DB."""
    parent = tmp_path / "databases"
    make_db_dir(
        clade_name="Andreesenella_1", clade_taxonomy=GENUS_TAX,
        clade_rank="g", clade_rank_name="genus", n_genomes=120,
        is_subclade=True, parent_taxon="Andreesenella", subclade_id=1,
        built=True, parent=parent,
    )
    make_db_dir(
        clade_name="Andreesenella_2", clade_taxonomy=GENUS_TAX,
        clade_rank="g", clade_rank_name="genus", n_genomes=90,
        is_subclade=True, parent_taxon="Andreesenella", subclade_id=2,
        built=True, parent=parent,
    )
    return parent


class TestSubcladeFieldsLoaded:
    def test_config_get_defaults_for_non_subclade(self, Router, make_db_dir, tmp_path):
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Escherichia",
            clade_taxonomy="d__Bacteria;p__Pseudomonadota;g__Escherichia",
            clade_rank="g", clade_rank_name="genus", parent=parent,
        )
        router = Router(database_dir=parent, output_dir=tmp_path / "out")
        db = router.databases[0]
        assert db["is_subclade"] is False
        assert db["built"] is True
        assert db["parent_taxon"] is None

    def test_subclade_fields_present(self, Router, two_subclades, tmp_path):
        router = Router(database_dir=two_subclades, output_dir=tmp_path / "out")
        by_name = {d["clade_name"]: d for d in router.databases}
        d1 = by_name["Andreesenella_1"]
        assert d1["is_subclade"] is True
        assert d1["parent_taxon"] == "Andreesenella"
        assert d1["subclade_id"] == 1
        assert Path(d1["sketch_file"]).exists()


class TestMashTieBreak:
    def test_picks_min_distance_subclade(self, Router, two_subclades, tmp_path, monkeypatch):
        router = Router(database_dir=two_subclades, output_dir=tmp_path / "out")

        # Mock mash: sketch is a no-op; dist returns near for _1, far for _2.
        def fake_mash(cmd):
            if cmd[1] == "sketch":
                return ""
            # mash dist query.msh <sketch>
            sketch = cmd[3]
            if "Andreesenella_1" in sketch:
                return ("ref\tq\t0.02\t0\t900/1000\n"
                        "ref2\tq\t0.05\t0\t800/1000\n")
            return ("ref\tq\t0.30\t0\t100/1000\n"
                    "ref2\tq\t0.40\t0\t50/1000\n")
        monkeypatch.setattr(router, "_run_mash", fake_mash)

        chosen = router._route_subclade_by_mash(
            _assembly(tmp_path),
            [d for d in router.databases if d.get("is_subclade")])
        assert chosen["clade_name"] == "Andreesenella_1"

    def test_route_assembly_uses_tiebreak(self, Router, two_subclades, tmp_path, monkeypatch):
        router = Router(database_dir=two_subclades, output_dir=tmp_path / "out")

        def fake_mash(cmd):
            if cmd[1] == "sketch":
                return ""
            sketch = cmd[3]
            d = "0.02" if "Andreesenella_2" in sketch else "0.50"
            return f"ref\tq\t{d}\t0\t900/1000\n"
        monkeypatch.setattr(router, "_run_mash", fake_mash)

        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        assert decision["pipeline"] == "ReLeaf"
        assert decision["matched_database"] == "Andreesenella_2"

    def test_no_tiebreak_for_single_subclade(self, Router, make_db_dir, tmp_path, monkeypatch):
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Andreesenella_1", clade_taxonomy=GENUS_TAX,
            clade_rank="g", clade_rank_name="genus", n_genomes=50,
            is_subclade=True, parent_taxon="Andreesenella", subclade_id=1,
            built=True, parent=parent,
        )
        router = Router(database_dir=parent, output_dir=tmp_path / "out")

        called = {"n": 0}
        monkeypatch.setattr(router, "_run_mash",
                            lambda cmd: called.__setitem__("n", called["n"] + 1) or "")
        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        # Only one subclade -> no MASH comparison needed.
        assert called["n"] == 0
        assert decision["matched_database"] == "Andreesenella_1"


class TestBuiltVsUnbuilt:
    def test_unbuilt_yields_subclade_build(self, Router, make_db_dir, tmp_path, monkeypatch):
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Andreesenella_1", clade_taxonomy=GENUS_TAX,
            clade_rank="g", clade_rank_name="genus", n_genomes=120,
            is_subclade=True, parent_taxon="Andreesenella", subclade_id=1,
            built=True, parent=parent,
        )
        make_db_dir(
            clade_name="Andreesenella_2", clade_taxonomy=GENUS_TAX,
            clade_rank="g", clade_rank_name="genus", n_genomes=90,
            is_subclade=True, parent_taxon="Andreesenella", subclade_id=2,
            built=False, source_genome_dir=str(tmp_path / "raw"), parent=parent,
        )
        router = Router(database_dir=parent, output_dir=tmp_path / "out")

        def fake_mash(cmd):
            if cmd[1] == "sketch":
                return ""
            sketch = cmd[3]
            d = "0.02" if "Andreesenella_2" in sketch else "0.50"
            return f"ref\tq\t{d}\t0\t900/1000\n"
        monkeypatch.setattr(router, "_run_mash", fake_mash)

        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        assert decision["pipeline"] == "OrthoPhyl_subclade_build"
        assert decision["subclade_name"] == "Andreesenella_2"
        assert decision["parent_taxon"] == "Andreesenella"

    def test_built_yields_releaf(self, Router, two_subclades, tmp_path, monkeypatch):
        router = Router(database_dir=two_subclades, output_dir=tmp_path / "out")
        monkeypatch.setattr(router, "_run_mash",
                            lambda cmd: "" if cmd[1] == "sketch"
                            else "ref\tq\t0.02\t0\t900/1000\n")
        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        assert decision["pipeline"] == "ReLeaf"


class TestSummaryHandlesNewPipeline:
    def test_batch_summary_no_keyerror(self, Router, make_db_dir, tmp_path, monkeypatch):
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Andreesenella_1", clade_taxonomy=GENUS_TAX,
            clade_rank="g", clade_rank_name="genus", n_genomes=120,
            is_subclade=True, parent_taxon="Andreesenella", subclade_id=1,
            built=False, source_genome_dir=str(tmp_path / "raw"), parent=parent,
        )
        make_db_dir(
            clade_name="Andreesenella_2", clade_taxonomy=GENUS_TAX,
            clade_rank="g", clade_rank_name="genus", n_genomes=90,
            is_subclade=True, parent_taxon="Andreesenella", subclade_id=2,
            built=False, source_genome_dir=str(tmp_path / "raw"), parent=parent,
        )
        router = Router(database_dir=parent, output_dir=tmp_path / "out")
        monkeypatch.setattr(router, "_run_mash",
                            lambda cmd: "" if cmd[1] == "sketch"
                            else "ref\tq\t0.02\t0\t900/1000\n")
        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        # _save_decision already ran inside route_assembly without raising.
        # Exercise the batch summary path too.
        router._generate_batch_summary([decision])
        summary = (tmp_path / "out" / "batch_routing_summary.txt").read_text()
        assert "subclade build" in summary.lower()
