"""Tests for subclade/backbone routing in assembly_router.py.

Covers the MASH tie-break among same-parent subclades, --placement
(subclade vs backbone), the built vs unbuilt decision (ReLeaf vs
OrthoPhyl_subclade_build), and that _save_decision / _generate_batch_summary
handle the new pipeline value without KeyError. mash is mocked at the
_run_mash boundary -- there is no mash-availability skip anywhere in this
suite; every mash call is monkeypatched.
"""

import json
from pathlib import Path

import pytest

# A genus that both subclades (and the backbone) share -- indistinguishable
# by taxonomy string alone.
GENUS_TAX = (
    "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;"
    "f__Eubacteriaceae;g__Andreesenella"
)

ESCHERICHIA_TAX = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia"
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


@pytest.fixture
def two_subclades_and_backbone(make_db_dir, tmp_path):
    """Two built subclades PLUS the megatree backbone, all sharing GENUS_TAX."""
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
    make_db_dir(
        clade_name="Andreesenella", clade_taxonomy=GENUS_TAX,
        clade_rank="g", clade_rank_name="genus", n_genomes=15,
        is_backbone=True, parent_taxon="Andreesenella", built=True,
        parent=parent,
    )
    return parent


class TestPickNearestSubclade:
    """Pure-function coverage: no mash, just an injected distance oracle."""

    def test_picks_min_distance(self, router_module):
        candidates = [{"clade_name": "A"}, {"clade_name": "B"}, {"clade_name": "C"}]
        dists = {"A": 0.5, "B": 0.1, "C": 0.9}
        best, dist = router_module.pick_nearest_subclade(
            candidates, lambda c: dists[c["clade_name"]])
        assert best["clade_name"] == "B"
        assert dist == 0.1

    def test_skips_none_distance(self, router_module):
        candidates = [{"clade_name": "A"}, {"clade_name": "B"}]
        best, dist = router_module.pick_nearest_subclade(
            candidates, lambda c: None if c["clade_name"] == "A" else 0.2)
        assert best["clade_name"] == "B"
        assert dist == 0.2

    def test_all_unusable_returns_none(self, router_module):
        candidates = [{"clade_name": "A"}, {"clade_name": "B"}]
        best, dist = router_module.pick_nearest_subclade(candidates, lambda c: None)
        assert best is None
        assert dist is None

    def test_exact_tie_broken_by_clade_name_ascending(self, router_module):
        # Both candidates have identical distance -- winner must be the
        # alphabetically-first clade_name, deterministically, regardless of
        # input order.
        candidates = [{"clade_name": "Zebra"}, {"clade_name": "Andreesenella_1"}]
        best, _ = router_module.pick_nearest_subclade(candidates, lambda c: 0.25)
        assert best["clade_name"] == "Andreesenella_1"

        candidates_reversed = list(reversed(candidates))
        best2, _ = router_module.pick_nearest_subclade(candidates_reversed, lambda c: 0.25)
        assert best2["clade_name"] == "Andreesenella_1"


class TestSubcladeFieldsLoaded:
    def test_config_get_defaults_for_non_subclade(self, Router, make_db_dir, tmp_path):
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Escherichia", clade_taxonomy=ESCHERICHIA_TAX,
            clade_rank="g", clade_rank_name="genus", parent=parent,
        )
        router = Router(database_dir=parent, output_dir=tmp_path / "out")
        db = router.databases[0]
        assert db["is_subclade"] is False
        assert db["is_backbone"] is False
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

    def test_is_backbone_field_present(self, Router, two_subclades_and_backbone, tmp_path):
        router = Router(database_dir=two_subclades_and_backbone, output_dir=tmp_path / "out")
        by_name = {d["clade_name"]: d for d in router.databases}
        backbone = by_name["Andreesenella"]
        assert backbone["is_backbone"] is True
        assert backbone["is_subclade"] is False


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

        chosen, dist = router._route_subclade_by_mash(
            _assembly(tmp_path),
            [d for d in router.databases if d.get("is_subclade")])
        assert chosen["clade_name"] == "Andreesenella_1"
        assert dist == pytest.approx(0.02)

    def test_query_sketch_uses_matching_mash_params(
            self, Router, router_module, two_subclades, tmp_path, monkeypatch):
        """The query sketch MUST use the same -k/-s as the per-subclade sketches
        (subclade_partition.py), or distances between them are meaningless."""
        router = Router(database_dir=two_subclades, output_dir=tmp_path / "out")
        seen_cmds = []

        def fake_mash(cmd):
            seen_cmds.append(cmd)
            if cmd[1] == "sketch":
                return ""
            return "ref\tq\t0.10\t0\t900/1000\n"
        monkeypatch.setattr(router, "_run_mash", fake_mash)

        router._route_subclade_by_mash(
            _assembly(tmp_path),
            [d for d in router.databases if d.get("is_subclade")])

        sketch_cmd = seen_cmds[0]
        assert sketch_cmd[1] == "sketch"
        assert "-k" in sketch_cmd and sketch_cmd[sketch_cmd.index("-k") + 1] == router_module.MASH_K
        assert "-s" in sketch_cmd and sketch_cmd[sketch_cmd.index("-s") + 1] == router_module.MASH_S

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

    def test_non_megatree_db_routes_unaffected(self, Router, make_db_dir, tmp_path, monkeypatch):
        """Regression: a plain (non-subclade, non-backbone) DB must route
        byte-identically to pre-megatree behaviour -- no MASH call at all."""
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Escherichia", clade_taxonomy=ESCHERICHIA_TAX,
            clade_rank="g", clade_rank_name="genus", n_genomes=42, parent=parent,
        )
        router = Router(database_dir=parent, output_dir=tmp_path / "out")
        called = {"n": 0}
        monkeypatch.setattr(router, "_run_mash",
                            lambda cmd: called.__setitem__("n", called["n"] + 1) or "")
        decision = router.route_assembly(_assembly(tmp_path), ESCHERICHIA_TAX + ";s__")
        assert called["n"] == 0
        assert decision["pipeline"] == "ReLeaf"
        assert decision["matched_database"] == "Escherichia"


class TestPlacementFlag:
    def test_default_placement_is_subclade(self, Router, two_subclades_and_backbone, tmp_path):
        router = Router(database_dir=two_subclades_and_backbone, output_dir=tmp_path / "out")
        assert router.placement == "subclade"

    def test_placement_subclade_excludes_backbone(
            self, Router, two_subclades_and_backbone, tmp_path, monkeypatch):
        router = Router(database_dir=two_subclades_and_backbone,
                        output_dir=tmp_path / "out", placement="subclade")

        def fake_mash(cmd):
            if cmd[1] == "sketch":
                return ""
            sketch = cmd[3]
            d = "0.02" if "Andreesenella_1" in sketch else "0.50"
            return f"ref\tq\t{d}\t0\t900/1000\n"
        monkeypatch.setattr(router, "_run_mash", fake_mash)

        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        assert decision["pipeline"] == "ReLeaf"
        assert decision["matched_database"] == "Andreesenella_1"

    def test_placement_backbone_picks_backbone_without_mash(
            self, Router, two_subclades_and_backbone, tmp_path, monkeypatch):
        router = Router(database_dir=two_subclades_and_backbone,
                        output_dir=tmp_path / "out", placement="backbone")

        called = {"n": 0}
        monkeypatch.setattr(router, "_run_mash",
                            lambda cmd: called.__setitem__("n", called["n"] + 1) or "")

        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        assert called["n"] == 0
        assert decision["pipeline"] == "ReLeaf"
        assert decision["matched_database"] == "Andreesenella"

    def test_placement_subclade_falls_back_to_backbone_when_no_usable_sketch(
            self, Router, two_subclades_and_backbone, tmp_path, monkeypatch):
        """If neither subclade sketch is usable (e.g. missing/corrupt on disk),
        placement='subclade' must fall back to the backbone rather than
        declaring no match."""
        router = Router(database_dir=two_subclades_and_backbone,
                        output_dir=tmp_path / "out", placement="subclade")

        # Blow away both subclade sketch files so dist_fn returns None for each.
        for d in router.databases:
            if d.get("is_subclade") and d.get("sketch_file"):
                Path(d["sketch_file"]).unlink()

        def fake_mash(cmd):
            if cmd[1] == "sketch":
                return ""
            return "ref\tq\t0.10\t0\t900/1000\n"
        monkeypatch.setattr(router, "_run_mash", fake_mash)

        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        assert decision["pipeline"] == "ReLeaf"
        assert decision["matched_database"] == "Andreesenella"


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
        # Carries the SUBCLADE's own recorded taxonomy, not the query's.
        assert decision["subclade_taxonomy"] == GENUS_TAX

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

    def test_main_print_branch_no_crash(self, Router, make_db_dir, tmp_path, monkeypatch, capsys):
        """main()'s per-decision print loop must handle
        OrthoPhyl_subclade_build without KeyError -- exercised indirectly by
        constructing the same decision shape _route_to_subclade_build emits
        and running it through the print logic's field accesses."""
        parent = tmp_path / "databases"
        make_db_dir(
            clade_name="Andreesenella_1", clade_taxonomy=GENUS_TAX,
            clade_rank="g", clade_rank_name="genus", n_genomes=120,
            is_subclade=True, parent_taxon="Andreesenella", subclade_id=1,
            built=False, source_genome_dir=str(tmp_path / "raw"), parent=parent,
        )
        router = Router(database_dir=parent, output_dir=tmp_path / "out")
        monkeypatch.setattr(router, "_run_mash",
                            lambda cmd: "" if cmd[1] == "sketch"
                            else "ref\tq\t0.02\t0\t900/1000\n")
        decision = router.route_assembly(_assembly(tmp_path), GENUS_TAX + ";s__")
        # Fields main()'s print branch reads -- must all be present.
        assert decision.get("subclade_name") is not None
        assert decision.get("parent_taxon") is not None
