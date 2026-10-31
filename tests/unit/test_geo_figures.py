"""Tests for python_scripts/annotate_tree_by_country.py (--geo figures).

top_k_cross_country_pairs/pair_alphas are pure functions of an ete3 Tree +
leaf->country map, so they are tested directly with small newick strings and
cross-checked against a brute-force ete3 get_distance reference. Country
resolution / TSV parsing / leaf matching are tested on the real bundled data
files (python_scripts/data/) since those files ARE the behavior under test.
Render functions need matplotlib/ete3 and are only smoke-tested (file
written, non-empty) since pixel content isn't worth asserting on.
"""

import itertools

import pytest

ete3 = pytest.importorskip("ete3")
from ete3 import Tree  # noqa: E402


@pytest.fixture
def atbc(geo_figures_module):
    return geo_figures_module


def _tree(nw):
    return Tree(nw, format=0)


# --------------------------------------------------------------------------- #
# top_k_cross_country_pairs                                                   #
# --------------------------------------------------------------------------- #

class TestTopKCrossCountryPairs:
    def test_hand_checked_distances(self, atbc):
        nwk = "((A1:1,A2:2):1,(B1:3,B2:1):2,(C1:1,C2:1):5);"
        tree = _tree(nwk)
        leaf_country = {"A1": "USA", "A2": "USA", "B1": "Canada", "B2": "Canada",
                         "C1": "Mexico", "C2": "Mexico"}
        pairs = atbc.top_k_cross_country_pairs(tree, leaf_country, 100)
        # Hand-computed via patristic distance = sum of branch lengths on path.
        assert pairs[0] == ("A1", "B2", 5.0)
        assert pairs[1] == ("A2", "B2", 6.0)
        assert pairs[-1] == ("B1", "C2", 11.0)
        # Sorted ascending by distance.
        dists = [d for _, _, d in pairs]
        assert dists == sorted(dists)

    def test_excludes_same_country_and_unknown(self, atbc):
        nwk = "((A1:1,A2:1):1,(A3:1,U1:1):1);"
        tree = _tree(nwk)
        leaf_country = {"A1": "USA", "A2": "USA", "A3": "USA"}  # U1 absent -> Unknown
        pairs = atbc.top_k_cross_country_pairs(tree, leaf_country, 100)
        assert pairs == []

    def test_k_larger_than_available_pairs(self, atbc):
        nwk = "(A1:1,B1:1);"
        tree = _tree(nwk)
        leaf_country = {"A1": "USA", "B1": "Canada"}
        pairs = atbc.top_k_cross_country_pairs(tree, leaf_country, 1000)
        assert len(pairs) == 1
        assert pairs[0][:2] == ("A1", "B1")

    def test_k_zero_returns_empty(self, atbc):
        nwk = "(A1:1,B1:1);"
        tree = _tree(nwk)
        leaf_country = {"A1": "USA", "B1": "Canada"}
        assert atbc.top_k_cross_country_pairs(tree, leaf_country, 0) == []

    def test_only_top_k_returned(self, atbc):
        nwk = "((A1:1,A2:2):1,(B1:3,B2:1):2,(C1:1,C2:1):5);"
        tree = _tree(nwk)
        leaf_country = {"A1": "USA", "A2": "USA", "B1": "Canada", "B2": "Canada",
                         "C1": "Mexico", "C2": "Mexico"}
        pairs = atbc.top_k_cross_country_pairs(tree, leaf_country, 3)
        assert len(pairs) == 3
        assert pairs == [("A1", "B2", 5.0), ("A2", "B2", 6.0), ("A1", "B1", 7.0)]

    def test_matches_brute_force_random_tree(self, atbc):
        import random
        random.seed(12345)

        n = 40
        names = [f"L{i}" for i in range(n)]
        countries = ["USA", "Canada", "Mexico", "Brazil"]
        leaf_country = {name: random.choice(countries) for name in names}

        # Build a random-ish balanced binary tree via successive pairings.
        nodes = [f"{name}:{random.uniform(0.1, 5.0):.3f}" for name in names]
        random.shuffle(nodes)
        while len(nodes) > 1:
            new_nodes = []
            for i in range(0, len(nodes) - 1, 2):
                branch = random.uniform(0.1, 5.0)
                new_nodes.append(f"({nodes[i]},{nodes[i+1]}):{branch:.3f}")
            if len(nodes) % 2 == 1:
                new_nodes.append(nodes[-1])
            nodes = new_nodes
        newick = nodes[0] + ";"
        tree = _tree(newick)

        pairs = atbc.top_k_cross_country_pairs(tree, leaf_country, 1000000,
                                                row_chunk=7)

        brute = []
        for a, b in itertools.combinations(names, 2):
            if leaf_country[a] != leaf_country[b]:
                brute.append((a, b, tree.get_distance(a, b)))
        brute.sort(key=lambda t: t[2])

        assert len(pairs) == len(brute)
        got_dists = sorted(round(d, 6) for _, _, d in pairs)
        want_dists = sorted(round(d, 6) for _, _, d in brute)
        assert got_dists == want_dists


# --------------------------------------------------------------------------- #
# pair_alphas                                                                 #
# --------------------------------------------------------------------------- #

class TestPairAlphas:
    def test_empty(self, atbc):
        assert atbc.pair_alphas([]) == []

    def test_kth_pair_is_zero_closest_is_one(self, atbc):
        alphas = atbc.pair_alphas([1.0, 2.0, 4.0])
        assert alphas[0] == pytest.approx(1.0 - 1.0 / 4.0)
        assert alphas[-1] == pytest.approx(0.0)

    def test_monotonic(self, atbc):
        dists = [0.0, 1.0, 2.0, 3.0, 10.0]
        alphas = atbc.pair_alphas(dists)
        assert alphas == sorted(alphas, reverse=True)

    def test_all_zero_distances_gives_alpha_one(self, atbc):
        assert atbc.pair_alphas([0.0, 0.0, 0.0]) == [1.0, 1.0, 1.0]

    def test_single_nonzero_pair_gives_alpha_zero(self, atbc):
        # A lone pair is simultaneously the closest AND the most-distant
        # (K-th) pair in the set, so the "K-th pair is transparent" rule
        # makes it alpha 0 -- not 1. Only an all-zero-distance set (d_max<=0)
        # gets the alpha-1 special case (see test_all_zero_distances_...).
        assert atbc.pair_alphas([5.0]) == [0.0]


# --------------------------------------------------------------------------- #
# normalize_country / country resolution                                      #
# --------------------------------------------------------------------------- #

class TestNormalizeCountry:
    def test_plain_name(self, atbc):
        assert atbc.normalize_country("Viet Nam") == "Viet Nam"

    def test_underscore_region_encoding(self, atbc):
        assert atbc.normalize_country("USA:_Montana") == "USA"
        assert atbc.normalize_country("USA: Montana") == "USA"

    def test_unknown_values(self, atbc):
        for v in ("", "NA", "na", "Unknown", "none", "NULL", "-"):
            assert atbc.normalize_country(v) == atbc.UNKNOWN

    def test_quoted_value(self, atbc):
        assert atbc.normalize_country('"Viet Nam"') == "Viet Nam"


class TestCountryCoordResolution:
    def test_every_alias_row_resolves(self, atbc):
        """Every alias in country_aliases.tsv should resolve to some lat/lon,
        either directly (latlon rows) or indirectly via a country name that
        exists in the bundled Natural Earth label points."""
        _, name_to_latlon = atbc.load_country_coords()
        aliases = atbc._load_aliases()
        unresolved = []
        for alias_name, (kind, value) in aliases.items():
            coord = atbc.resolve_country_latlon(alias_name, name_to_latlon)
            if coord is None:
                unresolved.append((alias_name, kind, value))
        assert unresolved == [], f"unresolved aliases: {unresolved}"

    def test_usa_resolves(self, atbc):
        _, name_to_latlon = atbc.load_country_coords()
        assert atbc.resolve_country_latlon("USA", name_to_latlon) is not None

    def test_unknown_returns_none(self, atbc):
        _, name_to_latlon = atbc.load_country_coords()
        assert atbc.resolve_country_latlon(atbc.UNKNOWN, name_to_latlon) is None
        assert atbc.resolve_country_latlon("", name_to_latlon) is None

    def test_nonexistent_country_returns_none(self, atbc):
        _, name_to_latlon = atbc.load_country_coords()
        assert atbc.resolve_country_latlon("Narnia", name_to_latlon) is None


# --------------------------------------------------------------------------- #
# parse_geo_tsv                                                               #
# --------------------------------------------------------------------------- #

class TestParseGeoTsv:
    def test_with_header(self, atbc, tmp_path):
        tsv = tmp_path / "geo.tsv"
        tsv.write_text("genome_id\tlocation\nGCF_000001\tUSA\nGCF_000002\tCanada\n")
        result = atbc.parse_geo_tsv(str(tsv))
        assert result == {"GCF_000001": "USA", "GCF_000002": "Canada"}

    def test_with_comments_and_blank_lines(self, atbc, tmp_path):
        tsv = tmp_path / "geo.tsv"
        tsv.write_text("# a comment\n\nGCF_000001\tUSA\n# another comment\nGCF_000002\tCanada\n")
        result = atbc.parse_geo_tsv(str(tsv))
        assert result == {"GCF_000001": "USA", "GCF_000002": "Canada"}

    def test_blank_location_is_unknown(self, atbc, tmp_path):
        tsv = tmp_path / "geo.tsv"
        tsv.write_text("GCF_000001\t\n")
        result = atbc.parse_geo_tsv(str(tsv))
        assert result["GCF_000001"] == atbc.UNKNOWN

    def test_malformed_lines_skipped(self, atbc, tmp_path):
        tsv = tmp_path / "geo.tsv"
        tsv.write_text("GCF_000001\tUSA\nno_tab_here\nGCF_000002\tCanada\n")
        result = atbc.parse_geo_tsv(str(tsv))
        assert result == {"GCF_000001": "USA", "GCF_000002": "Canada"}


class TestBuildLeafCountryMap:
    def test_tsv_overrides_metadata(self, atbc):
        acc_country = {"GCF_000001": "USA"}
        geo_tsv_country = {"GCF_000001": "Canada"}
        leaf_country, missing = atbc.build_leaf_country_map(
            ["GCF_000001"], acc_country, geo_tsv_country)
        assert leaf_country["GCF_000001"] == "Canada"
        assert missing == []

    def test_metadata_only(self, atbc):
        acc_country = {"GCF_000001": "USA"}
        leaf_country, missing = atbc.build_leaf_country_map(
            ["GCF_000001"], acc_country, {})
        assert leaf_country["GCF_000001"] == "USA"

    def test_unmatched_leaf_is_unknown_and_missing(self, atbc):
        leaf_country, missing = atbc.build_leaf_country_map(
            ["GCF_999999"], {}, {})
        assert leaf_country["GCF_999999"] == atbc.UNKNOWN
        assert missing == ["GCF_999999"]

    def test_stem_matching(self, atbc):
        """Leaf name carries a .fna extension; metadata key doesn't."""
        acc_country = {"GCF_000001": "USA"}
        leaf_country, missing = atbc.build_leaf_country_map(
            ["GCF_000001.fna"], acc_country, {})
        assert leaf_country["GCF_000001.fna"] == "USA"
        assert missing == []

    def test_version_stripped_matching(self, atbc):
        """Leaf name carries an NCBI version suffix; metadata key doesn't."""
        acc_country = {"GCF_000001": "USA"}
        leaf_country, missing = atbc.build_leaf_country_map(
            ["GCF_000001.1"], acc_country, {})
        assert leaf_country["GCF_000001.1"] == "USA"
        assert missing == []

    def test_reverse_version_matching(self, atbc):
        """Metadata key carries a version suffix; leaf name doesn't."""
        acc_country = {"GCF_000001.1": "USA"}
        leaf_country, missing = atbc.build_leaf_country_map(
            ["GCF_000001"], acc_country, {})
        assert leaf_country["GCF_000001"] == "USA"
        assert missing == []


# --------------------------------------------------------------------------- #
# parse_metadata                                                              #
# --------------------------------------------------------------------------- #

class TestParseMetadata:
    def test_basic_columns(self, atbc, tmp_path):
        meta = tmp_path / "all_asm_acc_metadata"
        meta.write_text(
            "GCF_000001 species1 strain1 biosample1 1234 ok ok USA:_Montana\n"
            "GCF_000002 species2 strain2 biosample2 5678 ok ok Canada\n"
        )
        acc_country, acc_species = atbc.parse_metadata(
            str(meta), acc_col=1, country_col=8, species_col=3, want_species=False)
        assert acc_country == {"GCF_000001": "USA", "GCF_000002": "Canada"}
        assert acc_species == {}

    def test_want_species(self, atbc, tmp_path):
        meta = tmp_path / "all_asm_acc_metadata"
        meta.write_text("GCF_000001 species1 strain1 biosample1 1234 ok ok USA\n")
        acc_country, acc_species = atbc.parse_metadata(
            str(meta), acc_col=1, country_col=8, species_col=3, want_species=True)
        assert acc_species == {"GCF_000001": "strain1"}

    def test_short_lines_skipped(self, atbc, tmp_path):
        meta = tmp_path / "all_asm_acc_metadata"
        meta.write_text("GCF_000001 too short\nGCF_000002 species2 strain2 bios2 5678 ok ok Canada\n")
        acc_country, _ = atbc.parse_metadata(
            str(meta), acc_col=1, country_col=8, species_col=3, want_species=False)
        assert acc_country == {"GCF_000002": "Canada"}


# --------------------------------------------------------------------------- #
# _jitter                                                                     #
# --------------------------------------------------------------------------- #

class TestJitter:
    def test_offsets_stay_within_disk(self, atbc):
        import math

        radius = 1.5
        for i in range(200):
            name = f"genome_{i}"
            lat, lon = atbc._jitter(name, 10.0, 20.0, max_radius=radius)
            dist = math.hypot(lat - 10.0, lon - 20.0)
            assert dist <= radius + 1e-9

    def test_deterministic_per_name(self, atbc):
        a1 = atbc._jitter("GCF_000001", 10.0, 20.0, max_radius=1.5)
        a2 = atbc._jitter("GCF_000001", 10.0, 20.0, max_radius=1.5)
        assert a1 == a2

    def test_different_names_differ(self, atbc):
        a = atbc._jitter("GCF_000001", 10.0, 20.0, max_radius=1.5)
        b = atbc._jitter("GCF_000002", 10.0, 20.0, max_radius=1.5)
        assert a != b

    def test_smaller_radius_shrinks_offsets(self, atbc):
        import math

        for i in range(50):
            name = f"genome_{i}"
            lat_small, lon_small = atbc._jitter(name, 0.0, 0.0, max_radius=0.1)
            dist_small = math.hypot(lat_small, lon_small)
            assert dist_small <= 0.1 + 1e-9


# --------------------------------------------------------------------------- #
# Render smoke tests                                                          #
# --------------------------------------------------------------------------- #

class TestRenderSmoke:
    def test_render_map_writes_nonempty_file(self, atbc, tmp_path):
        pytest.importorskip("matplotlib")
        land_polygons, name_to_latlon = atbc.load_country_coords()
        color_map = atbc.build_color_map({"USA", "Canada"})
        leaf_country = {"A1": "USA", "B1": "Canada"}
        counts = {"USA": 1, "Canada": 1}
        pairs = [("A1", "B1", 5.0)]
        alphas = [1.0]
        out = tmp_path / "map.png"
        atbc.render_map(str(out), leaf_country, counts, pairs, alphas, color_map,
                        land_polygons, name_to_latlon, "png", 72, title="test")
        assert out.exists()
        assert out.stat().st_size > 0

    def test_render_tree_writes_nonempty_file(self, atbc, tmp_path):
        tree_path = tmp_path / "t.nwk"
        tree_path.write_text("(A1:1,B1:1);")
        color_map = atbc.build_color_map({"USA", "Canada"})
        out = tmp_path / "tree.png"
        counts, missing = atbc.render_tree(
            str(tree_path), str(out), {"A1": "USA", "B1": "Canada"}, {}, color_map,
            False, "png", 72, title="test")
        assert out.exists()
        assert out.stat().st_size > 0
        assert counts == {"USA": 1, "Canada": 1}
        assert missing == []
