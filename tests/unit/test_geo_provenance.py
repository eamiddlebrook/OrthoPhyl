"""Tests for python_scripts/geo_provenance.py.

geo_provenance.py flags cases where all_asm_acc_metadata's col-8 geo_loc_name
country might actually be the sequencing/submitting lab's own country, by
comparing it against a country resolved from sequenced_by / asm_submitter /
biosample_owner (cols 9-12), and separately flags cases where the BioSample
lat_lon (col 13) disagrees with col 8 or lands on the submitter org's own ROR
coordinates. It is flag-only: these tests pin the exact flag matrices
(geo_vs_submitter, latlon_vs_geo, latlon_vs_submitter) and the resolution
heuristics (override table, free-text scanning, point-in-polygon, ROR), not
just that "something" gets written.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "python_scripts" / "geo_provenance.py"


@pytest.fixture
def gp(geo_provenance_module):
    return geo_provenance_module


@pytest.fixture(scope="module")
def index():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "geo_provenance_for_index", str(SCRIPT))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.load_country_index()


@pytest.fixture
def overrides(gp):
    return gp.load_overrides()


# A 13-column all_asm_acc_metadata row, matching project_metadata()'s layout:
#   1 acc 2 species 3 species.strain 4 biosample 5 taxid 6 taxcheck 7 excl
#   8 geo  9 collected_by 10 sequenced_by 11 biosample_owner 12 asm_submitter
#   13 lat_lon
def row(geo="NA", collected_by="NA", sequenced_by="NA", biosample_owner="NA",
        asm_submitter="NA", lat_lon="NA"):
    return [
        "GCF_1.1", "Esch_coli", "Esch_coli.K12", "SAMN001", "562", "OK", "na",
        geo, collected_by, sequenced_by, biosample_owner, asm_submitter,
        lat_lon,
    ]


# --------------------------------------------------------------------------- #
# Country index / canonicalization                                           #
# --------------------------------------------------------------------------- #

class TestLoadCountryIndex:
    def test_builds_without_error(self, gp):
        idx = gp.load_country_index()
        assert "United Kingdom" in idx.names
        assert idx.name_lookup["japan"] == "Japan"

    def test_alias_chain_resolves_to_a_canonical_name(self, gp):
        """South Georgia... etc are latlon-kind and must not appear as a
        resolvable country alias; USA/United Kingdom (country-kind) must."""
        idx = gp.load_country_index()
        assert idx.alias_lookup["usa"] == "United States of America"
        assert idx.alias_lookup["uk"] == "United Kingdom"


class TestCanonicalCountry:
    def test_known_aliases_canonicalize_to_the_same_name(self, gp, index):
        assert gp.canonical_country("USA", index) == "United States of America"
        assert (
            gp.canonical_country("United_States", index)
            == "United States of America"
        )
        assert gp.canonical_country("United_Kingdom", index) == "United Kingdom"

    def test_direct_ne_name_with_underscore(self, gp, index):
        assert gp.canonical_country("Japan", index) == "Japan"

    def test_missing_unknown_and_none_return_none(self, gp, index):
        assert gp.canonical_country(None, index) is None
        assert gp.canonical_country("", index) is None
        assert gp.canonical_country("NA", index) is None
        assert gp.canonical_country("Unknown", index) is None

    def test_unrecognized_country_returns_none(self, gp, index):
        assert gp.canonical_country("Neverland", index) is None


# --------------------------------------------------------------------------- #
# Free-text country extraction                                               #
# --------------------------------------------------------------------------- #

class TestCountryFromText:
    def test_finds_country_at_end_of_affiliation(self, gp, index):
        text = "Dept of Microbiology, University of Tokyo, Japan"
        assert gp.country_from_text(text, index) == "Japan"

    def test_prefers_last_occurrence(self, gp, index):
        text = "Formerly in Japan, now based in France"
        assert gp.country_from_text(text, index) == "France"

    def test_georgia_is_too_ambiguous_to_match(self, gp, index):
        """Georgia collides with the US state; must never match from text."""
        assert gp.country_from_text("University of Georgia", index) is None

    def test_georgia_does_not_block_a_real_country_later_in_the_text(self, gp, index):
        text = "CDC, Atlanta, Georgia, USA"
        assert gp.country_from_text(text, index) == "United States of America"

    def test_new_mexico_guard(self, gp, index):
        """New Mexico State University must not match the country Mexico."""
        assert gp.country_from_text("New Mexico State University", index) is None

    def test_mexico_matches_without_the_guard_word(self, gp, index):
        assert gp.country_from_text("Instituto in Mexico", index) == "Mexico"

    def test_empty_or_none_returns_none(self, gp, index):
        assert gp.country_from_text("", index) is None
        assert gp.country_from_text(None, index) is None
        assert gp.country_from_text("NA", index) is None

    def test_no_country_mentioned_returns_none(self, gp, index):
        assert gp.country_from_text("Department of Biology", index) is None


# --------------------------------------------------------------------------- #
# Override table                                                             #
# --------------------------------------------------------------------------- #

class TestLoadOverrides:
    def test_default_table_has_known_entries(self, gp, overrides):
        assert gp.match_overrides("PulseNet USA network", overrides) == "USA"
        assert gp.match_overrides("Centers for Disease Control", overrides) == "USA"
        assert gp.match_overrides("Wellcome Sanger Institute", overrides) == "United Kingdom"

    def test_custom_path_parses_tabs_and_skips_comments(self, gp, tmp_path):
        f = tmp_path / "overrides.tsv"
        f.write_text("# a comment\nFoo Lab\tGermany\nbad line with no tab\n")
        overrides = gp.load_overrides(str(f))
        assert overrides == [("foo lab", "Germany")]

    def test_missing_path_returns_empty_list(self, gp, tmp_path):
        assert gp.load_overrides(str(tmp_path / "nope.tsv")) == []


class TestMatchOverrides:
    def test_longest_substring_wins(self, gp):
        overrides = [("abcdef", "Y"), ("abc", "X")]
        assert gp.match_overrides("xx abcdef xx", overrides) == "Y"
        assert gp.match_overrides("xx abc only xx", overrides) == "X"

    def test_no_match_returns_none(self, gp):
        assert gp.match_overrides("some random lab", [("abc", "X")]) is None

    def test_none_text_returns_none(self, gp):
        assert gp.match_overrides(None, [("abc", "X")]) is None


# --------------------------------------------------------------------------- #
# is_usable_text / resolve_submitter_country                                 #
# --------------------------------------------------------------------------- #

class TestIsUsableText:
    @pytest.mark.parametrize("value", ["NA", "", "Unknown", "none", "null", "-", None])
    def test_unusable_values(self, gp, value):
        assert gp.is_usable_text(value) is False

    def test_usable_value(self, gp):
        assert gp.is_usable_text("CDC") is True


class TestResolveSubmitterCountry:
    def test_sequenced_by_override_wins_first(self, gp, index, overrides):
        country, source = gp.resolve_submitter_country(
            "PulseNet", "NA", "NA", index, overrides)
        assert country == "United States of America"
        assert source == "sequenced_by"

    def test_falls_back_to_asm_submitter(self, gp, index, overrides):
        country, source = gp.resolve_submitter_country(
            "NA", "University of Tokyo, Japan", "NA", index, overrides)
        assert country == "Japan"
        assert source == "asm_submitter"

    def test_falls_back_to_biosample_owner(self, gp, index, overrides):
        country, source = gp.resolve_submitter_country(
            "NA", "NA", "National Institute, Tokyo, Japan", index, overrides)
        assert country == "Japan"
        assert source == "biosample_owner"

    def test_nothing_resolves_returns_none_none(self, gp, index, overrides):
        assert gp.resolve_submitter_country(
            "NA", "NA", "NA", index, overrides) == (None, None)

    def test_free_text_scan_used_when_no_override_matches(self, gp, index, overrides):
        country, source = gp.resolve_submitter_country(
            "Institut Jacques Monod, France", "NA", "NA", index, overrides)
        assert country == "France"
        assert source == "sequenced_by"

    def test_priority_order_sequenced_by_beats_later_fields(self, gp, index, overrides):
        """If sequenced_by resolves, asm_submitter/biosample_owner are ignored
        even if they'd resolve to a different country."""
        country, source = gp.resolve_submitter_country(
            "PulseNet", "Some Lab, Japan", "Some Owner, France", index, overrides)
        assert country == "United States of America"
        assert source == "sequenced_by"


# --------------------------------------------------------------------------- #
# geo_vs_submitter flag matrix                                               #
# --------------------------------------------------------------------------- #

class TestGeoVsSubmitter:
    def test_match(self, gp, index):
        assert gp.geo_vs_submitter(
            "USA", "United States of America", index) == "match"

    def test_differ(self, gp, index):
        assert gp.geo_vs_submitter(
            "Japan", "United States of America", index) == "differ"

    def test_no_geo(self, gp, index):
        assert gp.geo_vs_submitter(
            "NA", "United States of America", index) == "no_geo"
        assert gp.geo_vs_submitter("Unknown", "Japan", index) == "no_geo"

    def test_unresolved(self, gp, index):
        assert gp.geo_vs_submitter("USA", None, index) == "unresolved"


# --------------------------------------------------------------------------- #
# lat_lon parsing                                                            #
# --------------------------------------------------------------------------- #

class TestParseLatLon:
    def test_insdc_format(self, gp):
        assert gp.parse_lat_lon("38.98_N_77.11_W") == (38.98, -77.11)

    def test_decimal_comma_format(self, gp):
        assert gp.parse_lat_lon("38.98,_-77.11") == (38.98, -77.11)

    def test_decimal_space_format(self, gp):
        assert gp.parse_lat_lon("38.98_-77.11") == (38.98, -77.11)

    def test_southern_western_hemisphere(self, gp):
        assert gp.parse_lat_lon("33.87_S_151.21_E") == (-33.87, 151.21)

    @pytest.mark.parametrize("value", ["NA", "", "Unknown", None, "not_collected"])
    def test_missing_sentinels_return_none(self, gp, value):
        assert gp.parse_lat_lon(value) is None

    def test_out_of_range_lat_returns_none(self, gp):
        assert gp.parse_lat_lon("999_N_77_W") is None

    def test_out_of_range_lon_returns_none(self, gp):
        assert gp.parse_lat_lon("38_N_999_W") is None

    def test_garbage_returns_none(self, gp):
        assert gp.parse_lat_lon("not_a_coordinate") is None


# --------------------------------------------------------------------------- #
# haversine_km                                                               #
# --------------------------------------------------------------------------- #

class TestHaversineKm:
    def test_same_point_is_zero(self, gp):
        assert gp.haversine_km(35.0, 139.0, 35.0, 139.0) == 0.0

    def test_known_distance_tokyo_to_gifu(self, gp):
        # Tokyo (35.6895, 139.6917) to Gifu Univ (35.42291, 136.76039).
        d = gp.haversine_km(35.6895, 139.6917, 35.42291, 136.76039)
        assert 260 < d < 275


# --------------------------------------------------------------------------- #
# country_at (point-in-polygon + coastal fallback)                          #
# --------------------------------------------------------------------------- #

class TestCountryAt:
    def test_tokyo_resolves_to_japan(self, gp):
        assert gp.country_at(35.6895, 139.6917) == "Japan"

    def test_pune_resolves_to_india(self, gp):
        assert gp.country_at(18.5204, 73.8567) == "India"

    def test_mid_atlantic_is_no_country(self, gp):
        assert gp.country_at(30.0, -40.0) is None

    def test_mid_pacific_is_no_country(self, gp):
        assert gp.country_at(0.0, -160.0) is None

    def test_coastal_point_resolves_via_fallback(self, gp):
        """(35.0, 140.3) sits just offshore Japan -- misses every 110m
        polygon directly, but resolves via the <=50km nearest-vertex
        fallback rather than being called open ocean."""
        countries = gp.load_country_polygons()
        direct_hit = any(
            gp._point_in_polygon(35.0, 140.3, poly)
            for c in countries for poly in c["polygons"]
        )
        assert direct_hit is False
        assert gp.country_at(35.0, 140.3) == "Japan"

    def test_far_offshore_stays_no_country(self, gp):
        assert gp.country_at(35.0, 145.0) is None


# --------------------------------------------------------------------------- #
# resolve_latlon_country / latlon_vs_geo                                    #
# --------------------------------------------------------------------------- #

class TestResolveLatlonCountry:
    def test_resolves_coordinates_to_country(self, gp):
        assert gp.resolve_latlon_country("35.6895_N_139.6917_E") == (
            35.6895, 139.6917, "Japan")

    def test_resolves_to_no_country(self, gp):
        assert gp.resolve_latlon_country("30_N_40_W") == (30.0, -40.0, None)

    def test_unparseable_returns_all_none(self, gp):
        assert gp.resolve_latlon_country("NA") == (None, None, None)


class TestLatlonVsGeo:
    def test_match(self, gp, index):
        assert gp.latlon_vs_geo(
            "35.6895_N_139.6917_E", "Japan", index) == "match"

    def test_differ(self, gp, index):
        assert gp.latlon_vs_geo(
            "35.6895_N_139.6917_E", "USA", index) == "differ"

    def test_no_country_when_latlon_is_open_ocean(self, gp, index):
        assert gp.latlon_vs_geo("30_N_40_W", "USA", index) == "no_country"

    def test_no_geo_when_col8_is_unusable(self, gp, index):
        assert gp.latlon_vs_geo(
            "35.6895_N_139.6917_E", "NA", index) == "no_geo"

    def test_NA_when_latlon_is_missing(self, gp, index):
        assert gp.latlon_vs_geo("NA", "USA", index) == "NA"


# --------------------------------------------------------------------------- #
# ROR lookup                                                                 #
# --------------------------------------------------------------------------- #

def _ror_payload(country, lat, lng, chosen=True):
    return {
        "items": [{
            "chosen": chosen,
            "organization": {
                "locations": [{
                    "geonames_details": {
                        "country_name": country, "lat": lat, "lng": lng,
                    }
                }]
            },
        }]
    }


class TestRorLookup:
    def test_chosen_match_resolves(self, gp):
        fetch = lambda org: _ror_payload("Japan", 35.42291, 136.76039)
        result = gp.ror_lookup("Gifu University", {}, fetch=fetch)
        assert result == ("Japan", 35.42291, 136.76039)

    def test_non_chosen_match_is_ignored(self, gp):
        fetch = lambda org: _ror_payload("Japan", 1.0, 1.0, chosen=False)
        assert gp.ror_lookup("Mystery Org", {}, fetch=fetch) == (None, None, None)

    def test_no_items_resolves_to_none(self, gp):
        assert gp.ror_lookup("Mystery Org", {}, fetch=lambda org: {"items": []}) == (
            None, None, None)

    def test_fetch_exception_is_swallowed(self, gp):
        def boom(org):
            raise RuntimeError("network down")
        assert gp.ror_lookup("Some Org", {}, fetch=boom) == (None, None, None)

    def test_unusable_org_text_skips_fetch_entirely(self, gp):
        calls = []
        fetch = lambda org: calls.append(org) or _ror_payload("Japan", 1, 1)
        assert gp.ror_lookup("NA", {}, fetch=fetch) == (None, None, None)
        assert calls == []

    def test_cache_hit_avoids_a_second_fetch(self, gp):
        calls = []
        def fetch(org):
            calls.append(org)
            return _ror_payload("Japan", 35.42291, 136.76039)
        cache = {}
        gp.ror_lookup("Gifu University", cache, fetch=fetch)
        gp.ror_lookup("Gifu University", cache, fetch=fetch)
        assert len(calls) == 1

    def test_cached_miss_is_also_not_refetched(self, gp):
        calls = []
        def fetch(org):
            calls.append(org)
            return {"items": []}
        cache = {}
        gp.ror_lookup("Mystery Org", cache, fetch=fetch)
        gp.ror_lookup("Mystery Org", cache, fetch=fetch)
        assert len(calls) == 1


class TestResolveSubmitterRor:
    def test_prefers_asm_submitter_over_biosample_owner_and_sequenced_by(self, gp):
        fetch = lambda org: (
            _ror_payload("Japan", 35.42291, 136.76039) if org == "asm submitter org"
            else _ror_payload("France", 1.0, 1.0)
        )
        country, lat, lng, source = gp.resolve_submitter_ror(
            "sequenced by org", "asm_submitter_org", "biosample_owner_org",
            {}, fetch=fetch)
        assert source == "asm_submitter"
        assert country == "Japan"

    def test_falls_back_when_asm_submitter_is_unusable(self, gp):
        fetch = lambda org: _ror_payload("France", 1.0, 1.0)
        country, lat, lng, source = gp.resolve_submitter_ror(
            "NA", "NA", "biosample_owner_org", {}, fetch=fetch)
        assert source == "biosample_owner"
        assert country == "France"

    def test_nothing_resolves_returns_all_none(self, gp):
        fetch = lambda org: {"items": []}
        assert gp.resolve_submitter_ror("NA", "NA", "NA", {}, fetch=fetch) == (
            None, None, None, None)


# --------------------------------------------------------------------------- #
# latlon_vs_submitter flag matrix                                           #
# --------------------------------------------------------------------------- #

class TestLatlonVsSubmitter:
    def test_at_submitter_when_within_threshold(self, gp):
        dist, flag = gp.latlon_vs_submitter(
            35.42291, 136.76039, "Japan", 35.42291, 136.76039, "Japan", 25.0)
        assert flag == "at_submitter"
        assert dist == 0.0

    def test_same_country_when_beyond_threshold_but_same_country(self, gp):
        dist, flag = gp.latlon_vs_submitter(
            35.6895, 139.6917, "Japan", 35.42291, 136.76039, "Japan", 25.0)
        assert flag == "same_country"
        assert dist > 25.0

    def test_differ_when_different_country(self, gp):
        dist, flag = gp.latlon_vs_submitter(
            35.6895, 139.6917, "Japan", 38.98, -77.11,
            "United States of America", 25.0)
        assert flag == "differ"
        assert dist > 25.0

    def test_unresolved_when_latlon_missing(self, gp):
        assert gp.latlon_vs_submitter(
            None, None, None, 38.98, -77.11, "United States of America", 25.0
        ) == (None, "unresolved")

    def test_unresolved_when_submitter_coords_missing(self, gp):
        assert gp.latlon_vs_submitter(
            35.0, 140.0, "Japan", None, None, None, 25.0
        ) == (None, "unresolved")

    def test_custom_threshold(self, gp):
        # Same two points as the "same_country" case above (~267 km apart),
        # but with a threshold wide enough to swallow it into at_submitter.
        dist, flag = gp.latlon_vs_submitter(
            35.6895, 139.6917, "Japan", 35.42291, 136.76039, "Japan", 300.0)
        assert flag == "at_submitter"


# --------------------------------------------------------------------------- #
# annotate_rows                                                              #
# --------------------------------------------------------------------------- #

class TestAnnotateRows:
    def test_match_row(self, gp, index, overrides):
        r = row(geo="USA", sequenced_by="PulseNet")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert len(extended) == 19
        assert extended[:13] == r
        assert extended[13] == "United_States_of_America"
        assert extended[14] == "match"
        assert rec["submitter_source"] == "sequenced_by"

    def test_differ_row(self, gp, index, overrides):
        r = row(geo="Japan", sequenced_by="PulseNet")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert extended[13] == "United_States_of_America"
        assert extended[14] == "differ"

    def test_submitter_country_is_sanitized_for_the_whitespace_split_output(
            self, gp, index, overrides):
        """A multi-word canonical country name must not create extra columns
        in the whitespace-split output row."""
        r = row(geo="USA", sequenced_by="PulseNet")
        extended, _ = next(gp.annotate_rows([r], index, overrides))
        assert " " not in extended[13]
        assert len(" ".join(extended).split()) == 19

    def test_no_geo_row(self, gp, index, overrides):
        r = row(geo="NA", sequenced_by="PulseNet")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert extended[14] == "no_geo"

    def test_unresolved_row(self, gp, index, overrides):
        r = row(geo="USA")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert extended[13] == "NA"
        assert extended[14] == "unresolved"

    def test_short_row_is_padded_with_NA(self, gp, index, overrides):
        """A row missing columns 9-13 (e.g. an older all_asm_acc_metadata) is
        padded rather than raising."""
        short = ["GCF_1.1", "Esch_coli", "Esch_coli.K12", "SAMN001",
                 "562", "OK", "na", "USA"]
        extended, rec = next(gp.annotate_rows([short], index, overrides))
        assert len(extended) == 19
        assert rec["collected_by"] == "NA"
        assert rec["sequenced_by"] == "NA"
        assert rec["biosample_owner"] == "NA"
        assert rec["asm_submitter"] == "NA"
        assert rec["lat_lon"] == "NA"
        assert extended[14] == "unresolved"
        assert extended[16] == "NA"  # latlon_vs_geo

    def test_multiple_rows_preserve_order(self, gp, index, overrides):
        rows_in = [row(geo="USA", sequenced_by="PulseNet"),
                   row(geo="Japan", sequenced_by="PulseNet")]
        out = list(gp.annotate_rows(rows_in, index, overrides))
        assert len(out) == 2
        assert out[0][0][14] == "match"
        assert out[1][0][14] == "differ"

    def test_latlon_columns_populated_without_ror(self, gp, index, overrides):
        r = row(geo="Japan", lat_lon="35.6895_N_139.6917_E")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert extended[15] == "Japan"      # latlon_country
        assert extended[16] == "match"      # latlon_vs_geo
        assert extended[17] == "NA"         # submitter_km (no --ror)
        assert extended[18] == "unresolved"  # latlon_vs_submitter

    def test_ror_disabled_never_calls_fetch(self, gp, index, overrides):
        calls = []
        r = row(geo="Japan", sequenced_by="PulseNet",
                lat_lon="35.6895_N_139.6917_E")
        extended, rec = next(gp.annotate_rows(
            [r], index, overrides, ror_enabled=False,
            ror_fetch=lambda org: calls.append(org)))
        assert calls == []
        assert extended[17] == "NA"
        assert extended[18] == "unresolved"

    def test_ror_enabled_computes_submitter_km_and_flag(self, gp, index, overrides):
        fetch = lambda org: _ror_payload("Japan", 35.42291, 136.76039)
        r = row(geo="Japan", biosample_owner="University_of_Tokyo,_Japan",
                lat_lon="35.6895_N_139.6917_E")
        extended, rec = next(gp.annotate_rows(
            [r], index, overrides, ror_enabled=True, ror_cache={},
            ror_fetch=fetch))
        assert rec["submitter_km"] == "266.8"
        assert rec["latlon_vs_submitter"] == "same_country"

    def test_custom_at_submitter_km_threshold(self, gp, index, overrides):
        fetch = lambda org: _ror_payload("Japan", 35.42291, 136.76039)
        r = row(geo="Japan", biosample_owner="University_of_Tokyo,_Japan",
                lat_lon="35.6895_N_139.6917_E")
        extended, rec = next(gp.annotate_rows(
            [r], index, overrides, ror_enabled=True, ror_cache={},
            ror_fetch=fetch, at_submitter_km=300.0))
        assert rec["latlon_vs_submitter"] == "at_submitter"


# --------------------------------------------------------------------------- #
# CLI                                                                        #
# --------------------------------------------------------------------------- #

def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + list(args),
        capture_output=True, text=True,
    )


class TestCLI:
    def test_help(self):
        assert run_cli("--help").returncode == 0

    def test_annotate_end_to_end(self, tmp_path):
        meta = tmp_path / "all_asm_acc_metadata"
        meta.write_text("\n".join([
            " ".join(row(geo="USA", sequenced_by="PulseNet")),
            " ".join(row(geo="Japan", sequenced_by="PulseNet")),
            " ".join(row(geo="NA")),
        ]) + "\n")
        out = tmp_path / "out"
        report = tmp_path / "report.tsv"

        r = run_cli("annotate", "--meta", str(meta),
                    "--out", str(out), "--report", str(report))
        assert r.returncode == 0, r.stderr
        assert "19 columns" in r.stderr
        assert "geo_vs_submitter match: 1" in r.stderr
        assert "geo_vs_submitter differ: 1" in r.stderr

        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3
        for ln in lines:
            assert len(ln.split()) == 19

        report_lines = report.read_text().strip().split("\n")
        assert len(report_lines) == 4  # header + 3 rows
        header = report_lines[0].split("\t")
        assert header[-1] == "latlon_vs_submitter"
        assert header[-2] == "submitter_km"
        assert header[-3] == "latlon_vs_geo"
        assert header[-4] == "latlon_country"
        assert header[-5] == "geo_vs_submitter"
        assert header[-6] == "submitter_source"
        assert header[-7] == "submitter_country"

    def test_lab_default_fingerprint_reported(self, tmp_path):
        """5+ samples from the same submitter, all matching their own
        country, are surfaced as a possible lab-default geolocation."""
        meta = tmp_path / "all_asm_acc_metadata"
        meta.write_text("\n".join(
            " ".join(row(geo="USA", sequenced_by="PulseNet")) for _ in range(5)
        ) + "\n")
        out = tmp_path / "out"
        report = tmp_path / "report.tsv"

        r = run_cli("annotate", "--meta", str(meta),
                    "--out", str(out), "--report", str(report))
        assert r.returncode == 0, r.stderr
        assert "lab-default" in r.stderr
        assert "PulseNet" in r.stderr
        assert "sequenced_by" in r.stderr

    def test_custom_overrides_file(self, tmp_path):
        meta = tmp_path / "all_asm_acc_metadata"
        # sequenced_by arrives already sanitize_field()-clean (underscores,
        # not spaces) -- match_overrides/country_from_text undo that before
        # matching free text.
        meta.write_text(
            " ".join(row(geo="Germany", sequenced_by="Foo_Lab")) + "\n")
        overrides = tmp_path / "overrides.tsv"
        overrides.write_text("Foo Lab\tGermany\n")
        out = tmp_path / "out"
        report = tmp_path / "report.tsv"

        r = run_cli("annotate", "--meta", str(meta), "--out", str(out),
                    "--report", str(report), "--overrides", str(overrides))
        assert r.returncode == 0, r.stderr
        fields = out.read_text().strip().split()
        assert fields[13] == "Germany"
        assert fields[14] == "match"

    def test_missing_meta_file_is_not_a_crash(self, tmp_path):
        out = tmp_path / "out"
        report = tmp_path / "report.tsv"
        r = run_cli("annotate", "--meta", str(tmp_path / "nope"),
                    "--out", str(out), "--report", str(report))
        assert r.returncode == 0
        assert out.read_text() == ""
        assert report.read_text().strip() != ""  # header only

    def test_latlon_vs_geo_summary_reported(self, tmp_path):
        meta = tmp_path / "all_asm_acc_metadata"
        meta.write_text(
            " ".join(row(geo="Japan", lat_lon="35.6895_N_139.6917_E")) + "\n")
        out = tmp_path / "out"
        report = tmp_path / "report.tsv"

        r = run_cli("annotate", "--meta", str(meta),
                    "--out", str(out), "--report", str(report))
        assert r.returncode == 0, r.stderr
        assert "latlon_vs_geo match: 1" in r.stderr

    def test_ror_flag_accepted_without_network_for_unusable_submitters(
            self, tmp_path):
        """--ror only triggers a network fetch for usable submitter text;
        a row with no usable submitter field never calls out, so this stays
        a pure offline test while still exercising the --ror/--ror-cache/
        --at-submitter-km CLI plumbing end to end."""
        meta = tmp_path / "all_asm_acc_metadata"
        meta.write_text(
            " ".join(row(geo="USA", lat_lon="38.98_N_77.11_W")) + "\n")
        out = tmp_path / "out"
        report = tmp_path / "report.tsv"
        cache = tmp_path / "ror_cache.json"

        r = run_cli("annotate", "--meta", str(meta), "--out", str(out),
                    "--report", str(report), "--ror", "--ror-cache", str(cache),
                    "--at-submitter-km", "10")
        assert r.returncode == 0, r.stderr
        fields = out.read_text().strip().split()
        assert len(fields) == 19
        assert fields[17] == "NA"          # submitter_km: nothing to look up
        assert fields[18] == "unresolved"  # latlon_vs_submitter
        assert "latlon_vs_submitter" in r.stderr
        assert cache.exists()
