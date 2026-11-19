"""Tests for python_scripts/geo_provenance.py.

geo_provenance.py flags cases where all_asm_acc_metadata's col-8 geo_loc_name
country might actually be the sequencing/submitting lab's own country, by
comparing it against a country resolved from sequenced_by / asm_submitter /
biosample_owner. It is flag-only: these tests pin the exact flag matrix
(match/differ/no_geo/unresolved) and the country-resolution heuristics
(override table, free-text scanning with its ambiguity guards), not just that
"something" gets written.
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


# A 12-column all_asm_acc_metadata row, matching project_metadata()'s layout:
#   1 acc 2 species 3 species.strain 4 biosample 5 taxid 6 taxcheck 7 excl
#   8 geo  9 collected_by 10 sequenced_by 11 biosample_owner 12 asm_submitter
def row(geo="NA", collected_by="NA", sequenced_by="NA", biosample_owner="NA",
        asm_submitter="NA"):
    return [
        "GCF_1.1", "Esch_coli", "Esch_coli.K12", "SAMN001", "562", "OK", "na",
        geo, collected_by, sequenced_by, biosample_owner, asm_submitter,
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
# annotate_rows                                                              #
# --------------------------------------------------------------------------- #

class TestAnnotateRows:
    def test_match_row(self, gp, index, overrides):
        r = row(geo="USA", sequenced_by="PulseNet")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert len(extended) == 14
        assert extended[:12] == r
        assert extended[12] == "United_States_of_America"
        assert extended[13] == "match"
        assert rec["submitter_source"] == "sequenced_by"

    def test_differ_row(self, gp, index, overrides):
        r = row(geo="Japan", sequenced_by="PulseNet")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert extended[12] == "United_States_of_America"
        assert extended[13] == "differ"

    def test_submitter_country_is_sanitized_for_the_whitespace_split_output(
            self, gp, index, overrides):
        """A multi-word canonical country name must not create extra columns
        in the whitespace-split output row."""
        r = row(geo="USA", sequenced_by="PulseNet")
        extended, _ = next(gp.annotate_rows([r], index, overrides))
        assert " " not in extended[12]
        assert len(" ".join(extended).split()) == 14

    def test_no_geo_row(self, gp, index, overrides):
        r = row(geo="NA", sequenced_by="PulseNet")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert extended[13] == "no_geo"

    def test_unresolved_row(self, gp, index, overrides):
        r = row(geo="USA")
        extended, rec = next(gp.annotate_rows([r], index, overrides))
        assert extended[12] == "NA"
        assert extended[13] == "unresolved"

    def test_short_row_is_padded_with_NA(self, gp, index, overrides):
        """A row missing columns 9-12 (e.g. an older all_asm_acc_metadata) is
        padded rather than raising."""
        short = ["GCF_1.1", "Esch_coli", "Esch_coli.K12", "SAMN001",
                 "562", "OK", "na", "USA"]
        extended, rec = next(gp.annotate_rows([short], index, overrides))
        assert len(extended) == 14
        assert rec["collected_by"] == "NA"
        assert rec["sequenced_by"] == "NA"
        assert rec["biosample_owner"] == "NA"
        assert rec["asm_submitter"] == "NA"
        assert extended[13] == "unresolved"

    def test_multiple_rows_preserve_order(self, gp, index, overrides):
        rows_in = [row(geo="USA", sequenced_by="PulseNet"),
                   row(geo="Japan", sequenced_by="PulseNet")]
        out = list(gp.annotate_rows(rows_in, index, overrides))
        assert len(out) == 2
        assert out[0][0][13] == "match"
        assert out[1][0][13] == "differ"


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
        assert "14 columns" in r.stderr
        assert "match: 1" in r.stderr
        assert "differ: 1" in r.stderr

        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3
        for ln in lines:
            assert len(ln.split()) == 14

        report_lines = report.read_text().strip().split("\n")
        assert len(report_lines) == 4  # header + 3 rows
        header = report_lines[0].split("\t")
        assert header[-1] == "geo_vs_submitter"
        assert header[-2] == "submitter_source"
        assert header[-3] == "submitter_country"

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
        assert fields[12] == "Germany"
        assert fields[13] == "match"

    def test_missing_meta_file_is_not_a_crash(self, tmp_path):
        out = tmp_path / "out"
        report = tmp_path / "report.tsv"
        r = run_cli("annotate", "--meta", str(tmp_path / "nope"),
                    "--out", str(out), "--report", str(report))
        assert r.returncode == 0
        assert out.read_text() == ""
        assert report.read_text().strip() != ""  # header only
