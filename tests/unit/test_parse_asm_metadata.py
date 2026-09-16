"""Tests for python_scripts/parse_asm_metadata.py.

This module replaces three slow bash constructs in utils/gather_filter_asms.sh
(a nested O(n*m) `while read` join and two ~100k-pattern `grep -f` filters). The
tests therefore focus on *behavioural equivalence* with the bash they replace,
including the quirks that downstream column indices depend on:

  - merge-geoloc is an INNER join (unmatched metadata rows are dropped),
  - it emits one row per matching geoloc entry,
  - it collapses the metadata tabs to single spaces (the old unquoted echo did).

The transforms are pure functions over (line, fields) tuples, so they are tested
directly; the CLI layer is tested through tmp_path files.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "python_scripts" / "parse_asm_metadata.py"


@pytest.fixture
def mp(metadata_parser_module):
    return metadata_parser_module


def rows(*lines):
    """Build the (raw, fields) tuples the transforms consume."""
    return [(ln, ln.split()) for ln in lines]


# A realistic BS_to_meta row: 10 xtract columns.
#   1 BioSample 2 RefSeq 3 Genbank 4 Species 5 Sub_value
#   6 FtpGenBank 7 FtpRefSeq 8 Taxid 9 taxcheck 10 ExclFromRefSeq
META_A = "SAMN001\tGCF_000000001.1\tGCA_000000001.1\tEscherichia_coli\tK12\tftp://gb/A\tftp://rs/A\t562\tOK\tna"
META_B = "SAMN002\tGCF_000000002.1\tGCA_000000002.1\tEscherichia_coli\tW3110\tftp://gb/B\tftp://rs/B\t562\tInconclusive\tna"
META_C = "SAMN003\tNA\tGCA_000000003.1\tShigella_flexneri\t2a\tftp://gb/C\tNA\t623\tOK\tna"


# --------------------------------------------------------------------------- #
# merge_geoloc -- the hash join replacing the O(n*m) nested bash loop          #
# --------------------------------------------------------------------------- #

class TestMergeGeoloc:
    def test_joins_on_biosample(self, mp):
        out = list(mp.merge_geoloc(
            rows(META_A, META_B),
            {"SAMN001": ["USA"], "SAMN002": ["Japan"]},
        ))
        assert len(out) == 2
        assert out[0].startswith("SAMN001\t")
        assert out[0].endswith("\tUSA")
        assert out[1].endswith("\tJapan")

    def test_inner_join_drops_unmatched(self, mp):
        """A metadata row with no geoloc entry is dropped -- matches the bash."""
        out = list(mp.merge_geoloc(rows(META_A, META_B), {"SAMN001": ["USA"]}))
        assert len(out) == 1
        assert "SAMN001" in out[0]
        assert "SAMN002" not in out[0]

    def test_empty_geoloc_map_yields_nothing(self, mp):
        assert list(mp.merge_geoloc(rows(META_A, META_B), {})) == []

    def test_multiple_geoloc_matches_emit_multiple_rows(self, mp):
        """The bash loop emitted one line per matching geoloc row; so do we."""
        out = list(mp.merge_geoloc(rows(META_A), {"SAMN001": ["USA", "Canada"]}))
        assert len(out) == 2
        assert out[0].endswith("\tUSA")
        assert out[1].endswith("\tCanada")

    def test_collapses_metadata_whitespace(self, mp):
        """Metadata tabs collapse to single spaces (old unquoted `echo $BLAH`).

        Downstream awk splits on any whitespace, so the joined row must remain
        field-addressable: geoloc lands in $11 after the 10 metadata columns.
        """
        out = list(mp.merge_geoloc(rows(META_A), {"SAMN001": ["USA"]}))[0]
        fields = out.split()
        assert len(fields) == 11
        assert fields[0] == "SAMN001"
        assert fields[8] == "OK"       # taxonomy-check-status, awk $9
        assert fields[10] == "USA"     # geolocation, awk $11

    def test_geoloc_row_order_is_input_order(self, mp):
        """Determinism: output follows metadata order, not dict order."""
        out = list(mp.merge_geoloc(
            rows(META_C, META_A, META_B),
            {"SAMN001": ["USA"], "SAMN002": ["Japan"], "SAMN003": ["Egypt"]},
        ))
        assert [ln.split()[0] for ln in out] == ["SAMN003", "SAMN001", "SAMN002"]


class TestReadGeolocMap:
    def test_reads_tab_table(self, mp, tmp_path):
        f = tmp_path / "geo"
        f.write_text("SAMN001\tUSA\nSAMN002\tJapan\n")
        assert mp.read_geoloc_map(str(f)) == {
            "SAMN001": ["USA"], "SAMN002": ["Japan"],
        }

    def test_missing_file_is_empty(self, mp, tmp_path):
        assert mp.read_geoloc_map(str(tmp_path / "nope")) == {}

    def test_duplicate_accession_accumulates(self, mp, tmp_path):
        f = tmp_path / "geo"
        f.write_text("SAMN001\tUSA\nSAMN001\tCanada\n")
        assert mp.read_geoloc_map(str(f)) == {"SAMN001": ["USA", "Canada"]}

    def test_blank_lines_skipped(self, mp, tmp_path):
        f = tmp_path / "geo"
        f.write_text("SAMN001\tUSA\n\n   \nSAMN002\tJapan\n")
        assert len(mp.read_geoloc_map(str(f))) == 2

    def test_accession_with_no_location(self, mp, tmp_path):
        """A lone accession yields an empty location, not a crash."""
        f = tmp_path / "geo"
        f.write_text("SAMN001\n")
        assert mp.read_geoloc_map(str(f)) == {"SAMN001": [""]}


# --------------------------------------------------------------------------- #
# filter_by_acc / project_metadata -- replacing `grep -f` + awk                #
# --------------------------------------------------------------------------- #

# Merged rows (11 columns) as they appear in BS_to_all_meta.
MERGED_A = META_A.replace("\t", " ") + " USA"
MERGED_B = META_B.replace("\t", " ") + " Japan"
MERGED_C = META_C.replace("\t", " ") + " Egypt"


class TestProjectMetadata:
    def test_refseq_projection_matches_awk(self, mp):
        """awk '{print $2,$4,$4"."$5,$1,$8,$9,$10,$11}' for GCF rows."""
        fields = MERGED_A.split()
        assert mp.project_metadata(fields, acc_index=2) == " ".join([
            "GCF_000000001.1", "Escherichia_coli", "Escherichia_coli.K12",
            "SAMN001", "562", "OK", "na", "USA",
        ])

    def test_genbank_projection_uses_column_3(self, mp):
        fields = MERGED_C.split()
        assert mp.project_metadata(fields, acc_index=3).split()[0] == "GCA_000000003.1"

    def test_short_row_pads_like_awk(self, mp):
        """Fields past NF are empty strings, as in awk -- not an IndexError."""
        out = mp.project_metadata(["SAMN9", "GCF_9"], acc_index=2)
        assert out.split()[0] == "GCF_9"
        assert out.endswith(" ")  # trailing empty projected columns

    def test_taxcheck_lands_in_column_6(self, mp):
        """filter_asm_by_taxCheck reads $6 of the projected row."""
        out = mp.project_metadata(MERGED_B.split(), acc_index=2)
        assert out.split()[5] == "Inconclusive"

    def test_country_lands_in_column_8(self, mp):
        """annotate_tree_by_country.py reads column 8 of all_asm_acc_metadata."""
        out = mp.project_metadata(MERGED_A.split(), acc_index=2)
        assert out.split()[7] == "USA"


class TestFilterByAcc:
    def test_keeps_only_requested_refseq(self, mp):
        out = list(mp.filter_by_acc(
            rows(MERGED_A, MERGED_B), {"GCF_000000001.1"}, acc_index=2))
        assert len(out) == 1
        assert out[0].startswith("GCF_000000001.1 ")

    def test_genbank_column(self, mp):
        out = list(mp.filter_by_acc(
            rows(MERGED_C), {"GCA_000000003.1"}, acc_index=3))
        assert len(out) == 1

    def test_empty_acc_set_keeps_nothing(self, mp):
        assert list(mp.filter_by_acc(rows(MERGED_A, MERGED_B), set(), 2)) == []

    def test_exact_match_not_substring(self, mp):
        """Fixes a `grep -f` false positive: .1 must not match .10."""
        out = list(mp.filter_by_acc(
            rows(MERGED_A), {"GCF_000000001.10"}, acc_index=2))
        assert out == []

    def test_accession_in_ftp_path_does_not_match(self, mp):
        """Another `grep -f` false positive: match the column, not the line."""
        line = ("SAMN9 NA GCA_9.1 Sp st "
                "ftp://gb/GCF_000000001.1 NA 562 OK na USA")
        out = list(mp.filter_by_acc(rows(line), {"GCF_000000001.1"}, acc_index=2))
        assert out == []


class TestSplitAccessions:
    def test_splits_by_prefix(self, mp):
        gcf, gca = mp.split_accessions(
            {"GCF_1.1", "GCA_2.1", "GCF_3.1"})
        assert gcf == {"GCF_1.1", "GCF_3.1"}
        assert gca == {"GCA_2.1"}

    def test_ignores_other_prefixes(self, mp):
        gcf, gca = mp.split_accessions({"weird_name", "GCF_1.1"})
        assert gcf == {"GCF_1.1"}
        assert gca == set()


# --------------------------------------------------------------------------- #
# exclude_downloaded -- replacing `grep -vf`                                   #
# --------------------------------------------------------------------------- #

class TestExcludeDownloaded:
    def test_drops_rows_matching_refseq(self, mp):
        out = list(mp.exclude_downloaded(
            rows(META_A, META_B), {"GCF_000000001.1"}))
        assert len(out) == 1
        assert "SAMN002" in out[0]

    def test_drops_rows_matching_genbank(self, mp):
        """A GenBank-only assembly is excluded via column 3."""
        out = list(mp.exclude_downloaded(rows(META_C), {"GCA_000000003.1"}))
        assert out == []

    def test_keeps_everything_when_nothing_downloaded(self, mp):
        out = list(mp.exclude_downloaded(rows(META_A, META_B, META_C), set()))
        assert len(out) == 3

    def test_passes_lines_through_verbatim(self, mp):
        """The caller re-parses these with its own field layout, so preserve tabs."""
        out = list(mp.exclude_downloaded(rows(META_A), set()))
        assert out[0] == META_A


# --------------------------------------------------------------------------- #
# CLI -- what gather_filter_asms.sh actually invokes                           #
# --------------------------------------------------------------------------- #

def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + list(args),
        capture_output=True, text=True,
    )


class TestCLI:
    def test_help(self):
        assert run_cli("--help").returncode == 0

    def test_merge_geoloc_end_to_end(self, tmp_path):
        meta = tmp_path / "meta"
        meta.write_text(META_A + "\n" + META_B + "\n")
        geo = tmp_path / "geo"
        geo.write_text("SAMN001\tUSA\nSAMN002\tJapan\n")
        out = tmp_path / "merged"

        r = run_cli("merge-geoloc", "--meta", str(meta),
                    "--geoloc", str(geo), "--out", str(out))
        assert r.returncode == 0, r.stderr
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2
        assert all(len(ln.split()) == 11 for ln in lines)

    def test_merge_geoloc_empty_geoloc_warns_and_succeeds(self, tmp_path):
        """Geolocation is decorative; a missing table must not fail the run."""
        meta = tmp_path / "meta"
        meta.write_text(META_A + "\n")
        geo = tmp_path / "geo"
        geo.write_text("")
        out = tmp_path / "merged"

        r = run_cli("merge-geoloc", "--meta", str(meta),
                    "--geoloc", str(geo), "--out", str(out))
        assert r.returncode == 0
        assert out.read_text() == ""
        assert "WARNING" in r.stderr

    def test_filter_by_acc_end_to_end(self, tmp_path):
        merged = tmp_path / "all_meta"
        merged.write_text("\n".join([MERGED_A, MERGED_B, MERGED_C]) + "\n")
        accs = tmp_path / "all_asm_acc"
        accs.write_text("GCF_000000001.1\nGCA_000000003.1\n")
        out = tmp_path / "out"
        gcf = tmp_path / "out.GCF"
        gca = tmp_path / "out.GCA"

        r = run_cli("filter-by-acc", "--meta", str(merged),
                    "--acc-list", str(accs), "--out", str(out),
                    "--out-gcf", str(gcf), "--out-gca", str(gca))
        assert r.returncode == 0, r.stderr

        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2
        # GCF rows first, then GCA -- same layout the bash appended in.
        assert lines[0].split()[0] == "GCF_000000001.1"
        assert lines[1].split()[0] == "GCA_000000003.1"
        # Intermediate accession lists preserved for inspection.
        assert gcf.read_text().strip() == "GCF_000000001.1"
        assert gca.read_text().strip() == "GCA_000000003.1"

    def test_exclude_downloaded_end_to_end(self, tmp_path):
        meta = tmp_path / "meta"
        meta.write_text(META_A + "\n" + META_B + "\n")
        dl = tmp_path / "downloaded"
        dl.write_text("GCF_000000001.1\n")
        out = tmp_path / "out"

        r = run_cli("exclude-downloaded", "--meta", str(meta),
                    "--downloaded", str(dl), "--out", str(out))
        assert r.returncode == 0, r.stderr
        assert out.read_text().strip() == META_B

    def test_missing_input_file_is_not_a_crash(self, tmp_path):
        """gather_filter_asms.sh can legitimately pass an empty/absent table."""
        out = tmp_path / "out"
        r = run_cli("merge-geoloc", "--meta", str(tmp_path / "nope"),
                    "--geoloc", str(tmp_path / "nope2"), "--out", str(out))
        assert r.returncode == 0
        assert out.read_text() == ""
