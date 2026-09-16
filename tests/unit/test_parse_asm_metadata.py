"""Tests for python_scripts/parse_asm_metadata.py.

This module replaces slow and/or incorrect bash in utils/gather_filter_asms.sh:
a nested O(n*m) `while read` join, two ~100k-pattern `grep -f` filters, and two
`xtract` pipelines that shifted columns. The tests therefore focus on
*behavioural equivalence* with the bash they replace, including the quirks that
downstream column indices depend on:

  - merge-geoloc is an INNER join (unmatched metadata rows are dropped),
  - it emits one row per matching geoloc entry,
  - it collapses the metadata tabs to single spaces (the old unquoted echo did).

The XML tests pin the column-shift bug that motivated parse-assembly-xml: a
repeated <Sub_value> element made xtract emit extra tokens, pushing
taxonomy-check-status into the geolocation column ("geolocation == OK"). They
also cover the multi-root layout of real esearch output.

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
# XML parsing -- the column-shift fix                                          #
# --------------------------------------------------------------------------- #

def asm_record(biosample="SAMN001", refseq="GCF_1.1", genbank="GCA_1.1",
               species="Pseudomonas aeruginosa", sub_values=("K12",),
               taxid="287", taxcheck="OK", excl=None, ftp=True):
    """Build one assembly <DocumentSummary>, mirroring real esummary layout.

    sub_values may hold several entries -- that is the bug trigger: the real XML
    nests one <Infraspecie> per entry (e.g. culture-collection + strain).
    """
    parts = ["  <DocumentSummary>"]
    parts.append("    <BioSampleAccn>%s</BioSampleAccn>" % biosample)
    parts.append("    <Synonym>")
    if refseq:
        parts.append("      <RefSeq>%s</RefSeq>" % refseq)
    parts.append("      <Genbank>%s</Genbank>" % genbank)
    parts.append("    </Synonym>")
    parts.append("    <SpeciesName>%s</SpeciesName>" % species)
    if sub_values:
        parts.append("    <Biosource>")
        parts.append("      <InfraspeciesList>")
        for sv in sub_values:
            parts.append("        <Infraspecie>")
            parts.append("          <Sub_type>strain</Sub_type>")
            parts.append("          <Sub_value>%s</Sub_value>" % sv)
            parts.append("        </Infraspecie>")
        parts.append("      </InfraspeciesList>")
        parts.append("    </Biosource>")
    if ftp:
        parts.append("    <FtpPath_GenBank>ftp://gb/%s</FtpPath_GenBank>" % genbank)
        if refseq:
            parts.append("    <FtpPath_RefSeq>ftp://rs/%s</FtpPath_RefSeq>" % refseq)
    parts.append("    <Taxid>%s</Taxid>" % taxid)
    if excl:
        parts.append("    <ExclFromRefSeq>%s</ExclFromRefSeq>" % excl)
    parts.append("    <Meta>")
    if taxcheck:
        parts.append("      <taxonomy-check-status>%s</taxonomy-check-status>"
                     % taxcheck)
    parts.append("    </Meta>")
    parts.append("  </DocumentSummary>")
    return "\n".join(parts)


def asm_xml(*records, **kw):
    """Wrap records in <DocumentSummarySet> root(s).

    n_roots>1 reproduces the real concatenated-batch layout: `esearch | esummary`
    appends one root PER BATCH, so a 75k-assembly file has 76 roots and is not a
    single XML document.
    """
    n_roots = kw.get("n_roots", 1)
    chunks = []
    per_root = max(1, len(records) // n_roots) if n_roots > 1 else len(records)
    groups = [records[i:i + per_root] for i in range(0, len(records), per_root)] \
        if records else [()]
    for group in groups:
        chunks.append('<?xml version="1.0" encoding="UTF-8" ?>')
        chunks.append("<!DOCTYPE DocumentSummarySet>")
        chunks.append('<DocumentSummarySet status="OK">')
        chunks.append("  <DbBuild>Build260914-2040.1</DbBuild>")
        chunks.extend(group)
        chunks.append("</DocumentSummarySet>")
    return "\n".join(chunks) + "\n"


def bios_record(accession="SAMN001", geoloc="USA", other_attrs=True):
    parts = ["  <DocumentSummary>"]
    parts.append("    <Accession>%s</Accession>" % accession)
    parts.append("    <SampleData>")
    parts.append("      <BioSample accession=\"%s\">" % accession)
    parts.append("        <Attributes>")
    if other_attrs:
        parts.append('          <Attribute attribute_name="strain" '
                     'harmonized_name="strain">someStrain</Attribute>')
    if geoloc is not None:
        parts.append('          <Attribute attribute_name="geo_loc_name" '
                     'harmonized_name="geo_loc_name">%s</Attribute>' % geoloc)
    parts.append("        </Attributes>")
    parts.append("      </BioSample>")
    parts.append("    </SampleData>")
    parts.append("  </DocumentSummary>")
    return "\n".join(parts)


class TestAssemblyRow:
    def test_ten_columns_for_simple_record(self, mp):
        import xml.etree.ElementTree as ET
        row = mp.assembly_row(ET.fromstring(asm_record()))
        assert len(row) == mp.EXPECTED_META_COLUMNS

    def test_multi_subvalue_does_not_shift_columns(self, mp):
        """THE BUG. Extra <Sub_value> entries must not move later columns."""
        import xml.etree.ElementTree as ET
        for n in (1, 2, 3, 4):
            rec = asm_record(sub_values=tuple("s%d" % i for i in range(n)))
            row = mp.assembly_row(ET.fromstring(rec))
            assert len(row) == 10, "n=%d gave %d columns" % (n, len(row))
            # taxcheck stays in column 9 and taxid in column 8 (1-based).
            assert row[7] == "287"
            assert row[8] == "OK"

    def test_only_first_subvalue_is_kept(self, mp):
        import xml.etree.ElementTree as ET
        row = mp.assembly_row(ET.fromstring(
            asm_record(sub_values=("ATCC:19660", "Xen5"))))
        assert row[4] == "ATCC:19660"

    def test_missing_fields_become_NA(self, mp):
        import xml.etree.ElementTree as ET
        row = mp.assembly_row(ET.fromstring(
            asm_record(refseq=None, sub_values=(), taxcheck=None, ftp=False)))
        assert len(row) == 10
        assert row[1] == "NA"   # RefSeq
        assert row[4] == "NA"   # Sub_value
        assert row[8] == "NA"   # taxcheck

    def test_species_spaces_become_underscores(self, mp):
        """Otherwise a two-word species name would split into two columns."""
        import xml.etree.ElementTree as ET
        row = mp.assembly_row(ET.fromstring(asm_record()))
        assert row[3] == "Pseudomonas_aeruginosa"

    def test_row_never_contains_whitespace(self, mp):
        import xml.etree.ElementTree as ET
        row = mp.assembly_row(ET.fromstring(
            asm_record(species="Genus species subsp thing",
                       sub_values=("has space",))))
        for field in row:
            assert " " not in field and "\t" not in field


class TestNormalizeGeoloc:
    @pytest.mark.parametrize("raw,expected", [
        ("USA", "USA"),
        ("India: Chhatrapati Sambhajinagar", "India"),   # 75,609 real rows
        ("United Kingdom", "United_Kingdom"),
        ("USA: CA: San Diego", "USA"),
        ("missing", "Unknown"),
        ("Missing: control sample", "Unknown"),
        ("unknown", "Unknown"),
        ("not collected", "Unknown"),
        ("not_collected", "Unknown"),
        ("not applicable", "Unknown"),
        ("NONE", "Unknown"),
        ("not_provided", "Unknown"),   # the 296 rows the old seds missed
        ("", "NA"),
        (None, "NA"),
    ])
    def test_normalization(self, mp, raw, expected):
        assert mp.normalize_geoloc(raw) == expected

    def test_case_insensitive_unknown(self, mp):
        for spelling in ("MISSING", "Unknown", "uNkNoWn", "Not Collected"):
            assert mp.normalize_geoloc(spelling) == "Unknown"


class TestBiosampleGeoloc:
    def test_extracts_geoloc_attribute(self, mp):
        import xml.etree.ElementTree as ET
        acc, geo = mp.biosample_geoloc(ET.fromstring(bios_record()))
        assert (acc, geo) == ("SAMN001", "USA")

    def test_picks_geoloc_not_other_attributes(self, mp):
        """Must select by harmonized_name, not by position."""
        import xml.etree.ElementTree as ET
        acc, geo = mp.biosample_geoloc(
            ET.fromstring(bios_record(geoloc="China", other_attrs=True)))
        assert geo == "China"

    def test_absent_geoloc_is_NA(self, mp):
        import xml.etree.ElementTree as ET
        acc, geo = mp.biosample_geoloc(ET.fromstring(bios_record(geoloc=None)))
        assert (acc, geo) == ("SAMN001", "NA")


class TestIterDocumentSummaries:
    def test_reads_single_root(self, mp, tmp_path):
        f = tmp_path / "a.xml"
        f.write_text(asm_xml(asm_record("SAMN001"), asm_record("SAMN002")))
        assert len(list(mp.iter_document_summaries(str(f)))) == 2

    def test_reads_concatenated_roots(self, mp, tmp_path):
        """Real esearch output is many <DocumentSummarySet> roots concatenated.

        A strict whole-document parser stops after the first root -- that bug
        silently returned 1,000 of 75,155 records on the real Pseudomonas file.
        """
        f = tmp_path / "a.xml"
        recs = [asm_record("SAMN%03d" % i) for i in range(10)]
        f.write_text(asm_xml(*recs, n_roots=5))
        got = list(mp.iter_document_summaries(str(f)))
        assert len(got) == 10

    def test_missing_file_yields_nothing(self, mp, tmp_path):
        assert list(mp.iter_document_summaries(str(tmp_path / "nope"))) == []

    def test_truncated_tail_keeps_earlier_records(self, mp, tmp_path):
        """A dropped connection truncates mid-record; keep what parsed."""
        f = tmp_path / "a.xml"
        text = asm_xml(asm_record("SAMN001"), asm_record("SAMN002"))
        f.write_text(text[:len(text) - 40])
        assert len(list(mp.iter_document_summaries(str(f)))) >= 1

    def test_memory_flat_over_many_records(self, mp, tmp_path):
        """Generator must not accumulate: only one record held at a time."""
        f = tmp_path / "a.xml"
        f.write_text(asm_xml(*[asm_record("SAMN%04d" % i) for i in range(500)]))
        n = sum(1 for _ in mp.iter_document_summaries(str(f)))
        assert n == 500


class TestXMLCLI:
    def test_parse_assembly_xml_end_to_end(self, tmp_path):
        xml = tmp_path / "asm.xml"
        xml.write_text(asm_xml(
            asm_record("SAMN001", sub_values=("K12",)),
            asm_record("SAMN002", sub_values=("ATCC:19660", "Xen5")),
            asm_record("SAMN003", sub_values=("a", "b", "c")),
            n_roots=2,
        ))
        out = tmp_path / "BS_to_meta"
        r = run_cli("parse-assembly-xml", "--xml", str(xml), "--out", str(out))
        assert r.returncode == 0, r.stderr

        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3
        # Every row has exactly 10 columns regardless of Sub_value count.
        for ln in lines:
            assert len(ln.split("\t")) == 10
            assert len(ln.split()) == 10
        # The multi-Sub_value records are reported so the log shows the impact.
        assert "2 record(s) had multiple <Sub_value>" in r.stderr

    def test_parsed_output_feeds_merge_geoloc_correctly(self, tmp_path):
        """End-to-end: the fixed columns must survive the join.

        This is the actual bug's blast radius -- a shifted row put "OK" where the
        geolocation belongs. Chain both stages and assert the geolocation column
        really holds the country.
        """
        xml = tmp_path / "asm.xml"
        xml.write_text(asm_xml(
            asm_record("SAMN001", sub_values=("a", "b", "c"), taxcheck="OK")))
        meta = tmp_path / "BS_to_meta"
        assert run_cli("parse-assembly-xml", "--xml", str(xml),
                       "--out", str(meta)).returncode == 0

        geo = tmp_path / "geo"
        geo.write_text("SAMN001\tSweden\n")
        merged = tmp_path / "merged"
        r = run_cli("merge-geoloc", "--meta", str(meta),
                    "--geoloc", str(geo), "--out", str(merged))
        assert r.returncode == 0, r.stderr

        fields = merged.read_text().strip().split()
        assert len(fields) == 11
        assert fields[8] == "OK"       # taxcheck stayed in column 9
        assert fields[10] == "Sweden"  # geolocation is a country, not "OK"
        assert "WARNING" not in r.stderr

    def test_merge_geoloc_warns_on_shifted_legacy_input(self, tmp_path):
        """Recurrence guard: an 11-column legacy row must be flagged."""
        meta = tmp_path / "BS_to_meta"
        meta.write_text(
            "SAMN001\tGCF_1.1\tGCA_1.1\tSp\tsv1\tEXTRA\tftp://gb\tftp://rs"
            "\t287\tOK\tna\n")
        geo = tmp_path / "geo"
        geo.write_text("SAMN001\tSweden\n")
        r = run_cli("merge-geoloc", "--meta", str(meta),
                    "--geoloc", str(geo), "--out", str(tmp_path / "m"))
        assert r.returncode == 0
        assert "did not have the expected 10 columns" in r.stderr

    def test_parse_biosample_xml_end_to_end(self, tmp_path):
        xml = tmp_path / "bios.xml"
        # Reuse the biosample root wrapper via asm_xml (same envelope).
        xml.write_text(asm_xml(
            bios_record("SAMN001", "USA"),
            bios_record("SAMN002", "India: Pune"),
            bios_record("SAMN003", "missing"),
            bios_record("SAMN004", None),
        ))
        out = tmp_path / "geo"
        r = run_cli("parse-biosample-xml", "--xml", str(xml), "--out", str(out))
        assert r.returncode == 0, r.stderr

        got = dict(ln.split("\t") for ln in out.read_text().strip().split("\n"))
        assert got == {
            "SAMN001": "USA",
            "SAMN002": "India",
            "SAMN003": "Unknown",
            "SAMN004": "NA",
        }

    def test_empty_xml_warns_and_succeeds(self, tmp_path):
        xml = tmp_path / "empty.xml"
        xml.write_text("")
        out = tmp_path / "out"
        r = run_cli("parse-assembly-xml", "--xml", str(xml), "--out", str(out))
        assert r.returncode == 0
        assert out.read_text() == ""
        assert "WARNING" in r.stderr


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
