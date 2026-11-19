"""Tests for utils/gather_viral_annots.py (NCBI viral annotation download,
Task 2 of the viral-support plan).

Fixtures are deliberately small/synthetic but mirror real `datasets download
virus genome taxon ... --include cds,protein,annotation` output confirmed by
hand against real packages (Norwalk virus, Lassa mammarenavirus) -- see the
module docstring for the exact format facts. In particular:
  - the ambisense/minus-strand case (Lassa S segment: GPC forward, NP only on
    the reverse complement, no strand field anywhere in the metadata), and
  - the ribosomal-frameshift frame-offset case (Norovirus ORF1/ORF2 junction:
    the recorded CDS range starts 2 nt before the true translation start)
are both reproduced here with minimal synthetic sequences, since downloading
real packages in the test suite would require live network access.
"""

import json

import pytest

from conftest import load_module, REPO_ROOT


@pytest.fixture(scope="module")
def gva():
    return load_module(REPO_ROOT / "utils" / "gather_viral_annots.py",
                        name="gather_viral_annots_test")


# --------------------------------------------------------------------------- #
# parse_jsonlines / load_fasta                                                #
# --------------------------------------------------------------------------- #

class TestParseJsonlines:
    def test_parses_each_line(self, gva):
        text = '{"a": 1}\n{"a": 2}\n'
        assert gva.parse_jsonlines(text) == [{"a": 1}, {"a": 2}]

    def test_skips_blank_lines(self, gva):
        text = '{"a": 1}\n\n\n{"a": 2}\n'
        assert gva.parse_jsonlines(text) == [{"a": 1}, {"a": 2}]

    def test_skips_unparseable_lines(self, gva):
        text = '{"a": 1}\nnot json\n{"a": 2}\n'
        assert gva.parse_jsonlines(text) == [{"a": 1}, {"a": 2}]

    def test_empty_text(self, gva):
        assert gva.parse_jsonlines("") == []


class TestLoadFasta:
    def test_basic_single_line(self, gva, tmp_path):
        f = tmp_path / "x.fna"
        f.write_text(">PV275052.1:5-5104 some title\nACGTACGT\n")
        recs = gva.load_fasta(f)
        assert recs["PV275052.1:5-5104"] == (
            "PV275052.1:5-5104 some title", "ACGTACGT")

    def test_multiline_sequence_is_joined(self, gva, tmp_path):
        f = tmp_path / "x.fna"
        f.write_text(">id1\nACGT\nACGT\n>id2\nTTTT\n")
        recs = gva.load_fasta(f)
        assert recs["id1"] == ("id1", "ACGTACGT")
        assert recs["id2"] == ("id2", "TTTT")

    def test_seqid_is_first_whitespace_token(self, gva, tmp_path):
        f = tmp_path / "x.faa"
        f.write_text(">BHY91758.1:1-1698 nonstructural polyprotein [x]\nMAAA\n")
        recs = gva.load_fasta(f)
        assert "BHY91758.1:1-1698" in recs


# --------------------------------------------------------------------------- #
# genome_group_key / sanitize_genome_id / dedupe_prefer_refseq               #
# --------------------------------------------------------------------------- #

class TestGenomeGroupKey:
    def test_uses_isolate_name_when_present(self, gva):
        rec = {"accession": "LC936410.1",
               "isolate": {"name": "Hu/GII.17[P17]/TKY250177/2025/JP"}}
        assert gva.genome_group_key(rec) == "Hu/GII.17[P17]/TKY250177/2025/JP"

    def test_falls_back_to_accession_when_isolate_absent(self, gva):
        rec = {"accession": "PV275052.1"}
        assert gva.genome_group_key(rec) == "PV275052.1"

    def test_falls_back_to_accession_when_isolate_name_blank(self, gva):
        rec = {"accession": "PV275052.1", "isolate": {"name": "  "}}
        assert gva.genome_group_key(rec) == "PV275052.1"

    def test_shared_isolate_groups_segments(self, gva):
        s = {"accession": "LASSA_S", "isolate": {"name": "LF2350"}}
        l = {"accession": "LASSA_L", "isolate": {"name": "LF2350"}}
        assert gva.genome_group_key(s) == gva.genome_group_key(l)


class TestSanitizeGenomeId:
    def test_leaves_simple_names_unchanged(self, gva):
        assert gva.sanitize_genome_id("PV275052.1") == "PV275052.1"

    def test_replaces_unsafe_characters(self, gva):
        result = gva.sanitize_genome_id("Hu/GII.17[P17]/TKY250177/2025/JP")
        assert "/" not in result and "[" not in result and "]" not in result

    def test_empty_string_falls_back(self, gva):
        assert gva.sanitize_genome_id("") == "unnamed"

    def test_strips_leading_trailing_dots_underscores(self, gva):
        assert gva.sanitize_genome_id("..weird__") == "weird"


class TestDedupePreferRefseq:
    def test_keeps_refseq_over_genbank_duplicate(self, gva):
        genbank = {"accession": "MK000001.1", "sourceDatabase": "GenBank",
                   "nucleotide": {"sequenceHash": "H1"}}
        refseq = {"accession": "NC_000001.1", "sourceDatabase": "RefSeq",
                  "nucleotide": {"sequenceHash": "H1"}}
        kept = gva.dedupe_prefer_refseq([genbank, refseq])
        assert len(kept) == 1
        assert kept[0]["accession"] == "NC_000001.1"

    def test_order_independent_refseq_still_wins(self, gva):
        genbank = {"accession": "MK000001.1", "sourceDatabase": "GenBank",
                   "nucleotide": {"sequenceHash": "H1"}}
        refseq = {"accession": "NC_000001.1", "sourceDatabase": "RefSeq",
                  "nucleotide": {"sequenceHash": "H1"}}
        kept = gva.dedupe_prefer_refseq([refseq, genbank])
        assert len(kept) == 1
        assert kept[0]["accession"] == "NC_000001.1"

    def test_no_hash_always_kept(self, gva):
        a = {"accession": "A", "nucleotide": {}}
        b = {"accession": "B", "nucleotide": {}}
        assert gva.dedupe_prefer_refseq([a, b]) == [a, b]

    def test_distinct_hashes_both_kept(self, gva):
        a = {"accession": "A", "nucleotide": {"sequenceHash": "H1"}}
        b = {"accession": "B", "nucleotide": {"sequenceHash": "H2"}}
        assert gva.dedupe_prefer_refseq([a, b]) == [a, b]


# --------------------------------------------------------------------------- #
# extract_cds_entries / build_genomes / filter_genomes_by_gene_count         #
# --------------------------------------------------------------------------- #

def _noro_annotation_rec():
    """Mirrors the real annotation_report.jsonl record shape (trimmed),
    including the maturePeptide-bearing polyprotein case."""
    return {
        "accession": "LC936410.1",
        "genes": [
            {
                "name": "CDS: ",
                "cds": [{
                    "name": "nonstructural polyprotein",
                    "nucleotide": {"seqId": "LC936410.1:1-5099"},
                    "protein": {"seqId": "BHY91758.1:1-1698"},
                    "maturePeptide": [
                        {"name": "p48",
                         "nucleotide": {"seqId": "LC936410.1:1-1000"},
                         "protein": {"seqId": "BHY91758.1:1-333"}},
                    ],
                }],
            },
            {
                "name": "CDS: ",
                "cds": [{
                    "name": "capsid protein VP1",
                    "nucleotide": {"seqId": "LC936410.1:5080-6696"},
                    "protein": {"seqId": "BHY91759.1:1-538"},
                }],
            },
        ],
    }


class TestExtractCdsEntries:
    def test_yields_one_entry_per_top_level_cds(self, gva):
        entries = list(gva.extract_cds_entries(_noro_annotation_rec()))
        assert len(entries) == 2
        nt_ids = {e["nt_seqid"] for e in entries}
        assert nt_ids == {"LC936410.1:1-5099", "LC936410.1:5080-6696"}

    def test_mature_peptides_are_not_extracted(self, gva):
        entries = list(gva.extract_cds_entries(_noro_annotation_rec()))
        nt_ids = {e["nt_seqid"] for e in entries}
        assert "LC936410.1:1-1000" not in nt_ids

    def test_cds_missing_seqid_is_skipped(self, gva):
        rec = {"accession": "A", "genes": [
            {"name": "CDS: ", "cds": [{"nucleotide": {}, "protein": {}}]},
        ]}
        assert list(gva.extract_cds_entries(rec)) == []

    def test_no_genes_yields_nothing(self, gva):
        assert list(gva.extract_cds_entries({"accession": "A"})) == []


class TestBuildGenomes:
    def test_single_accession_genome(self, gva):
        data_report = [{"accession": "PV275052.1"}]
        genomes = gva.build_genomes(data_report, [_noro_annotation_rec_for("PV275052.1")])
        assert len(genomes) == 1
        genome = next(iter(genomes.values()))
        assert genome["representative_accession"] == "PV275052.1"
        assert len(genome["entries"]) == 2

    def test_segmented_genome_groups_by_isolate(self, gva):
        data_report = [
            {"accession": "LASSA_S", "isolate": {"name": "LF2350"}},
            {"accession": "LASSA_L", "isolate": {"name": "LF2350"}},
        ]
        annotation = [
            {"accession": "LASSA_S", "genes": [
                {"name": "CDS: ", "cds": [{
                    "nucleotide": {"seqId": "LASSA_S:1-100"},
                    "protein": {"seqId": "P_S:1-30"}}]}]},
            {"accession": "LASSA_L", "genes": [
                {"name": "CDS: ", "cds": [{
                    "nucleotide": {"seqId": "LASSA_L:1-200"},
                    "protein": {"seqId": "P_L:1-60"}}]}]},
        ]
        genomes = gva.build_genomes(data_report, annotation)
        assert len(genomes) == 1
        genome = next(iter(genomes.values()))
        assert set(genome["accessions"]) == {"LASSA_S", "LASSA_L"}
        assert len(genome["entries"]) == 2

    def test_accession_absent_from_annotation_report_is_dropped(self, gva):
        data_report = [{"accession": "NO_ANNOT"}]
        genomes = gva.build_genomes(data_report, [])
        assert genomes == {}

    def test_duplicate_isolate_names_get_unique_ids(self, gva):
        data_report = [
            {"accession": "A1", "isolate": {"name": "Strain1"}},
            {"accession": "A2", "isolate": {"name": "Strain1"}},
        ]
        ann = [
            {"accession": "A1", "genes": [{"name": "CDS: ", "cds": [{
                "nucleotide": {"seqId": "A1:1-9"}, "protein": {"seqId": "P1:1-3"}}]}]},
        ]
        # A2 has no annotation entry -> only A1's group (keyed "Strain1")
        # should materialize; this also exercises that a shared raw_key with
        # a partially-missing accession still produces exactly one genome.
        genomes = gva.build_genomes(data_report, ann)
        assert len(genomes) == 1

    def test_distinct_isolate_name_collision_gets_suffixed(self, gva):
        # Two *different* raw keys that sanitize to the same filename stem.
        data_report = [
            {"accession": "A1", "isolate": {"name": "Strain/1"}},
            {"accession": "A2", "isolate": {"name": "Strain_1"}},
        ]
        ann = [
            {"accession": "A1", "genes": [{"name": "CDS: ", "cds": [{
                "nucleotide": {"seqId": "A1:1-9"}, "protein": {"seqId": "P1:1-3"}}]}]},
            {"accession": "A2", "genes": [{"name": "CDS: ", "cds": [{
                "nucleotide": {"seqId": "A2:1-9"}, "protein": {"seqId": "P2:1-3"}}]}]},
        ]
        genomes = gva.build_genomes(data_report, ann)
        assert len(genomes) == 2
        assert sorted(genomes.keys()) == ["Strain_1", "Strain_1_2"]


def _noro_annotation_rec_for(accession):
    rec = _noro_annotation_rec()
    rec["accession"] = accession
    rec["genes"][0]["cds"][0]["nucleotide"]["seqId"] = f"{accession}:1-5099"
    rec["genes"][1]["cds"][0]["nucleotide"]["seqId"] = f"{accession}:5080-6696"
    return rec


class TestFilterGenomesByGeneCount:
    def test_drops_low_count_genome(self, gva):
        genomes = {
            "g1": {"entries": [1, 2, 3, 4]},
            "g2": {"entries": [1, 2, 3, 4]},
            "g3": {"entries": [1]},
        }
        kept, dropped = gva.filter_genomes_by_gene_count(genomes, min_frac=0.5)
        assert "g3" in dropped
        assert set(kept.keys()) == {"g1", "g2"}

    def test_empty_input(self, gva):
        kept, dropped = gva.filter_genomes_by_gene_count({})
        assert kept == {} and dropped == []

    def test_uniform_counts_keeps_all(self, gva):
        genomes = {"g1": {"entries": [1, 2]}, "g2": {"entries": [1, 2]}}
        kept, dropped = gva.filter_genomes_by_gene_count(genomes)
        assert set(kept.keys()) == {"g1", "g2"}
        assert dropped == []


# --------------------------------------------------------------------------- #
# has_internal_stop / ambiguous_base_fraction                                #
# --------------------------------------------------------------------------- #

class TestHasInternalStop:
    def test_no_stop_until_the_end(self, gva):
        assert gva.has_internal_stop("ATGAAATAA") is False

    def test_stop_before_the_end_is_flagged(self, gva):
        assert gva.has_internal_stop("ATGTAAAAA") is True

    def test_empty_sequence(self, gva):
        assert gva.has_internal_stop("") is False


class TestAmbiguousBaseFraction:
    def test_pure_acgt_is_zero(self, gva):
        assert gva.ambiguous_base_fraction("ACGTACGT") == 0.0

    def test_counts_non_acgt(self, gva):
        assert gva.ambiguous_base_fraction("ACGTNNNN") == 0.5

    def test_empty_sequence(self, gva):
        assert gva.ambiguous_base_fraction("") == 0.0


# --------------------------------------------------------------------------- #
# reverse_complement / translate / resolve_cds_strand                        #
# (the real-data-driven orientation and frame-offset fix)                    #
# --------------------------------------------------------------------------- #

class TestReverseComplement:
    def test_basic(self, gva):
        assert gva.reverse_complement("ACGT") == "ACGT"
        assert gva.reverse_complement("ATTG") == "CAAT"

    def test_preserves_n(self, gva):
        assert gva.reverse_complement("ACGN") == "NCGT"


class TestTranslate:
    def test_standard_codons(self, gva):
        assert gva.translate("ATGAAATAA") == "MK*"

    def test_unknown_codon_is_x(self, gva):
        assert gva.translate("NNNATG") == "XM"

    def test_trailing_partial_codon_is_ignored(self, gva):
        assert gva.translate("ATGAA") == "M"


# Shared synthetic coding sequence: M A K E *
_CORRECT_CDS = "ATGGCTAAAGAATAA"
_CORRECT_AA = "MAKE*"


class TestResolveCdsStrand:
    def test_already_correct_forward_frame0(self, gva):
        resolved = gva.resolve_cds_strand(_CORRECT_CDS, _CORRECT_AA)
        assert resolved == _CORRECT_CDS

    def test_ambisense_minus_strand_is_recovered(self, gva):
        """Mirrors the real Lassa S-segment finding: cds.fna gives the
        un-reverse-complemented span for a minus-strand CDS."""
        given = gva.reverse_complement(_CORRECT_CDS)
        resolved = gva.resolve_cds_strand(given, _CORRECT_AA)
        assert resolved == _CORRECT_CDS
        assert gva.translate(resolved).rstrip("*") == _CORRECT_AA.rstrip("*")

    def test_ribosomal_frameshift_offset_is_recovered(self, gva):
        """Mirrors the real Norovirus ORF1/ORF2 junction finding: the
        recorded range can start a couple of bases before the true
        translation start."""
        given = "CC" + _CORRECT_CDS
        resolved = gva.resolve_cds_strand(given, _CORRECT_AA)
        assert resolved == _CORRECT_CDS

    def test_no_orientation_matches_returns_none(self, gva):
        garbage = "TTTTTTTTTTTTTTTT"
        assert gva.resolve_cds_strand(garbage, _CORRECT_AA) is None

    def test_min_score_threshold_is_respected(self, gva):
        # A sequence that matches perfectly should pass even a strict threshold.
        resolved = gva.resolve_cds_strand(_CORRECT_CDS, _CORRECT_AA, min_score=0.99)
        assert resolved == _CORRECT_CDS


# --------------------------------------------------------------------------- #
# write_genome_fastas / write_metadata_tsv                                   #
# --------------------------------------------------------------------------- #

class TestWriteGenomeFastas:
    def _genomes_cds_prot(self):
        genomes = {
            "g1": {"entries": [
                {"nt_seqid": "nt1", "prot_seqid": "aa1"},
                {"nt_seqid": "nt2", "prot_seqid": "aa2"},
            ]},
        }
        cds_seqs = {
            "nt1": ("nt1", _CORRECT_CDS),
            "nt2": ("nt2", gva_reverse(_CORRECT_CDS)),
        }
        prot_seqs = {
            "aa1": ("aa1", "MAKE"),
            "aa2": ("aa2", "MAKE"),
        }
        return genomes, cds_seqs, prot_seqs

    def test_writes_resolved_sequences(self, gva, tmp_path):
        genomes = {
            "g1": {"entries": [{"nt_seqid": "nt1", "prot_seqid": "aa1"}]},
        }
        cds_seqs = {"nt1": ("nt1", gva.reverse_complement(_CORRECT_CDS))}
        prot_seqs = {"aa1": ("aa1", "MAKE")}
        nucl_dir = tmp_path / "nucl"
        prot_dir = tmp_path / "prot"
        written = gva.write_genome_fastas(genomes, cds_seqs, prot_seqs,
                                           nucl_dir, prot_dir)
        assert written == {"g1": 1}
        nt_out = (nucl_dir / "g1.fna").read_text()
        assert _CORRECT_CDS in nt_out
        assert (prot_dir / "g1.faa").read_text().strip() == ">g1_1\nMAKE"

    def test_missing_sequence_is_skipped_not_fatal(self, gva, tmp_path):
        genomes = {
            "g1": {"entries": [
                {"nt_seqid": "missing", "prot_seqid": "aa1"},
                {"nt_seqid": "nt1", "prot_seqid": "aa1"},
            ]},
        }
        cds_seqs = {"nt1": ("nt1", _CORRECT_CDS)}
        prot_seqs = {"aa1": ("aa1", "MAKE")}
        written = gva.write_genome_fastas(genomes, cds_seqs, prot_seqs,
                                           tmp_path / "nucl", tmp_path / "prot")
        assert written == {"g1": 1}

    def test_genome_with_zero_usable_genes_is_dropped(self, gva, tmp_path):
        genomes = {
            "g1": {"entries": [{"nt_seqid": "nt1", "prot_seqid": "aa1"}]},
        }
        cds_seqs = {"nt1": ("nt1", "TTTTTTTTTTTTTTTT")}  # won't resolve
        prot_seqs = {"aa1": ("aa1", "MAKE")}
        written = gva.write_genome_fastas(genomes, cds_seqs, prot_seqs,
                                           tmp_path / "nucl", tmp_path / "prot")
        assert written == {}
        assert not (tmp_path / "nucl" / "g1.fna").exists()

    def test_internal_stop_gene_is_skipped(self, gva, tmp_path):
        # Resolves (matches aa) but has an earlier in-frame stop too --
        # use a dedicated aa/nt pair with a premature stop baked in.
        nt_with_early_stop = "ATGTAAGCTAAAGAATAA"  # M * A K E *
        genomes = {
            "g1": {"entries": [{"nt_seqid": "nt1", "prot_seqid": "aa1"}]},
        }
        cds_seqs = {"nt1": ("nt1", nt_with_early_stop)}
        # aa_seq matches the full translation (incl. internal stop) so that
        # resolve_cds_strand accepts it before the has_internal_stop check runs.
        prot_seqs = {"aa1": ("aa1", gva.translate(nt_with_early_stop))}
        written = gva.write_genome_fastas(genomes, cds_seqs, prot_seqs,
                                           tmp_path / "nucl", tmp_path / "prot")
        assert written == {}

    def test_too_ambiguous_gene_is_skipped(self, gva, tmp_path):
        nt_ambiguous = "ATGNNNNNNNNNNNNNN"
        genomes = {
            "g1": {"entries": [{"nt_seqid": "nt1", "prot_seqid": "aa1"}]},
        }
        cds_seqs = {"nt1": ("nt1", nt_ambiguous)}
        prot_seqs = {"aa1": ("aa1", gva.translate(nt_ambiguous))}
        written = gva.write_genome_fastas(genomes, cds_seqs, prot_seqs,
                                           tmp_path / "nucl", tmp_path / "prot")
        assert written == {}

    def test_gene_ids_renumbered_sequentially(self, gva, tmp_path):
        genomes = {
            "g1": {"entries": [
                {"nt_seqid": "nt1", "prot_seqid": "aa1"},
                {"nt_seqid": "nt2", "prot_seqid": "aa2"},
            ]},
        }
        cds_seqs = {
            "nt1": ("nt1", _CORRECT_CDS),
            "nt2": ("nt2", _CORRECT_CDS),
        }
        prot_seqs = {
            "aa1": ("aa1", "MAKE"),
            "aa2": ("aa2", "MAKE"),
        }
        gva.write_genome_fastas(genomes, cds_seqs, prot_seqs,
                                 tmp_path / "nucl", tmp_path / "prot")
        nt_out = (tmp_path / "nucl" / "g1.fna").read_text()
        assert ">g1_1" in nt_out and ">g1_2" in nt_out


def gva_reverse(seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[b] for b in reversed(seq))


class TestWriteMetadataTsv:
    def test_columns_and_order(self, gva, tmp_path):
        genomes = {
            "g1": {"representative_accession": "LC936410.1"},
        }
        data_report_recs = [{
            "accession": "LC936410.1",
            "virus": {"organismName": "Norovirus GII", "taxId": 142786},
            "location": {"geographicLocation": "Japan"},
        }]
        out = tmp_path / "all_asm_acc_metadata"
        gva.write_metadata_tsv(["g1"], genomes, data_report_recs, out)
        line = out.read_text().strip()
        fields = line.split(" ")
        # Exactly 13 columns (8 original + 5 NA submitter-provenance/lat_lon
        # placeholders, see parse_asm_metadata.py's project_metadata) -- a
        # raw space in species/geo would silently produce more, shifting
        # every later column for a whitespace-split reader (see
        # annotate_tree_by_country.py's parse_metadata).
        assert len(fields) == 13
        assert fields[0] == "g1"  # accession == genome_id
        assert fields[1] == "Norovirus_GII"  # space collapsed to underscore
        assert fields[3] == "g1"  # biosample == genome_id (no real biosample)
        assert "142786" in fields
        assert fields[7] == "Japan"  # geo_loc, column 8
        assert fields[8:13] == ["NA", "NA", "NA", "NA", "NA"]

    def test_species_and_geo_spaces_are_sanitized(self, gva, tmp_path):
        genomes = {"g1": {"representative_accession": "ACC1"}}
        data_report_recs = [{
            "accession": "ACC1",
            "virus": {"organismName": "Some Virus Name", "taxId": 1},
            "location": {"geographicLocation": "USA: California"},
        }]
        out = tmp_path / "all_asm_acc_metadata"
        gva.write_metadata_tsv(["g1"], genomes, data_report_recs, out)
        fields = out.read_text().strip().split(" ")
        assert len(fields) == 13
        assert fields[1] == "Some_Virus_Name"
        assert fields[7] == "USA:_California"

    def test_missing_data_report_entry_writes_blank_fields(self, gva, tmp_path):
        genomes = {"g1": {"representative_accession": "UNKNOWN_ACC"}}
        out = tmp_path / "all_asm_acc_metadata"
        gva.write_metadata_tsv(["g1"], genomes, [], out)
        line = out.read_text().rstrip("\n")
        fields = line.split(" ")
        assert len(fields) == 13
        assert fields[0] == "g1"
        assert fields[4] == ""  # taxid blank (no virus-package equivalent)
        assert fields[5:7] == ["NA", "NA"]  # taxcheck/excl_from_refseq
        assert fields[8:13] == ["NA", "NA", "NA", "NA", "NA"]
