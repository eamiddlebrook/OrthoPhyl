"""Unit tests for python_scripts/viral_annotate.py.

Fixtures under tests/fixtures/viral_annotate/ are REAL captured output from
running v-annotate.pl (VADR 1.6.4, dengue flavi model library) and vigor4
(built from source, dengue_txid12637_db) on the same small dengue test FASTA
-- see CLAUDE.md's viral-annotation plan for how they were generated. No tool
subprocesses run in this file; only the pure parsing/translation/picking
functions are exercised.
"""
import os

import pytest

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "fixtures", "viral_annotate")


@pytest.fixture
def va():
    # viral_annotate.py is a standalone script, not an importable package --
    # load it by path the same way tests/conftest.py does for the other
    # root-level scripts.
    import importlib.util
    path = os.path.join(os.path.dirname(__file__), "..", "..", "python_scripts", "viral_annotate.py")
    spec = importlib.util.spec_from_file_location("viral_annotate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_fixture(*parts):
    with open(os.path.join(FIXTURES, *parts)) as fh:
        return fh.read()


class TestTranslate:
    def test_standard_table_known_codons(self, va):
        assert va.translate("ATG") == "M"
        assert va.translate("TAA") == "*"
        assert va.translate("TGG") == "W"
        assert va.translate("GGG") == "G"

    def test_truncates_incomplete_trailing_codon(self, va):
        assert va.translate("ATGAA") == "M"

    def test_unknown_codon_maps_to_x(self, va):
        assert va.translate("NNN") == "X"


class TestParseFasta:
    def test_basic_multi_record(self, va):
        text = ">a desc1\nACGT\nACGT\n>b\nTTTT\n"
        recs = va.parse_fasta(text)
        assert recs == [("a desc1", "ACGTACGT"), ("b", "TTTT")]


class TestVadrSqc:
    def test_parses_real_fixture(self, va):
        text = read_fixture("vadr_dengue", "va_test2.vadr.sqc")
        rows = va.parse_sqc(text)
        assert len(rows) == 5
        names = [r["name"] for r in rows]
        assert "gi|770146038|gb|KM031079.1|" in names
        assert all(r["pf"] == "PASS" for r in rows)
        assert all(va.is_classified(r) for r in rows)

    def test_is_classified_false_for_dash_model(self, va):
        assert va.is_classified({"model1": "-"}) is False
        assert va.is_classified({"model1": ""}) is False
        assert va.is_classified({"model1": "NC_001474"}) is True

    def test_is_classified_true_even_on_fail(self, va):
        # per CLAUDE.md: classified means "a model was assigned", not "passed QC"
        assert va.is_classified({"model1": "NC_001474", "pf": "FAIL"}) is True


class TestChooseVadrLibrary:
    def test_picks_first_classifying_library(self, va):
        sqc_pass = [{"model1": "NC_001474", "pf": "PASS"}]
        sqc_fail = [{"model1": "-", "pf": "FAIL"}]
        results = [("calici_dir", "calici", sqc_fail), ("flavi_dir", "flavi", sqc_pass)]
        winner = va.choose_vadr_library(results)
        assert winner == ("flavi_dir", "flavi", sqc_pass)

    def test_returns_none_when_nothing_classifies(self, va):
        sqc_fail = [{"model1": "-", "pf": "FAIL"}]
        results = [("calici_dir", "calici", sqc_fail), ("flavi_dir", "flavi", sqc_fail)]
        assert va.choose_vadr_library(results) is None

    def test_order_matters_first_classifier_wins_even_if_later_also_classifies(self, va):
        sqc_pass = [{"model1": "X", "pf": "PASS"}]
        results = [("a", "a", sqc_pass), ("b", "b", sqc_pass)]
        assert va.choose_vadr_library(results)[1] == "a"


class TestParseVadrFeatureFasta:
    def test_real_fixture_headers(self, va):
        text = read_fixture("vadr_dengue", "va_test2.vadr.NC_001474.CDS.1.fa")
        feats = va.parse_vadr_feature_fasta(text)
        assert len(feats) == 4
        seq_id, label, seq = feats[0]
        assert seq_id == "gi|770146038|gb|KM031079.1|"
        assert label == "CDS.1"
        assert len(seq) == 515

    def test_length_matches_header_coords_for_all_records(self, va):
        text = read_fixture("vadr_dengue", "va_test2.vadr.NC_001477.CDS.1.fa")
        feats = va.parse_vadr_feature_fasta(text)
        for seq_id, label, seq in feats:
            assert len(seq) == 1485  # the single NC_001477-assigned seq in this fixture


class TestParseAlt:
    def test_empty_alert_table(self, va):
        text = read_fixture("vadr_dengue", "va_test2.vadr.alt")
        rows = va.parse_alt(text)
        assert rows == []

    def test_synthetic_fatal_alert_row(self, va):
        # Real column layout confirmed against VADR 1.6.4's .alt header; no
        # genome in hand actually triggers a fatal alert, so this row is
        # hand-built in the confirmed format to exercise the parser.
        text = (
            "#     seq          ftr   ftr   ftr  alert        alert           seq  seq     mdl  mdl  alert \n"
            "#idx  name  model  type  name  idx  code   fail  description  coords  len  coords  len  detail\n"
            "#---  ----  -----  ----  ----  ---  -----  ----  -----------  ------  ---  ------  ---  ------\n"
            "1     gi|123|gb|AB1.1|  NC_001474  CDS  GP  1  mutstart  yes  MUTATION_AT_START  1..3  3393  1..3  3393  expected ATG\n"
        )
        rows = va.parse_alt(text)
        assert len(rows) == 1
        assert rows[0]["code"] == "mutstart"
        assert rows[0]["fail"] == "yes"


class TestVigorParsing:
    def test_cds_pep_ids_from_real_fixture(self, va):
        cds = va.parse_vigor_fasta(read_fixture("vigor4_dengue", "test.cds"))
        pep = va.parse_vigor_fasta(read_fixture("vigor4_dengue", "test.pep"))
        assert set(cds.keys()) == {"NC_001477.1", "NC_001474.1"}
        # .pep has the matching CDS-level ids PLUS 12 mat_peptide entries per genome
        assert len(pep) == len(cds) + 24

    def test_filter_vigor_pep_drops_mat_peptides(self, va):
        cds = va.parse_vigor_fasta(read_fixture("vigor4_dengue", "test.cds"))
        pep = va.parse_vigor_fasta(read_fixture("vigor4_dengue", "test.pep"))
        filtered = va.filter_vigor_pep(cds, pep)
        assert set(filtered.keys()) == set(cds.keys())
        for rid in cds:
            # pep is the CDS minus the stop codon
            assert len(filtered[rid]) == len(cds[rid]) // 3 - 1


class TestPickVigorDb:
    def test_picks_highest_score(self, va):
        assert va.pick_vigor_db({"dengue_txid12637_db": 500.0, "zika_db": 12.0}) == "dengue_txid12637_db"

    def test_ties_broken_alphabetically(self, va):
        assert va.pick_vigor_db({"zika_db": 100.0, "dengue_db": 100.0}) == "dengue_db"

    def test_empty_or_nonpositive_returns_none(self, va):
        assert va.pick_vigor_db({}) is None
        assert va.pick_vigor_db({"x": 0.0, "y": -5.0}) is None


class TestParseDiamondOutfmt6:
    def test_aggregates_by_db_prefix(self, va):
        text = (
            "q1\tdengue_txid12637_db|NP_1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
            "q1\tdengue_txid12637_db|NP_2\t90\t100\t0\t0\t1\t100\t1\t100\t1e-50\t150.0\n"
            "q1\tzika_db|NP_3\t50\t100\t0\t0\t1\t100\t1\t100\t1e-5\t20.0\n"
        )
        scores = va.parse_diamond_outfmt6(text)
        assert scores == {"dengue_txid12637_db": 350.0, "zika_db": 20.0}


class TestRunVadrOutputPaths:
    def test_out_dir_basename_equals_out_prefix(self, va, tmp_path, monkeypatch):
        # v-annotate.pl names every output file after the basename of the
        # output directory argument it's given (e.g. "<basename>.vadr.sqc").
        # out_dir must therefore be exactly workdir/out_prefix -- a prior bug
        # appended an extra ".vadr_<mkey>" suffix onto out_dir, so
        # annotate_via_vadr's "<out_prefix>.vadr.sqc" lookup never matched
        # anything VADR actually wrote, silently treating every real
        # classification as "not classified". Confirmed against a real
        # v-annotate.pl 1.6.4 run before fixing.
        monkeypatch.setattr(va, "_run", lambda cmd, env_run_prefix="": None)
        out_dir, _ = va.run_vadr("genome.fa", "mdir", "mkey", str(tmp_path), "stem.mkey")
        assert out_dir == os.path.join(str(tmp_path), "stem.mkey")
        assert os.path.basename(out_dir) == "stem.mkey"

    def test_annotate_via_vadr_finds_sqc_written_at_real_vadr_path(self, va, tmp_path, monkeypatch):
        def fake_run_vadr(genome_fa, mdir, mkey, workdir, out_prefix, vadr_options="", env_run_prefix=""):
            out_dir = os.path.join(workdir, out_prefix)
            os.makedirs(out_dir, exist_ok=True)
            # mirrors real VADR: file named after basename(out_dir), matching out_prefix
            with open(os.path.join(out_dir, f"{out_prefix}.vadr.sqc"), "w") as fh:
                fh.write("1  seq1  100  PASS  yes  modelA  grp  1  1  1  1  1  1  1  +  -  -  -  -  -  -\n")
            class R: returncode = 0
            return out_dir, R()
        monkeypatch.setattr(va, "run_vadr", fake_run_vadr)
        result = va.annotate_via_vadr("genome.fa", "stem", ["mdir:mkey"], str(tmp_path))
        assert result["classified"] is True
        assert result["mkey"] == "mkey"


class TestNormalizeRecords:
    def test_renames_in_order_and_drops_internal_stop(self, va):
        nt = [("a", "ATGAAATAA"), ("b", "ATGAAATAA")]
        pep = [("a", "MK*"), ("b", "M*K")]  # b has an internal stop
        trans_out, prots_out, n_dropped = va.normalize_records("g", nt, pep)
        assert n_dropped == 1
        assert trans_out == [("g_1", "ATGAAATAA")]
        assert prots_out == [("g_1", "MK")]  # trailing stop stripped

    def test_empty_input(self, va):
        trans_out, prots_out, n_dropped = va.normalize_records("g", [], [])
        assert trans_out == []
        assert prots_out == []
        assert n_dropped == 0
