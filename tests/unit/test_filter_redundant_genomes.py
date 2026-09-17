"""Tests for python_scripts/filter_redundant_genomes.py.

filter_for_redundancy in gather_filter_asms.sh de-versions genome stems with an
unanchored `grep $I | sort | tail -n 1` over a sed-truncated accession number. On
non-NCBI (user-supplied) filenames this silently collapses distinct genomes and is
locale-dependent. This module replaces that logic: only names that actually look like
NCBI accessions are de-versioned; everything else passes through untouched.
"""

import pytest


@pytest.fixture
def fr(filter_redundant_module):
    return filter_redundant_module


class TestPartitionByAccession:
    def test_splits_accession_like_from_passthrough(self, fr):
        names = ["GCF_000001405.40", "isolate_1", "GCA_000001405.29", "PA.contigs"]
        acc, other = fr.partition_by_accession(names)
        assert set(acc) == {"GCF_000001405.40", "GCA_000001405.29"}
        assert set(other) == {"isolate_1", "PA.contigs"}

    def test_bare_gc_prefix_without_full_pattern_is_passthrough(self, fr):
        # Not a full accession -- no version, or wrong prefix -- must not be treated
        # as one.
        names = ["GCF_000001405", "GCX_000001405.1", "GCF_abc.1"]
        acc, other = fr.partition_by_accession(names)
        assert acc == []
        assert set(other) == set(names)


class TestDedupeAccessions:
    def test_prefers_refseq_over_genbank(self, fr):
        kept = fr.dedupe_accessions(["GCF_000001405.40", "GCA_000001405.29"])
        assert kept == ["GCF_000001405.40"]

    def test_prefers_highest_version_within_same_prefix(self, fr):
        kept = fr.dedupe_accessions(["GCA_000001405.1", "GCA_000001405.29", "GCA_000001405.15"])
        assert kept == ["GCA_000001405.29"]

    def test_distinct_base_numbers_all_kept(self, fr):
        names = ["GCF_000001405.40", "GCF_000001635.27", "GCF_900004695.1"]
        kept = fr.dedupe_accessions(names)
        assert sorted(kept) == sorted(names)

    def test_no_refseq_falls_back_to_genbank(self, fr):
        kept = fr.dedupe_accessions(["GCA_000002.1", "GCA_000002.2"])
        assert kept == ["GCA_000002.2"]


class TestFilterRedundant:
    def test_reproduces_the_bug_scenario_without_the_bug(self, fr):
        # The exact scenario that broke the bash implementation: user-style filenames
        # with numeric substrings that an unanchored grep collapsed. All 6 distinct
        # genomes must survive here.
        names = [
            "isolate_1", "isolate_2", "isolate_11", "isolate_12",
            "PA.contigs", "PA.scaffolds",
        ]
        kept = fr.filter_redundant(names)
        assert sorted(kept) == sorted(names)

    def test_mixed_accession_and_user_names(self, fr):
        names = [
            "GCF_000001405.40", "GCA_000001405.29",  # redundant pair -> keep GCF only
            "GCF_000002.1",                            # distinct accession
            "isolate_1", "isolate_11",                 # must both survive
        ]
        kept = fr.filter_redundant(names)
        assert sorted(kept) == sorted([
            "GCF_000001405.40", "GCF_000002.1", "isolate_1", "isolate_11",
        ])

    def test_real_ncbi_accessions_dedupe_as_before(self, fr):
        names = [
            "GCF_000001405.40", "GCF_000001635.27", "GCF_900004695.1",
            "GCA_000001405.29",  # same base as the first GCF -> dropped
        ]
        kept = fr.filter_redundant(names)
        assert sorted(kept) == sorted([
            "GCF_000001405.40", "GCF_000001635.27", "GCF_900004695.1",
        ])

    def test_empty_input(self, fr):
        assert fr.filter_redundant([]) == []
