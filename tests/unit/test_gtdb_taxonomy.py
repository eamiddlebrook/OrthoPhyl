"""Tests for the GTDBTaxonomy class (plan section 1).

GTDBTaxonomy is defined twice -- once in the router and once in the database creator --
with near-identical logic. We parametrize over both copies so the pair cannot silently
drift apart. Methods that exist only on the router copy (is_within_clade,
get_taxonomy_string_at_rank) are tested against that copy alone.
"""

import pytest

from conftest import GAIELLALES_TAX, ESCHERICHIA_TAX


@pytest.fixture(params=["router", "db_creator"])
def GTDBTaxonomy(request, router_module, db_creator_module):
    """The GTDBTaxonomy class from each module that defines it."""
    module = router_module if request.param == "router" else db_creator_module
    return module.GTDBTaxonomy


class TestParsing:
    def test_parse_full_lineage(self, GTDBTaxonomy):
        tax = GTDBTaxonomy(
            "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S"
        )
        assert tax.get_rank("d") == "Bacteria"
        assert tax.get_rank("p") == "P"
        assert tax.get_rank("c") == "C"
        assert tax.get_rank("o") == "O"
        assert tax.get_rank("f") == "F"
        assert tax.get_rank("g") == "G"
        assert tax.get_rank("s") == "S"

    def test_parse_empty_species(self, GTDBTaxonomy):
        tax = GTDBTaxonomy(ESCHERICHIA_TAX + ";s__")
        assert tax.get_rank("s") is None
        assert tax.get_rank("g") == "Escherichia"

    def test_parse_strips_whitespace_after_semicolon(self, GTDBTaxonomy):
        tax = GTDBTaxonomy("d__Bacteria; p__Pseudomonadota")
        assert tax.get_rank("p") == "Pseudomonadota"

    def test_parse_ignores_malformed_parts(self, GTDBTaxonomy):
        tax = GTDBTaxonomy("garbage;p__Actinomycetota;also_garbage")
        assert tax.get_rank("p") == "Actinomycetota"
        assert tax.get_rank("d") is None

    def test_empty_string_yields_no_levels(self, GTDBTaxonomy):
        tax = GTDBTaxonomy("")
        assert tax.get_most_specific_rank() is None


class TestRankQueries:
    def test_get_most_specific_rank_genus_when_species_empty(self, GTDBTaxonomy):
        tax = GTDBTaxonomy(ESCHERICHIA_TAX + ";s__")
        assert tax.get_most_specific_rank() == "g"

    def test_get_most_specific_rank_species_when_defined(self, GTDBTaxonomy):
        tax = GTDBTaxonomy(ESCHERICHIA_TAX + ";s__Escherichia_coli")
        assert tax.get_most_specific_rank() == "s"

    def test_get_most_specific_rank_none_when_all_empty(self, GTDBTaxonomy):
        tax = GTDBTaxonomy("d__;p__;c__")
        assert tax.get_most_specific_rank() is None

    @pytest.mark.parametrize(
        "rank,expected",
        [
            ("d", "domain"),
            ("p", "phylum"),
            ("c", "class"),
            ("o", "order"),
            ("f", "family"),
            ("g", "genus"),
            ("s", "species"),
        ],
    )
    def test_get_rank_name(self, GTDBTaxonomy, rank, expected):
        tax = GTDBTaxonomy(GAIELLALES_TAX)
        assert tax.get_rank_name(rank) == expected

    def test_str_returns_raw(self, GTDBTaxonomy):
        tax = GTDBTaxonomy("  " + GAIELLALES_TAX + "  ")
        # __init__ strips the raw string.
        assert str(tax) == GAIELLALES_TAX


# --------------------------------------------------------------------------- #
# Router-only methods                                                         #
# --------------------------------------------------------------------------- #

@pytest.fixture
def RouterTaxonomy(router_module):
    return router_module.GTDBTaxonomy


class TestIsWithinClade:
    def test_query_within_family_clade(self, RouterTaxonomy):
        query = RouterTaxonomy(ESCHERICHIA_TAX + ";s__Escherichia_coli")
        clade = RouterTaxonomy(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae"
        )
        assert query.is_within_clade(clade, "f") is True

    def test_query_within_genus_clade(self, RouterTaxonomy):
        query = RouterTaxonomy(ESCHERICHIA_TAX + ";s__Escherichia_coli")
        clade = RouterTaxonomy(ESCHERICHIA_TAX)
        assert query.is_within_clade(clade, "g") is True

    def test_different_phylum_not_within(self, RouterTaxonomy):
        query = RouterTaxonomy(GAIELLALES_TAX)
        clade = RouterTaxonomy(ESCHERICHIA_TAX)
        assert query.is_within_clade(clade, "g") is False

    def test_query_shallower_than_clade_rank(self, RouterTaxonomy):
        # Query only defined to order; clade is at genus -> query lacks the genus.
        query = RouterTaxonomy(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales"
        )
        clade = RouterTaxonomy(ESCHERICHIA_TAX)
        assert query.is_within_clade(clade, "g") is False


class TestTaxonomyStringAtRank:
    def test_truncates_at_family(self, RouterTaxonomy):
        tax = RouterTaxonomy(ESCHERICHIA_TAX + ";s__Escherichia_coli")
        result = tax.get_taxonomy_string_at_rank("f")
        assert result.endswith("f__Enterobacteriaceae")
        assert "g__" not in result

    def test_renders_empty_ranks(self, RouterTaxonomy):
        tax = RouterTaxonomy("d__Bacteria;p__Pseudomonadota")
        result = tax.get_taxonomy_string_at_rank("c")
        # class is undefined -> rendered as bare "c__"
        assert result == "d__Bacteria;p__Pseudomonadota;c__"
