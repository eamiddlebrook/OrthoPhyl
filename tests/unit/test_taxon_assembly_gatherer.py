"""Tests for utils/taxon_assembly_gatherer.py, focused on the wrapper contract (B4).

Bug B4 was a mismatch between the wrapper's calls and the gatherer's real API. These
tests pin the reconciled contract:

  * constructor accepts taxon=, rank=, output_dir=
  * query_ncbi()          exists (alias for gather_assemblies)
  * download_assemblies()  exists
  * get_taxonomy_string()  exists
  * .taxid / .taxon_rank   attributes exist

Plus unit coverage of the two newly-added methods with network mocked out.
"""

import inspect
import gzip
from pathlib import Path

import pytest


# --------------------------------------------------------------------------- #
# A tiny in-memory NCBITaxonomy substitute so we don't touch the network or    #
# need a taxdump fixture just to exercise the new methods.                      #
# --------------------------------------------------------------------------- #

class FakeTaxonomy:
    def __init__(self):
        # 561 = Escherichia (genus) -> ... -> 2 Bacteria (superkingdom)
        self._lineage = {
            "genus": "Escherichia",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Pseudomonadota",
            "superkingdom": "Bacteria",
        }
        # Non-viral by default: a bacterial chain that never includes 10239.
        self._taxid_chain = ["561", "543", "91347", "1224", "3379134", "2"]

    def get_lineage(self, taxid):
        return dict(self._lineage)

    def get_taxid_chain(self, taxid):
        return list(self._taxid_chain)


@pytest.fixture
def gatherer(taxon_gatherer_module, monkeypatch, tmp_path):
    """A TaxonAssemblyGatherer with taxonomy DB bootstrap and resolution stubbed out."""
    Gatherer = taxon_gatherer_module.TaxonAssemblyGatherer

    # Avoid downloading/loading a real taxdump.
    monkeypatch.setattr(Gatherer, "_ensure_taxonomy_database", lambda self: None)
    monkeypatch.setattr(
        taxon_gatherer_module, "NCBITaxonomy", lambda taxdump_dir: FakeTaxonomy()
    )

    g = Gatherer.__new__(Gatherer)
    # Manually set the attributes the ctor would, without network I/O.
    g.taxon = "Escherichia"
    g.output_dir = tmp_path / "out"
    g.output_dir.mkdir(parents=True, exist_ok=True)
    g.taxdump_dir = tmp_path / "taxdump"
    g.taxonomy = FakeTaxonomy()
    g.taxid = "561"
    g.taxon_name = "Escherichia"
    g.taxon_rank = "genus"
    g.is_viral = False
    g.since_date = None
    g.assemblies = []
    return g


class TestWrapperContract:
    """B4: the reconciled interface the wrapper depends on."""

    def test_ctor_accepts_wrapper_kwargs(self, taxon_gatherer_module):
        sig = inspect.signature(taxon_gatherer_module.TaxonAssemblyGatherer.__init__)
        params = sig.parameters
        assert "taxon" in params
        assert "rank" in params
        assert "output_dir" in params
        # The old broken name must NOT be what the wrapper has to pass.
        assert "taxon_name" not in params

    @pytest.mark.parametrize(
        "method", ["query_ncbi", "download_assemblies", "get_taxonomy_string"]
    )
    def test_required_methods_exist(self, taxon_gatherer_module, method):
        assert hasattr(taxon_gatherer_module.TaxonAssemblyGatherer, method)

    def test_query_ncbi_aliases_gather(self, gatherer, monkeypatch):
        sentinel = [{"accession": "GCF_1.1"}]
        monkeypatch.setattr(gatherer, "gather_assemblies", lambda: sentinel)
        assert gatherer.query_ncbi() is sentinel


class TestGetTaxonomyString:
    def test_builds_gtdb_string(self, gatherer):
        s = gatherer.get_taxonomy_string()
        assert s == (
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia"
        )

    def test_empty_lineage_returns_empty(self, gatherer):
        gatherer.taxonomy._lineage = {}
        assert gatherer.get_taxonomy_string() == ""

    def test_stops_at_deepest_defined_rank(self, gatherer):
        gatherer.taxonomy._lineage = {
            "superkingdom": "Bacteria",
            "phylum": "Pseudomonadota",
        }
        s = gatherer.get_taxonomy_string()
        assert s == "d__Bacteria;p__Pseudomonadota"


class TestGetTaxidChain:
    """NCBITaxonomy.get_taxid_chain -- the TaxID-based path used for viral
    detection, as opposed to get_lineage's rank-name keying."""

    def _taxonomy(self, taxon_gatherer_module):
        tax = taxon_gatherer_module.NCBITaxonomy.__new__(
            taxon_gatherer_module.NCBITaxonomy
        )
        # norovirus (142786) -> Caliciviridae (11974) -> ... -> Viruses (10239) -> root (1)
        tax.nodes = {
            "142786": {"parent": "11974", "rank": "species"},
            "11974": {"parent": "10239", "rank": "family"},
            "10239": {"parent": "1", "rank": "superkingdom"},
            "1": {"parent": "1", "rank": "no rank"},
            "561": {"parent": "543", "rank": "genus"},
            "543": {"parent": "2", "rank": "family"},
            "2": {"parent": "2", "rank": "superkingdom"},
        }
        tax.names = {}
        tax.name_to_taxid = {}
        return tax

    def test_walks_full_chain_to_root(self, taxon_gatherer_module):
        tax = self._taxonomy(taxon_gatherer_module)
        assert tax.get_taxid_chain("142786") == ["142786", "11974", "10239", "1"]

    def test_bacterial_chain_excludes_virus_taxid(self, taxon_gatherer_module):
        tax = self._taxonomy(taxon_gatherer_module)
        chain = tax.get_taxid_chain("561")
        assert "10239" not in chain

    def test_unknown_taxid_returns_single_element_or_empty(self, taxon_gatherer_module):
        tax = self._taxonomy(taxon_gatherer_module)
        assert tax.get_taxid_chain("999999") == []


class TestIsViralTaxid:
    def test_true_when_virus_taxid_present(self, taxon_gatherer_module):
        assert taxon_gatherer_module.is_viral_taxid(["142786", "11974", "10239", "1"])

    def test_false_for_bacterial_chain(self, taxon_gatherer_module):
        assert not taxon_gatherer_module.is_viral_taxid(["561", "543", "2"])

    def test_false_for_empty_chain(self, taxon_gatherer_module):
        assert not taxon_gatherer_module.is_viral_taxid([])


class TestRenderGtdbLineageViral:
    def test_is_viral_forces_domain_viruses(self, taxon_gatherer_module):
        lineage = {
            "family": "Caliciviridae",
            "genus": "Norovirus",
            "species": "Norwalk virus",
        }
        s = taxon_gatherer_module.render_gtdb_lineage(lineage, is_viral=True)
        assert s.startswith("d__Viruses;")

    def test_overrides_whatever_rank_mapping_produced(self, taxon_gatherer_module):
        # Even if an NCBI rank happened to map to 'd' (e.g. a stray
        # 'superkingdom'/'domain' entry under Viruses), is_viral wins.
        lineage = {"superkingdom": "SomethingWeird", "genus": "Norovirus"}
        s = taxon_gatherer_module.render_gtdb_lineage(lineage, is_viral=True)
        assert s.startswith("d__Viruses;")
        assert "SomethingWeird" not in s

    def test_non_viral_unaffected(self, taxon_gatherer_module):
        lineage = {"superkingdom": "Bacteria", "genus": "Escherichia"}
        s = taxon_gatherer_module.render_gtdb_lineage(lineage, is_viral=False)
        assert s.startswith("d__Bacteria;")


class TestGetTaxonomyStringViral:
    def test_viral_taxon_renders_d_viruses(self, gatherer):
        gatherer.is_viral = True
        gatherer.taxonomy._lineage = {
            "family": "Caliciviridae",
            "genus": "Norovirus",
        }
        s = gatherer.get_taxonomy_string()
        assert s == "d__Viruses;p__;c__;o__;f__Caliciviridae;g__Norovirus"

    def test_is_viral_set_from_taxid_chain_in_ctor(
        self, taxon_gatherer_module, monkeypatch, tmp_path
    ):
        """End-to-end through __init__: a taxid chain containing 10239 must
        set self.is_viral, driving get_taxonomy_string's d__Viruses override."""
        Gatherer = taxon_gatherer_module.TaxonAssemblyGatherer

        class ViralFakeTaxonomy(FakeTaxonomy):
            def __init__(self):
                super().__init__()
                self._lineage = {"family": "Caliciviridae", "genus": "Norovirus"}
                self._taxid_chain = ["142786", "11974", "10239", "1"]

            def resolve_taxon(self, taxon):
                return "142786"

            def get_rank(self, taxid):
                return "species"

            def get_name(self, taxid):
                return "Norovirus"

        monkeypatch.setattr(Gatherer, "_ensure_taxonomy_database", lambda self: None)
        monkeypatch.setattr(
            taxon_gatherer_module, "NCBITaxonomy", lambda taxdump_dir: ViralFakeTaxonomy()
        )

        g = Gatherer(taxon="Norovirus", output_dir=tmp_path / "out")
        assert g.is_viral is True
        assert g.get_taxonomy_string().startswith("d__Viruses;")


class TestDownloadAssemblies:
    def test_downloads_and_gunzips(self, gatherer, monkeypatch, tmp_path):
        # Fake urlretrieve: write a gzipped FASTA to the requested destination.
        def fake_urlretrieve(url, dest):
            with gzip.open(dest, "wb") as fh:
                fh.write(b">contig1\nACGTACGT\n")

        import taxon_assembly_gatherer as tag_mod  # registered by conftest loader
        monkeypatch.setattr(tag_mod.urllib.request, "urlretrieve", fake_urlretrieve)

        assemblies = [
            {
                "accession": "GCF_000005845.2",
                "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2",
            }
        ]
        dl = tmp_path / "downloads"
        written = gatherer.download_assemblies(assemblies, dl)

        out = dl / "GCF_000005845.2.fna"
        assert out in written
        assert out.exists()
        assert out.read_text().startswith(">contig1")
        # The intermediate .gz is cleaned up.
        assert not (dl / "GCF_000005845.2.fna.gz").exists()

    def test_skips_assembly_without_ftp_path(self, gatherer, tmp_path):
        written = gatherer.download_assemblies(
            [{"accession": "GCA_x", "ftp_path": ""}], tmp_path / "d"
        )
        assert written == []

    def test_download_failure_is_logged_not_raised(
        self, gatherer, monkeypatch, tmp_path
    ):
        import taxon_assembly_gatherer as tag_mod

        def boom(url, dest):
            raise OSError("network down")

        monkeypatch.setattr(tag_mod.urllib.request, "urlretrieve", boom)
        written = gatherer.download_assemblies(
            [{"accession": "GCF_1.1", "ftp_path": "https://x/y/GCF_1.1_ASM"}],
            tmp_path / "d",
        )
        assert written == []
