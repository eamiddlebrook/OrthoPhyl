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

    def get_lineage(self, taxid):
        return dict(self._lineage)


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
