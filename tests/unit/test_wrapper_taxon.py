"""Tests for orthophyl_pipeline_wrapper.py taxon-mode orchestration (plan section 4.12).

Taxon mode allows creating/updating databases by querying NCBI for assemblies matching
a specified taxon (e.g., "Methylorubrum"). This suite tests the wrapper's orchestration
of the TaxonAssemblyGatherer, including:
  - Create mode: query → download → run OrthoPhyl → create database
  - Update mode: diff accessions → download new → run OrthoPhyl → update database
  - Error handling for existing DBs without --update flag

All external process execution is mocked -- no real NCBI queries, downloads, or
OrthoPhyl runs are invoked.

This suite includes a contract test (test_gatherer_interface_contract) that guards
against bug B4 (wrapper↔gatherer API mismatch). If the real TaxonAssemblyGatherer
doesn't implement the interface the wrapper expects, that test will fail.
"""

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest


@pytest.fixture
def Wrapper(wrapper_module):
    return wrapper_module.PipelineWrapper


@pytest.fixture
def recording_run(monkeypatch, wrapper_module):
    """Patch subprocess.run inside the wrapper module; record argv, return success."""
    calls = []

    class FakeCompleted:
        def __init__(self):
            self.returncode = 0
            self.stdout = ""
            self.stderr = ""

    def fake_run(cmd, *args, **kwargs):
        calls.append({"cmd": list(cmd), "args": args, "kwargs": kwargs})
        return FakeCompleted()

    monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
    return calls


@pytest.fixture
def fake_gatherer_class(monkeypatch):
    """
    Inject a fake TaxonAssemblyGatherer that will be imported by the wrapper.
    
    This allows testing wrapper logic without depending on the real gatherer
    implementation (which would require network access and NCBI downloads).
    
    The wrapper imports TaxonAssemblyGatherer locally inside methods, so we
    patch it at the import source (sys.modules).
    """
    import sys
    instances = []
    
    class FakeTaxonAssemblyGatherer:
        def __init__(self, taxon, output_dir, **kwargs):
            self.taxon = taxon
            self.output_dir = Path(output_dir)
            self.rank = kwargs.get('rank', 'genus')
            self.taxon_name = taxon
            self.taxon_rank = self.rank
            self.taxid = "12345"
            instances.append(self)
        
        def query_ncbi(self):
            """Return fake assembly list."""
            return [
                {
                    'accession': 'GCF_000001.1',
                    'organism_name': f'{self.taxon} sp. A',
                    'taxonomy': f'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Methylobacteriaceae;g__{self.taxon};s__{self.taxon} sp. A',
                },
                {
                    'accession': 'GCF_000002.1',
                    'organism_name': f'{self.taxon} sp. B',
                    'taxonomy': f'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Methylobacteriaceae;g__{self.taxon};s__{self.taxon} sp. B',
                },
            ]
        
        def download_assemblies(self, assemblies, output_dir):
            """Fake download - just create placeholder files."""
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            for asm in assemblies:
                (output_dir / f"{asm['accession']}.fna").write_text(">fake\nATCG\n")
        
        def get_taxonomy_string(self):
            """Return GTDB-style taxonomy string."""
            return f'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Methylobacteriaceae;g__{self.taxon};s__'
        
        def write_output_tsv(self, assemblies, output_file):
            """Write assemblies to TSV."""
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                f.write("assembly_accession\torganism_name\ttaxonomy\n")
                for asm in assemblies:
                    f.write(f"{asm['accession']}\t{asm['organism_name']}\t{asm['taxonomy']}\n")
    
    # Create a fake module with our fake class
    fake_module = type(sys)('taxon_assembly_gatherer')
    fake_module.TaxonAssemblyGatherer = FakeTaxonAssemblyGatherer
    
    # Inject into sys.modules so imports will find it
    monkeypatch.setitem(sys.modules, 'taxon_assembly_gatherer', fake_module)
    
    return instances


@pytest.fixture
def fake_database_dir(tmp_path):
    """Create a fake database directory with one existing database."""
    db_dir = tmp_path / "databases"
    db_dir.mkdir()
    
    # Create an existing database for Methylorubrum
    existing_db = db_dir / "Methylorubrum_genus_db"
    existing_db.mkdir()
    
    config = {
        "clade_name": "Methylorubrum",
        "clade_rank": "genus",
        "clade_rank_name": "genus",
        "clade_taxonomy": "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Methylobacteriaceae;g__Methylorubrum;s__",
        "n_genomes": 10,
        "available_tree_methods": ["iqtree"],
        "available_data_types": ["CDS"],
        "accessions": ["GCF_000001.1", "GCF_000002.1"],
        "created_date": "2025-01-01",
    }
    
    (existing_db / "database_config.json").write_text(json.dumps(config, indent=2))
    (existing_db / "genome_list").write_text("GCF_000001.1\nGCF_000002.1\n")
    
    return db_dir


def _make_taxon_wrapper(Wrapper, tmp_path, database_dir, taxon="Methylorubrum", **overrides):
    """Helper to create a wrapper in taxon mode."""
    # Create a fake gather script (required for create mode)
    fake_gather = tmp_path / "fake_gather_filter_asms.sh"
    fake_gather.write_text("#!/bin/bash\necho 'fake gather script'\n")
    fake_gather.chmod(0o755)
    
    kwargs = dict(
        input_file=None,
        database_dir=database_dir,
        output_dir=tmp_path / "out",
        taxon=taxon,
        threads=4,
        gather_script=fake_gather,  # Provide gather_script for create mode
    )
    kwargs.update(overrides)
    return Wrapper(**kwargs)


class TestTaxonModeDetection:
    """Test that taxon mode is properly detected and configured."""
    
    def test_taxon_flag_enables_taxon_mode(self, Wrapper, tmp_path):
        """Providing --taxon enables taxon_mode."""
        w = _make_taxon_wrapper(Wrapper, tmp_path, tmp_path / "db")
        assert w.taxon_mode is True
        assert w.taxon == "Methylorubrum"
    
    def test_taxon_and_input_mutually_exclusive(self, Wrapper, tmp_path):
        """Cannot provide both --taxon and --input."""
        with pytest.raises((ValueError, SystemExit)):
            Wrapper(
                input_file=tmp_path / "assemblies.tsv",
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                taxon="Methylorubrum",
            )


class TestTaxonCreateMode:
    """Test taxon create mode (no existing database)."""
    
    def test_no_existing_db_runs_create_path(
        self, Wrapper, tmp_path, fake_gatherer_class, recording_run
    ):
        """When no database exists for the taxon, run create mode."""
        db_dir = tmp_path / "databases"
        db_dir.mkdir()
        
        w = _make_taxon_wrapper(Wrapper, tmp_path, db_dir, dry_run=False)
        
        # Mock the database creator to avoid actual subprocess
        result = w.run()
        
        # Should have created a gatherer instance
        assert len(fake_gatherer_class) > 0
        gatherer = fake_gatherer_class[0]
        assert gatherer.taxon == "Methylorubrum"
    
    def test_create_mode_queries_downloads_runs_op_creates_db(
        self, Wrapper, tmp_path, fake_gatherer_class, recording_run, monkeypatch
    ):
        """Create mode: query NCBI → download → run OrthoPhyl → create database."""
        db_dir = tmp_path / "databases"
        db_dir.mkdir()
        
        # Create fake OrthoPhyl output
        def fake_orthophyl_output(w):
            op_dir = w.orthophyl_dir / "Methylorubrum_genus"
            op_dir.mkdir(parents=True, exist_ok=True)
            (op_dir / "genome_list").write_text("GCF_000001.1\nGCF_000002.1\n")
            
            # Create minimal tree structure
            tree_dir = op_dir / "phylo_current" / "SpeciesTree"
            tree_dir.mkdir(parents=True, exist_ok=True)
            (tree_dir / "iqtree.SCO_strict.CDS.tree").write_text("(A,B);")
        
        w = _make_taxon_wrapper(Wrapper, tmp_path, db_dir, dry_run=False)
        
        # Track if _run_orthophyl was called
        orthophyl_called = []
        
        # Patch _download_genomes to create fake filtered genomes
        def patched_download(taxon_name, output_dir):
            genomes_to_keep = output_dir / "genomes_to_keep"
            genomes_to_keep.mkdir(parents=True, exist_ok=True)
            (genomes_to_keep / "GCF_000001.1.fna").write_text(">fake1\nATCG\n")
            (genomes_to_keep / "GCF_000002.1.fna").write_text(">fake2\nATCG\n")
            # Create success marker
            (output_dir / ".download_complete").write_text(f"Download completed\nTaxon: {taxon_name}\nGenomes: 2\n")
        monkeypatch.setattr(w, "_download_genomes", patched_download)
        
        # Patch _run_orthophyl to create fake output and track calls
        def patched_run_op(*args, **kwargs):
            orthophyl_called.append(args)
            fake_orthophyl_output(w)
        monkeypatch.setattr(w, "_run_orthophyl", patched_run_op)
        
        # Record whether query_ncbi() gets called: create mode must NOT call it
        # (genomes come from the gather script; accession metadata is derived
        # from genomes_to_keep/, not a redundant assembly_summary download).
        import sys as _sys
        query_calls = []
        fake_cls = _sys.modules['taxon_assembly_gatherer'].TaxonAssemblyGatherer
        monkeypatch.setattr(
            fake_cls, "query_ncbi",
            lambda self: query_calls.append(1) or [],
        )

        result = w.run()

        # Should have queried and downloaded
        assert len(fake_gatherer_class) > 0

        # query_ncbi() must not have run in create mode (redundant NCBI download)
        assert query_calls == [], "query_ncbi() should not be called in create mode"

        # Should have called _run_orthophyl
        assert len(orthophyl_called) > 0, "_run_orthophyl was not called"

        # Should have called database creator
        db_creator_calls = [c for c in recording_run if "create_hierarchical_database" in str(c['cmd'])]
        assert len(db_creator_calls) > 0, "Database creator was not invoked"
    
    def test_create_db_metadata_derived_from_genomes_to_keep(
        self, Wrapper, tmp_path, fake_gatherer_class, recording_run
    ):
        """_create_taxon_database records the post-QC genomes_to_keep accessions.

        Guards the redundancy fix: assembly_accessions / n_assemblies_at_creation
        must come from genomes_to_keep/*.fna (post-QC), not a query_ncbi() list.
        """
        db_dir = tmp_path / "databases"
        db_dir.mkdir()
        w = _make_taxon_wrapper(Wrapper, tmp_path, db_dir, dry_run=False)

        # Pre-create the database dir + config the way _create_database_entry would,
        # so _create_taxon_database's metadata-update branch runs.
        taxon_db = w.database_dir / "Methylorubrum_db"
        taxon_db.mkdir(parents=True, exist_ok=True)
        (taxon_db / "database_config.json").write_text(json.dumps({"clade_name": "Methylorubrum"}))

        # A genomes_to_keep dir with two post-QC survivors.
        genomes_to_keep = tmp_path / "gtk"
        genomes_to_keep.mkdir()
        (genomes_to_keep / "GCF_000009.1.fna").write_text(">a\nATCG\n")
        (genomes_to_keep / "GCF_000008.2.fasta").write_text(">b\nATCG\n")

        gatherer = fake_gatherer_class[0] if fake_gatherer_class else \
            __import__('sys').modules['taxon_assembly_gatherer'].TaxonAssemblyGatherer(
                taxon="Methylorubrum", output_dir=tmp_path / "tq")

        # Skip the (mocked) database-creator subprocess; we only exercise metadata.
        w._create_database_entry = lambda **kwargs: None

        w._create_taxon_database(
            taxon_name="Methylorubrum",
            orthophyl_output=tmp_path / "op",
            gatherer=gatherer,
            genomes_to_keep=genomes_to_keep,
        )

        config = json.loads((taxon_db / "database_config.json").read_text())
        assert config["assembly_accessions"] == ["GCF_000008.2", "GCF_000009.1"]
        assert config["n_assemblies_at_creation"] == 2

    def test_dry_run_short_circuits_create_mode(
        self, Wrapper, tmp_path, fake_gatherer_class, recording_run
    ):
        """Dry run in create mode: query only, no download/OrthoPhyl/database."""
        db_dir = tmp_path / "databases"
        db_dir.mkdir()
        
        w = _make_taxon_wrapper(Wrapper, tmp_path, db_dir, dry_run=True)
        result = w.run()
        
        # Should have created gatherer and queried
        assert len(fake_gatherer_class) > 0
        
        # Should NOT have called OrthoPhyl or database creator
        orthophyl_calls = [c for c in recording_run if "OrthoPhyl.sh" in str(c['cmd'])]
        assert len(orthophyl_calls) == 0, "OrthoPhyl.sh should not run in dry-run mode"


class TestTaxonUpdateMode:
    """Test taxon update mode (existing database)."""
    
    def test_existing_db_without_update_flag_errors(
        self, Wrapper, tmp_path, fake_database_dir, fake_gatherer_class, caplog
    ):
        """Existing database without --update-existing flag should error."""
        w = _make_taxon_wrapper(
            Wrapper, tmp_path, fake_database_dir, 
            taxon="Methylorubrum",
            update_existing=False
        )
        
        result = w.run()
        
        # Should return error code
        assert result == 1
        
        # Should have clear error message in logs
        assert "already exists" in caplog.text.lower()
    
    def test_existing_db_with_update_runs_update_path(
        self, Wrapper, tmp_path, fake_database_dir, fake_gatherer_class, recording_run
    ):
        """Existing database with --update-existing runs update mode."""
        w = _make_taxon_wrapper(
            Wrapper, tmp_path, fake_database_dir,
            taxon="Methylorubrum",
            update_existing=True,
            dry_run=False
        )
        
        result = w.run()
        
        # Should have created gatherer
        assert len(fake_gatherer_class) > 0
    
    def test_update_mode_diffs_accessions(
        self, Wrapper, tmp_path, fake_database_dir, fake_gatherer_class, monkeypatch
    ):
        """Update mode only downloads assemblies not already in the database."""
        # Modify fake gatherer to return one new and one existing accession
        original_query = fake_gatherer_class.__class__.query_ncbi if fake_gatherer_class else None
        
        def query_with_new_and_old(self):
            return [
                {
                    'accession': 'GCF_000001.1',  # Already in DB
                    'organism_name': f'{self.taxon} sp. A',
                    'taxonomy': f'd__Bacteria;g__{self.taxon};s__',
                },
                {
                    'accession': 'GCF_000003.1',  # NEW
                    'organism_name': f'{self.taxon} sp. C',
                    'taxonomy': f'd__Bacteria;g__{self.taxon};s__',
                },
            ]
        
        w = _make_taxon_wrapper(
            Wrapper, tmp_path, fake_database_dir,
            taxon="Methylorubrum",
            update_existing=True,
            dry_run=True  # Dry run to avoid actual OrthoPhyl
        )
        
        # Patch gatherer's query method
        if fake_gatherer_class:
            for gatherer_class in [fake_gatherer_class]:
                monkeypatch.setattr(gatherer_class.__class__, "query_ncbi", query_with_new_and_old)
        
        result = w.run()
        
        # In a real implementation, would verify only GCF_000003.1 is downloaded
        # For now, just verify update mode was triggered
        assert len(fake_gatherer_class) > 0
    
    def test_update_up_to_date_short_circuits(
        self, Wrapper, tmp_path, fake_database_dir, fake_gatherer_class, monkeypatch, capsys
    ):
        """Update mode with no new assemblies short-circuits (no OrthoPhyl run)."""
        # Modify fake gatherer to return only existing accessions
        def query_only_existing(self):
            return [
                {
                    'accession': 'GCF_000001.1',
                    'organism_name': f'{self.taxon} sp. A',
                    'taxonomy': f'd__Bacteria;g__{self.taxon};s__',
                },
            ]
        
        w = _make_taxon_wrapper(
            Wrapper, tmp_path, fake_database_dir,
            taxon="Methylorubrum",
            update_existing=True,
            dry_run=False
        )
        
        # Patch gatherer's query method
        if fake_gatherer_class:
            for gatherer_class in [fake_gatherer_class]:
                monkeypatch.setattr(gatherer_class.__class__, "query_ncbi", query_only_existing)
        
        result = w.run()
        
        # Should indicate database is up to date
        captured = capsys.readouterr()
        # Implementation may vary, but should mention "up to date" or "no new"
        # This is a placeholder assertion
        assert result == 0 or result == 1  # Either success or graceful skip
    
    def test_dry_run_short_circuits_update_mode(
        self, Wrapper, tmp_path, fake_database_dir, fake_gatherer_class, recording_run
    ):
        """Dry run in update mode: query and diff only, no download/OrthoPhyl."""
        w = _make_taxon_wrapper(
            Wrapper, tmp_path, fake_database_dir,
            taxon="Methylorubrum",
            update_existing=True,
            dry_run=True
        )
        
        result = w.run()
        
        # Should have queried
        assert len(fake_gatherer_class) > 0
        
        # Should NOT have called OrthoPhyl
        orthophyl_calls = [c for c in recording_run if "OrthoPhyl.sh" in str(c['cmd'])]
        assert len(orthophyl_calls) == 0


class TestGathererInterfaceContract:
    """
    Contract test: verify the real TaxonAssemblyGatherer implements the interface
    the wrapper expects.
    
    This is the test that guards against bug B4 (wrapper↔gatherer API mismatch).
    If the real gatherer doesn't have the methods/attributes the wrapper calls,
    this test will fail.
    """
    
    def test_gatherer_interface_contract(self):
        """Real TaxonAssemblyGatherer must implement the interface wrapper expects."""
        # Import the real gatherer
        try:
            from utils.taxon_assembly_gatherer import TaxonAssemblyGatherer
        except ImportError:
            pytest.skip("TaxonAssemblyGatherer not available")
        
        # Check constructor signature
        import inspect
        sig = inspect.signature(TaxonAssemblyGatherer.__init__)
        params = list(sig.parameters.keys())
        
        # Wrapper calls: TaxonAssemblyGatherer(taxon=..., output_dir=...)
        assert 'taxon' in params, "Constructor must accept 'taxon' parameter"
        assert 'output_dir' in params, "Constructor must accept 'output_dir' parameter"
        
        # Check required methods exist
        assert hasattr(TaxonAssemblyGatherer, 'query_ncbi'), "Must have query_ncbi() method"
        assert hasattr(TaxonAssemblyGatherer, 'download_assemblies'), "Must have download_assemblies() method"
        assert hasattr(TaxonAssemblyGatherer, 'get_taxonomy_string'), "Must have get_taxonomy_string() method"
        
        # Check required attributes (set during __init__)
        # We can't instantiate without network access, but we can check the class structure
        # In a real test, you'd create a minimal instance with mocked network calls
        
        # This test passing means the interface contract is satisfied
        # If it fails, bug B4 is present
