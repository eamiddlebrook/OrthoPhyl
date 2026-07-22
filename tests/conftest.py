"""Shared fixtures and helpers for the OrthoPhyl wrapper/router test suite.

Several target modules have dotted filenames (e.g. ``assembly_router_multi.cmd_out3.py``,
``orthophyl_pipeline_wrapper.v2.py``) which cannot be imported with a normal ``import``
statement. We load them by path with :mod:`importlib` and expose the resulting module
objects as fixtures.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

# Repo root is the parent of the tests/ directory.
REPO_ROOT = Path(__file__).resolve().parent.parent


def load_module(path, name=None):
    """Load a Python module from an arbitrary file path.

    Works for files whose names are not valid module identifiers (dots, versions).
    The module is registered in ``sys.modules`` under a sanitized name so that
    dataclasses / pickling / logging that reference ``__module__`` behave.
    """
    path = Path(path)
    if name is None:
        # Sanitize: strip suffix, replace dots so it is a legal module name.
        name = "orthophyl_test_" + path.stem.replace(".", "_")
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create import spec for {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# Module fixtures (session-scoped: loading is cheap but avoids re-exec churn). #
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def router_module():
    return load_module(
        REPO_ROOT / "assembly_router" / "assembly_router_multi.cmd_out3.py",
        name="router_cmd_out3",
    )


@pytest.fixture(scope="session")
def db_creator_module():
    return load_module(
        REPO_ROOT / "assembly_router" / "create_hierarchical_database_v2.py",
        name="create_hierarchical_database_v2",
    )


@pytest.fixture(scope="session")
def wrapper_module():
    return load_module(
        REPO_ROOT / "orthophyl_pipeline_wrapper.v2.py",
        name="orthophyl_pipeline_wrapper_v2",
    )


@pytest.fixture(scope="session")
def add_releaf_version_module():
    return load_module(
        REPO_ROOT / "assembly_router" / "add_releaf_version.py",
        name="add_releaf_version",
    )


@pytest.fixture(scope="session")
def ncbi_stats_module():
    return load_module(
        REPO_ROOT / "utils" / "ncbi_assembly_stats_v2.py",
        name="ncbi_assembly_stats_v2",
    )


@pytest.fixture(scope="session")
def taxon_gatherer_module():
    return load_module(
        REPO_ROOT / "utils" / "taxon_assembly_gatherer.py",
        name="taxon_assembly_gatherer",
    )


# --------------------------------------------------------------------------- #
# Filesystem fixtures                                                         #
# --------------------------------------------------------------------------- #

# GTDB taxonomy strings reused across suites.
GAIELLALES_TAX = (
    "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales"
)
ESCHERICHIA_TAX = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia"
)


@pytest.fixture
def make_db_dir(tmp_path):
    """Factory building a fake ``{clade}_db`` directory with a database_config.json.

    Returns a callable::

        db_path = make_db_dir(
            clade_name="Escherichia",
            clade_taxonomy=ESCHERICHIA_TAX,
            clade_rank="g",
            clade_rank_name="genus",
            n_genomes=42,
            tree_methods=["iqtree", "fasttree"],
            data_types=["CDS", "PROT"],
            parent=<dir>,   # where to create it (default: tmp_path/databases)
        )
    """
    import json

    def _make(
        clade_name,
        clade_taxonomy,
        clade_rank,
        clade_rank_name,
        n_genomes=10,
        tree_methods=("iqtree",),
        data_types=("CDS",),
        parent=None,
    ):
        base = Path(parent) if parent else (tmp_path / "databases")
        base.mkdir(parents=True, exist_ok=True)
        safe = clade_name.replace(" ", "_").replace("/", "_")
        db_dir = base / f"{safe}_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        config = {
            "created": "2025-01-01T00:00:00",
            "version": "1.0",
            "database_type": "hierarchical_taxonomy",
            "clade_name": clade_name,
            "clade_taxonomy": clade_taxonomy,
            "clade_rank": clade_rank,
            "clade_rank_name": clade_rank_name,
            "orthophyl_source": str(db_dir / "orthophyl_run"),
            "n_genomes": n_genomes,
            "has_hmms": True,
            "has_trees": True,
            "available_tree_methods": list(tree_methods),
            "available_data_types": list(data_types),
        }
        (db_dir / "database_config.json").write_text(json.dumps(config, indent=2))
        # A stand-in orthophyl_run target so path joins resolve.
        (db_dir / "orthophyl_run").mkdir(exist_ok=True)
        return db_dir

    return _make


@pytest.fixture
def orthophyl_run_skeleton(tmp_path):
    """Create a minimal *valid* OrthoPhyl output directory.

    Contains HMMs, a trimmed protein alignment dir, an IQ-TREE species tree, and a
    genome_list -- enough to satisfy ``validate_orthophyl_run``.
    """
    def _make(n_genomes=3, with_tree=True, with_hmms=True):
        run = tmp_path / "orthophyl_run"
        run.mkdir(parents=True, exist_ok=True)

        if with_hmms:
            hmm_dir = run / "OG_alignmentsToHMM" / "hmms_final"
            hmm_dir.mkdir(parents=True, exist_ok=True)
            for i in range(3):
                (hmm_dir / f"OG{i:07d}.hmm").write_text("HMMER3/f placeholder\n")

        aln_dir = run / "phylo_current" / "AlignmentsProts.trm"
        aln_dir.mkdir(parents=True, exist_ok=True)
        for i in range(3):
            (aln_dir / f"OG{i:07d}.fa").write_text(">seq1\nMAA\n")

        if with_tree:
            tree_dir = run / "FINAL_SPECIES_TREES"
            tree_dir.mkdir(parents=True, exist_ok=True)
            (tree_dir / "SCO_strict.CDS.fasttree.tree").write_text("(a,b);\n")
            (tree_dir / "SCO_strict.CDS.iqtree.treefile").write_text("(a,b,c);\n")

        genomes = [f"genome_{i}" for i in range(n_genomes)]
        (run / "genome_list").write_text(
            "# genomes\n" + "\n".join(genomes) + "\n"
        )
        return run

    return _make
