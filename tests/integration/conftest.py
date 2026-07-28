"""Pytest configuration and fixtures for integration tests."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from .helpers.tree_compare import clean_taxon_name


def _collect_basenames(directory, patterns):
    """
    Collect cleaned taxon basenames from files matching patterns.
    
    Args:
        directory: Path object to search.
        patterns: List of glob patterns (e.g., ["*.fasta", "*.fa"]).
        
    Returns:
        Set of cleaned taxon names (basenames without extensions).
    """
    names = set()
    for pattern in patterns:
        for path in directory.glob(pattern):
            names.add(clean_taxon_name(path.stem))
    return names


def pytest_configure(config):
    """
    Configure pytest to skip integration tests by default.
    
    Integration tests are opt-in via:
    - Explicit marker selection: pytest -m integration
    - Environment variable: ORTHOPHYL_RUN_INTEGRATION=1
    """
    # Check if integration tests should run
    run_integration = (
        config.getoption("-m") == "integration" or
        os.getenv("ORTHOPHYL_RUN_INTEGRATION", "").lower() in ("1", "true", "yes")
    )
    
    # If not explicitly requested, exclude integration tests
    if not run_integration:
        markexpr = config.option.markexpr
        if markexpr:
            config.option.markexpr = f"({markexpr}) and not integration"
        else:
            config.option.markexpr = "not integration"


@pytest.fixture(scope="session", autouse=True)
def require_tools():
    """
    Check that required bioinformatics tools are available.
    
    Skips all integration tests if any required tool is missing.
    This allows unit tests to run in minimal environments while
    integration tests require the full OrthoPhyl conda environment.
    """
    required_tools = [
        "iqtree",
        "orthofinder",
        "mafft",
        "hmmbuild",
        "hmmsearch",
        "trimal",
    ]
    
    missing = [tool for tool in required_tools if shutil.which(tool) is None]
    
    if missing:
        pytest.skip(
            f"Integration tests require the following tools on PATH: {', '.join(missing)}\n"
            f"Activate the OrthoPhyl conda environment: conda activate OrthoPhyl"
        )


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Return the project root directory."""
    # conftest.py is in tests/integration/, so go up two levels
    return Path(__file__).parent.parent.parent


@pytest.fixture(scope="session")
def orthophyl_script(project_root) -> Path:
    """Return the path to OrthoPhyl.sh."""
    script = project_root / "OrthoPhyl.sh"
    if not script.exists():
        pytest.fail(f"OrthoPhyl.sh not found at {script}")
    return script


@pytest.fixture(scope="session")
def releaf_script(project_root) -> Path:
    """Return the path to ReLeaf.sh."""
    script = project_root / "ReLeaf.sh"
    if not script.exists():
        pytest.fail(f"ReLeaf.sh not found at {script}")
    return script


@pytest.fixture(scope="session")
def fasttest_genomes(project_root) -> Path:
    """Return the path to the fasttest genomes directory."""
    genomes_dir = project_root / "TESTER" / "genomes_fasttest"
    if not genomes_dir.exists():
        pytest.fail(f"Fasttest genomes not found at {genomes_dir}")
    return genomes_dir


@pytest.fixture(scope="session")
def fasttest_annots_nucls(project_root) -> Path:
    """Return the path to the fasttest nucleotide annotations."""
    annots_dir = project_root / "TESTER" / "annots_nucls_fasttest"
    if not annots_dir.exists():
        pytest.fail(f"Fasttest nucleotide annotations not found at {annots_dir}")
    return annots_dir


@pytest.fixture(scope="session")
def fasttest_annots_prots(project_root) -> Path:
    """Return the path to the fasttest protein annotations."""
    annots_dir = project_root / "TESTER" / "annots_prots_fasttest"
    if not annots_dir.exists():
        pytest.fail(f"Fasttest protein annotations not found at {annots_dir}")
    return annots_dir


@pytest.fixture(scope="session")
def fasttest_addasm_genomes(project_root) -> Path:
    """Return the path to the fasttest add-assembly genomes."""
    genomes_dir = project_root / "TESTER" / "genomes_fasttest_addasm"
    if not genomes_dir.exists():
        pytest.fail(f"Fasttest add-assembly genomes not found at {genomes_dir}")
    return genomes_dir


@pytest.fixture(scope="session")
def fasttest_addasm_annots_nucls(project_root) -> Path:
    """Return the path to the fasttest add-assembly nucleotide annotations."""
    annots_dir = project_root / "TESTER" / "annots_nucls_fasttest_addasm"
    if not annots_dir.exists():
        pytest.fail(f"Fasttest add-assembly nucleotide annotations not found at {annots_dir}")
    return annots_dir


@pytest.fixture(scope="session")
def fasttest_addasm_annots_prots(project_root) -> Path:
    """Return the path to the fasttest add-assembly protein annotations."""
    annots_dir = project_root / "TESTER" / "annots_prots_fasttest_addasm"
    if not annots_dir.exists():
        pytest.fail(f"Fasttest add-assembly protein annotations not found at {annots_dir}")
    return annots_dir


@pytest.fixture(scope="session")
def reference_trees(project_root) -> Path:
    """Return the path to the reference trees directory."""
    ref_dir = project_root / "TESTER" / "REFERENCE_TESTER_TREES"
    if not ref_dir.exists():
        pytest.fail(f"Reference trees not found at {ref_dir}")
    return ref_dir


@pytest.fixture(scope="session")
def control_file(project_root) -> Path:
    """Return the path to the control file."""
    control = project_root / "control_file.user"
    if not control.exists():
        pytest.fail(f"Control file not found at {control}")
    return control


@pytest.fixture(scope="session")
def orthophyl_run(
    tmp_path_factory,
    project_root,
    orthophyl_script,
    fasttest_genomes,
    fasttest_annots_nucls,
    fasttest_annots_prots,
    control_file,
):
    """
    Run OrthoPhyl.sh once for the entire test session.
    
    This is a session-scoped fixture that runs OrthoPhyl on the fasttest
    dataset and returns the output directory. All tests that need the
    OrthoPhyl output can depend on this fixture, avoiding redundant runs.
    
    Returns:
        Path to the OrthoPhyl output directory.
    """
    # Create a temporary output directory
    output_dir = tmp_path_factory.mktemp("orthophyl_fasttest")
    
    # Build the command (mirrors test.sh)
    cmd = [
        str(orthophyl_script),
        "-s", str(output_dir),
        "-a", f"{fasttest_annots_nucls},{fasttest_annots_prots}",
        "-o", "BOTH",
        "-R", "full",
        "-g", str(fasttest_genomes),
        "-t", "4",
        "-c", str(control_file),
        "-n", "5",  # Trigger MASH shortlist path to build HMMs (9 taxa > 5)
    ]
    
    print(f"\n{'='*80}")
    print(f"Running OrthoPhyl.sh (session fixture, runs once)")
    print(f"Output: {output_dir}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*80}\n")
    
    # Run OrthoPhyl
    result = subprocess.run(
        cmd,
        cwd=project_root,
        capture_output=True,
        text=True,
    )
    
    # Check for success
    if result.returncode != 0:
        pytest.fail(
            f"OrthoPhyl.sh failed with exit code {result.returncode}\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )
    
    print(f"\n{'='*80}")
    print(f"OrthoPhyl.sh completed successfully")
    print(f"Output directory: {output_dir}")
    print(f"{'='*80}\n")
    
    return output_dir


@pytest.fixture(scope="session")
def expected_orthophyl_taxa(fasttest_genomes, fasttest_annots_prots):
    """
    Cleaned taxon names for all OrthoPhyl inputs (genomes + pre-annotated prots).
    
    Dynamically derives expected taxa from the actual input directories,
    making tests self-updating if the test data changes.
    
    Returns:
        Set of cleaned taxon names (9 taxa for fasttest).
    """
    taxa = _collect_basenames(fasttest_genomes, ["*.fasta", "*.fa"])
    taxa |= _collect_basenames(fasttest_annots_prots, ["*.faa"])
    return taxa


@pytest.fixture(scope="session")
def expected_releaf_added_taxa(fasttest_addasm_genomes, fasttest_addasm_annots_prots):
    """
    Cleaned taxon names added during ReLeaf.
    
    Returns:
        Set of cleaned taxon names (5 taxa for fasttest addasm).
    """
    taxa = _collect_basenames(fasttest_addasm_genomes, ["*.fasta", "*.fa"])
    taxa |= _collect_basenames(fasttest_addasm_annots_prots, ["*.faa"])
    return taxa


@pytest.fixture(scope="session")
def expected_releaf_total_taxa(expected_orthophyl_taxa, expected_releaf_added_taxa):
    """
    All taxa expected in the final ReLeaf tree (original + added).
    
    Returns:
        Set of cleaned taxon names (14 taxa = 9 original + 5 added).
    """
    return expected_orthophyl_taxa | expected_releaf_added_taxa
