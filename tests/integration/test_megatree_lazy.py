"""Integration tests for --megatree --megatree-lazy end-to-end, on real chloroplast data.

Reproduces (and pins as a regression test) a scenario manually verified during
development: a novel taxon partitions into subclades, one subclade builds
eagerly (holds a query genome), the others register lazily (built=false). A
later query MASH-matches an unbuilt subclade, triggering an on-demand build,
then ReLeaf onto it.

This specifically guards the fix for a real bug found during that manual run:
_run_orthophyl never passed OrthoPhyl.sh's -n/--num_OF_annots, so any wrapper
build with <= ANI_shortlist (default 20) genomes skipped the MASH-shortlist
path and never produced OG_alignmentsToHMM/hmms_final/ -- meaning ReLeaf could
never graft onto it. --ani-shortlist (passed here as 5, matching the manual
repro) forces that path regardless of input size.
"""

import json
import subprocess
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.slow
class TestMegatreeLazyEndToEnd:
    """Drives orthophyl_pipeline_wrapper.py as a subprocess (not mocked) --
    catches bugs mocked-subprocess unit tests structurally cannot."""

    @pytest.fixture(scope="class")
    def database_dir(self, tmp_path_factory, project_root):
        """A database dir seeded with one throwaway registered DB.

        Batch (--input) mode's router refuses to start with zero databases
        loaded; a real OrthoPhyl run isn't needed for that -- a lazily
        registered placeholder (via --register-only) is enough, exactly as
        used in the manual reproduction of this scenario.
        """
        db_dir = tmp_path_factory.mktemp("megatree_lazy_db")
        db_tool = project_root / "assembly_router" / "OP_database_tool.py"
        result = subprocess.run(
            [
                "python3", str(db_tool),
                "--single-clade", "DecoyTaxon",
                "d__Bacteria;p__DecoyPhylum;g__DecoyTaxon",
                "/nonexistent/unused",
                "--output-dir", str(db_dir),
                "--is-subclade", "--parent-taxon", "DecoyTaxon",
                "--register-only", "--n-genomes", "0",
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"Decoy DB registration failed:\nSTDOUT:{result.stdout}\nSTDERR:{result.stderr}")
        return db_dir

    @pytest.fixture(scope="class")
    def step1_output_dir(self, tmp_path_factory):
        return tmp_path_factory.mktemp("megatree_lazy_step1")

    @pytest.fixture(scope="class")
    def step1_result(
        self, project_root, wrapper_script, gather_script, database_dir,
        step1_output_dir, cymbidium_megatree_lazy_files,
    ):
        """Run Step 1 once per class: create the Cymbidium megatree-lazy DB set
        from a batch containing the at-creation query."""
        raw_dir = cymbidium_megatree_lazy_files['raw_dir']
        at_creation_query = cymbidium_megatree_lazy_files['at_creation_query']

        # Stage the raw set under the wrapper's expected download path so
        # --skip-download finds it (mirrors the manual reproduction).
        staged_raw = (step1_output_dir / "02_orthophyl_novel" / "downloads"
                      / "Cymbidium" / "assemblies_all.TMP")
        staged_raw.mkdir(parents=True, exist_ok=True)
        for f in raw_dir.glob("*.fna"):
            (staged_raw / f.name).write_bytes(f.read_bytes())
        query_dst = staged_raw / at_creation_query.name
        query_dst.write_bytes(at_creation_query.read_bytes())

        batch_tsv = step1_output_dir / "batch1.tsv"
        batch_tsv.write_text(
            "assembly_path\ttaxonomy\tassembly_id\n"
            f"{query_dst}\t"
            "d__Eukaryota;p__Tracheophyta;c__Liliopsida;o__Asparagales;"
            "f__Orchidaceae;g__Cymbidium;s__\tQ1_at_creation\n"
        )

        cmd = [
            "python3", str(wrapper_script),
            "--input", str(batch_tsv),
            "--database-dir", str(database_dir),
            "--output-dir", str(step1_output_dir),
            "--gather-script", str(gather_script),
            "--skip-download",
            "--megatree", "--megatree-lazy",
            "--max-tree-genomes", "15",
            "--subclade-size", "6",
            "--ani-shortlist", "5",
            "--use-bbmap",
            "--threads", "4",
        ]
        result = subprocess.run(
            cmd, cwd=project_root, capture_output=True, text=True, timeout=900)
        return result

    def test_step1_completes_successfully(self, step1_result):
        assert step1_result.returncode == 0, (
            f"Step 1 (create) failed:\n"
            f"STDOUT:\n{step1_result.stdout}\n\nSTDERR:\n{step1_result.stderr}")

    def test_partitioned_into_three_subclades(self, step1_result, database_dir):
        for name in ("Cymbidium_1", "Cymbidium_2", "Cymbidium_3"):
            assert (database_dir / f"{name}_db" / "database_config.json").exists(), (
                f"Expected subclade DB {name}_db not created.\n"
                f"STDOUT:\n{step1_result.stdout}")

    def test_query_subclade_built_with_hmms(self, step1_result, database_dir):
        """The subclade holding the at-creation query must be built AND have
        HMMs -- this is the actual regression check for the bug: without
        --ani-shortlist forcing OrthoPhyl.sh's MASH-shortlist branch, a small
        subclade tree builds fine but never produces hmms_final/."""
        config = json.loads(
            (database_dir / "Cymbidium_1_db" / "database_config.json").read_text())
        assert config["built"] is True

        hmm_dir = (database_dir / "Cymbidium_1_db" / "orthophyl_run"
                   / "OG_alignmentsToHMM" / "hmms_final")
        assert hmm_dir.exists(), (
            f"hmms_final/ not found for the built subclade -- ReLeaf can never "
            f"graft onto this DB. STDOUT:\n{step1_result.stdout}")
        assert list(hmm_dir.glob("*.hmm")), "hmms_final/ exists but is empty"

    def test_other_subclades_registered_lazily(self, step1_result, database_dir):
        for name in ("Cymbidium_2", "Cymbidium_3"):
            config = json.loads(
                (database_dir / f"{name}_db" / "database_config.json").read_text())
            assert config["built"] is False, (
                f"{name}_db should be a lazy (built=false) placeholder -- it "
                f"held no query at partition time.")

    @pytest.fixture(scope="class")
    def step2_output_dir(self, tmp_path_factory):
        return tmp_path_factory.mktemp("megatree_lazy_step2")

    @pytest.fixture(scope="class")
    def step2_result(
        self, step1_result, project_root, wrapper_script, gather_script,
        database_dir, step2_output_dir, cymbidium_megatree_lazy_files,
    ):
        """Run Step 2 once per class (depends on step1_result via fixture
        ordering): route the held-out query, triggering the on-demand build
        of whichever lazy subclade it MASH-matches, then ReLeaf onto it."""
        assert step1_result.returncode == 0, "Step 1 must succeed before Step 2 runs"

        held_out_query = cymbidium_megatree_lazy_files['held_out_query']
        batch_tsv = step2_output_dir / "batch2.tsv"
        batch_tsv.write_text(
            "assembly_path\ttaxonomy\tassembly_id\n"
            f"{held_out_query}\t"
            "d__Eukaryota;p__Tracheophyta;c__Liliopsida;o__Asparagales;"
            "f__Orchidaceae;g__Cymbidium;s__\tQ2_heldout\n"
        )

        cmd = [
            "python3", str(wrapper_script),
            "--input", str(batch_tsv),
            "--database-dir", str(database_dir),
            "--output-dir", str(step2_output_dir),
            "--gather-script", str(gather_script),
            "--ani-shortlist", "5",
            "--use-bbmap",
            "--threads", "4",
        ]
        result = subprocess.run(
            cmd, cwd=project_root, capture_output=True, text=True, timeout=900)
        return result

    def test_step2_completes_successfully(self, step2_result):
        assert step2_result.returncode == 0, (
            f"Step 2 (route + on-demand build + ReLeaf) failed:\n"
            f"STDOUT:\n{step2_result.stdout}\n\nSTDERR:\n{step2_result.stderr}")

    def test_held_out_query_routed_to_nearest_lazy_subclade(self, step2_result):
        """Manually verified via direct MASH distance: the held-out genome is
        ~0.0003 from Cymbidium_2's nearest member vs ~0.006/0.01 for the
        other two subclades -- it must route to Cymbidium_2 specifically.

        assembly_router.py's own routing-decision print (which names the
        pipeline) lands on the CHILD subprocess's stdout, but the wrapper
        redirects that child's combined output to a log file rather than
        passing it through -- check the per-subclade log the wrapper writes
        instead of step2_result.stdout/stderr.
        """
        routing_log = list(Path(step2_result.args[step2_result.args.index("--output-dir") + 1])
                           .glob("00_routing/routing_decision_*.json"))
        assert routing_log, "No routing decision JSON found"
        decision = json.loads(routing_log[0].read_text())
        assert decision["pipeline"] == "OrthoPhyl_subclade_build"
        assert decision["subclade_name"] == "Cymbidium_2"

    def test_lazy_subclade_promoted_to_built(self, step2_result, database_dir):
        config = json.loads(
            (database_dir / "Cymbidium_2_db" / "database_config.json").read_text())
        assert config["built"] is True, (
            f"Cymbidium_2_db should be promoted built=false -> true after the "
            f"on-demand build.\nSTDOUT:\n{step2_result.stdout}\n\n"
            f"STDERR:\n{step2_result.stderr}")

    def test_releaf_succeeds_onto_freshly_built_subclade(self, step2_result, database_dir):
        """The actual fix under test: before it, this step failed with
        'No HMM files found' because the freshly-built subclade had no HMMs.

        logger.info's default StreamHandler writes to stderr (see
        orthophyl_pipeline_wrapper.py's _log_command docstring), so these
        messages land on stderr, not stdout.
        """
        assert "No HMM files found" not in step2_result.stderr
        assert "ReLeaf complete for Cymbidium_2" in step2_result.stderr

        new_trees_dir = (database_dir / "Cymbidium_2_db" / "orthophyl_run"
                         / "ReLeaf_dir" / "new_trees")
        assert new_trees_dir.exists(), (
            f"ReLeaf_dir/new_trees not produced.\nSTDOUT:\n{step2_result.stdout}")
        treefiles = list(new_trees_dir.glob("*.addasm.treefile"))
        assert treefiles, "No .addasm.treefile produced by ReLeaf"

    def test_held_out_query_present_in_final_tree(self, step2_result, database_dir):
        new_trees_dir = (database_dir / "Cymbidium_2_db" / "orthophyl_run"
                         / "ReLeaf_dir" / "new_trees")
        treefiles = list(new_trees_dir.glob("*.addasm.treefile"))
        assert treefiles, "No .addasm.treefile produced by ReLeaf"
        tree_text = treefiles[0].read_text()
        assert "Q2_heldout" in tree_text, (
            f"Query 'Q2_heldout' not found in final tree: {treefiles[0]}")
