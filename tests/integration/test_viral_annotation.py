"""
Real end-to-end test of --annotator viral (VADR -> VIGOR4 fallback).

Fixture genomes (TESTER/viral_test/, all real NCBI records, verified by header
before being committed):
  - 6 Norovirus genomes (Caliciviridae) -- classified by the real NCBI
    "calici" VADR model library (5 PASS, 1 FAIL-but-classified on a partial
    genome; VADR classification only requires a model1 assignment, not a
    passing QC status -- see viral_annotate.is_classified).
  - 2 MERS-CoV genomes (Coronaviridae) -- a real but UNRELATED family to the
    calici library, so VADR's classification scan assigns no model at all
    (model1="-", confirmed via manual `v-annotate.pl --mdir ... --mkey calici`
    runs: "no significant similarity detected", ~6s, no OOM risk). This is
    what correctly triggers viral_annotate.py's VIGOR4 fallback -- unlike a
    same-family-but-wrong-genus genome (e.g. Sapovirus against a
    Norovirus-only sub-library), which VADR still "classifies" with a
    garbage low-score/low-coverage model rather than refusing to assign one.
    VIGOR4 + VIGOR_DB's merbe_db was confirmed (manual run) to recover 11 real
    MERS genes (ORF1a, ORF1b, S, ORF3, ORF4a/4b/5, E, M, N, ORF8b) from these
    genomes -- a genuine multi-gene annotation, not the single-polyprotein
    case that broke OrthoFinder for flavivirus/SARS-CoV-2 fixtures earlier.
  - 1 additional Norovirus genome (genomes_releaf/) added via ReLeaf after
    the initial OrthoPhyl.sh build.

Both VADR (classification + coverage determination) and VIGOR4 (exonerate
alignment across every candidate reference model) are genuinely slow on real
genomes -- the MERS/merbe_db fallback alone took ~18 minutes per genome in
manual verification. This test is correspondingly slow (expect it to run
for upwards of 45-60 minutes) -- that cost is inherent to the real tools,
not this test's fixture design.

Requires (gated via skip, no network/download fallback):
  - the `viral_annot` conda env (v-annotate.pl, vigor4, diamond) -- see
    orthophyl_viral_env.yml / control_file.paths' viral_env_run.
  - ORTHOPHYL_VADR_MODELS_DIR pointing at a real NCBI VADR Caliciviridae
    model package (e.g. the "calici" key from
    https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/caliciviridae/).
  - ORTHOPHYL_VIGOR_DB_DIR pointing at VIGOR_DB's Reference_DBs/ (must
    contain merbe_db).
"""

import csv
import os
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration


def _conda_run_ok(cmd):
    """True if `conda run -n viral_annot <cmd> -h`-style invocation resolves."""
    viral_env_run = os.environ.get("ORTHOPHYL_VIRAL_ENV_RUN", "conda run -n viral_annot")
    try:
        result = subprocess.run(
            viral_env_run.split() + ["which", cmd],
            capture_output=True, text=True, timeout=30,
        )
        return result.returncode == 0 and result.stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        return False


@pytest.fixture(scope="session")
def viral_env_run_prefix():
    return os.environ.get("ORTHOPHYL_VIRAL_ENV_RUN", "conda run -n viral_annot")


@pytest.fixture(scope="session")
def require_viral_tools(viral_env_run_prefix):
    missing = [t for t in ("v-annotate.pl", "vigor4", "diamond") if not _conda_run_ok(t)]
    if missing:
        pytest.skip(
            f"--annotator viral integration test requires {missing} resolvable via "
            f"'{viral_env_run_prefix}' (the viral_annot conda env). See orthophyl_viral_env.yml."
        )


@pytest.fixture(scope="session")
def vadr_models_dir():
    path = os.environ.get("ORTHOPHYL_VADR_MODELS_DIR")
    if not path or not Path(path).is_dir():
        pytest.skip(
            "Set ORTHOPHYL_VADR_MODELS_DIR to a real NCBI VADR Caliciviridae model "
            "package directory (containing calici.cm/calici.minfo) to run this test. "
            "Not bundled/downloaded automatically -- see "
            "https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/caliciviridae/"
        )
    return Path(path)


@pytest.fixture(scope="session")
def vigor_db_dir():
    path = os.environ.get("ORTHOPHYL_VIGOR_DB_DIR")
    if not path or not Path(path).is_dir():
        pytest.skip(
            "Set ORTHOPHYL_VIGOR_DB_DIR to VIGOR_DB's Reference_DBs/ directory "
            "(must contain merbe_db) to run this test. See orthophyl_viral_env.yml."
        )
    if not (Path(path) / "merbe_db").exists():
        pytest.skip(f"ORTHOPHYL_VIGOR_DB_DIR={path} has no merbe_db -- required by this fixture set.")
    return Path(path)


@pytest.fixture(scope="session")
def viral_test_genomes_dir(project_root) -> Path:
    d = project_root / "TESTER" / "viral_test" / "genomes"
    if not d.exists():
        pytest.fail(f"Viral test genomes not found at {d}")
    return d


@pytest.fixture(scope="session")
def viral_test_releaf_genome(project_root) -> Path:
    d = project_root / "TESTER" / "viral_test" / "genomes_releaf"
    if not d.exists():
        pytest.fail(f"Viral ReLeaf test genome not found at {d}")
    return d


_VADR_PRIMARY_STEMS = {"NC_001959", "NC_029645", "NC_029646", "KC576915", "JX846925", "GU445325"}
_VIGOR4_FALLBACK_STEMS = {"NC_019843", "KT029139"}
_RELEAF_STEM = "LC952974"


@pytest.fixture(scope="session")
def viral_orthophyl_run(
    require_viral_tools,
    tmp_path_factory,
    project_root,
    orthophyl_script,
    viral_test_genomes_dir,
    vadr_models_dir,
    vigor_db_dir,
    control_file,
    viral_env_run_prefix,
):
    """Run OrthoPhyl.sh --annotator viral once for the whole module."""
    output_dir = tmp_path_factory.mktemp("viral_orthophyl")

    cmd = [
        str(orthophyl_script),
        "-g", str(viral_test_genomes_dir),
        "-s", str(output_dir),
        "-c", str(control_file),
        "-o", "CDS",
        "-p", "fasttree",
        "-t", "4",
        "--annotator", "viral",
        "--vadr-models", f"{vadr_models_dir}:calici",
        "--vigor-db-dir", str(vigor_db_dir),
        "--vigor-db", "auto",
    ]

    env = os.environ.copy()
    env["viral_env_run"] = viral_env_run_prefix

    print(f"\n{'='*80}\nRunning OrthoPhyl.sh --annotator viral (session fixture)\n"
          f"Output: {output_dir}\nCommand: {' '.join(cmd)}\n{'='*80}\n")

    result = subprocess.run(cmd, cwd=project_root, env=env, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(
            f"OrthoPhyl.sh --annotator viral failed with exit code {result.returncode}\n\n"
            f"STDOUT:\n{result.stdout[-8000:]}\n\nSTDERR:\n{result.stderr[-8000:]}"
        )

    return output_dir


def _read_annotation_sources(output_dir) -> dict:
    tsv_path = output_dir / "annots_details" / "annotation_sources.tsv"
    if not tsv_path.exists():
        pytest.fail(f"annotation_sources.tsv not found at {tsv_path}")
    rows = {}
    with open(tsv_path) as fh:
        for line in csv.reader(fh, delimiter="\t"):
            if not line:
                continue
            rows[line[0]] = line
    return rows


class TestViralAnnotationOrthoPhyl:
    def test_annotation_sources_rows_present(self, viral_orthophyl_run):
        rows = _read_annotation_sources(viral_orthophyl_run)
        expected = _VADR_PRIMARY_STEMS | _VIGOR4_FALLBACK_STEMS
        missing = expected - set(rows)
        assert not missing, f"Missing annotation_sources.tsv rows for: {missing}"

    def test_vadr_primary_genomes_classified_via_vadr(self, viral_orthophyl_run):
        rows = _read_annotation_sources(viral_orthophyl_run)
        for stem in _VADR_PRIMARY_STEMS:
            method = rows[stem][1]
            assert method == "vadr", (
                f"{stem}: expected VADR classification (real Norovirus genome, "
                f"covered by the calici library), got method={method!r}"
            )

    def test_vigor4_fallback_genomes_classified_via_vigor4(self, viral_orthophyl_run):
        rows = _read_annotation_sources(viral_orthophyl_run)
        for stem in _VIGOR4_FALLBACK_STEMS:
            method = rows[stem][1]
            assert method == "vigor4", (
                f"{stem}: expected VIGOR4 fallback (real MERS-CoV genome, unrelated "
                f"to the calici-only VADR library so model1 should be '-'), "
                f"got method={method!r}"
            )
            library_or_db = rows[stem][2]
            assert library_or_db == "merbe_db", (
                f"{stem}: expected VIGOR_DB auto-pick to land on merbe_db, got {library_or_db!r}"
            )

    def test_trans_and_prots_headers_match(self, viral_orthophyl_run):
        trans_dir = viral_orthophyl_run / "annots_nucls"
        prots_dir = viral_orthophyl_run / "annots_prots"
        for stem in _VADR_PRIMARY_STEMS | _VIGOR4_FALLBACK_STEMS:
            trans_fa = trans_dir / f"{stem}.fna"
            prots_fa = prots_dir / f"{stem}.faa"
            assert trans_fa.exists(), f"{trans_fa} missing"
            assert prots_fa.exists(), f"{prots_fa} missing"

            def headers(path):
                return {
                    line[1:].split()[0] for line in path.read_text().splitlines()
                    if line.startswith(">")
                }

            trans_ids, prots_ids = headers(trans_fa), headers(prots_fa)
            assert trans_ids, f"{trans_fa} has no records"
            assert trans_ids == prots_ids, (
                f"{stem}: CDS/protein header mismatch: "
                f"trans-only={trans_ids - prots_ids}, prots-only={prots_ids - trans_ids}"
            )

    def test_vigor4_fallback_genomes_have_multiple_genes(self, viral_orthophyl_run):
        """The whole point of picking MERS-CoV (not a single-polyprotein virus) as the
        fallback fixture: it must yield >1 real gene per genome, or OrthoFinder's
        Stats() crashes with StopIteration (confirmed with flavivirus/SARS-CoV-2
        single-ORF fixtures in earlier manual testing)."""
        prots_dir = viral_orthophyl_run / "annots_prots"
        for stem in _VIGOR4_FALLBACK_STEMS:
            prots_fa = prots_dir / f"{stem}.faa"
            n_records = prots_fa.read_text().count(">")
            assert n_records > 1, f"{stem}: expected multiple VIGOR4-called genes, got {n_records}"

    def test_final_tree_exists(self, viral_orthophyl_run):
        tree_dir = viral_orthophyl_run / "phylo_current" / "FINAL_SPECIES_TREES"
        assert tree_dir.exists(), f"{tree_dir} not found"
        trees = list(tree_dir.glob("*.nwk")) + list(tree_dir.glob("*.tree")) + list(tree_dir.glob("*.treefile"))
        assert trees, f"No tree files found in {tree_dir}"


class TestViralAnnotationReLeaf:
    """ReLeaf.sh has no separate output flag -- it extends the original OrthoPhyl
    run's `-s` storage dir in place, writing new work under `<storage>/ReLeaf_dir/`
    (mirrors tests/integration/test_releaf_fasttest.py's conventions)."""

    @pytest.fixture(scope="class")
    def viral_releaf_run(
        self,
        viral_orthophyl_run,
        require_viral_tools,
        project_root,
        releaf_script,
        viral_test_releaf_genome,
        vadr_models_dir,
        vigor_db_dir,
        control_file,
        viral_env_run_prefix,
    ):
        """Add one more Norovirus genome to the OrthoPhyl run via ReLeaf (in place)."""
        cmd = [
            str(releaf_script),
            "-g", str(viral_test_releaf_genome),
            "-s", str(viral_orthophyl_run),
            "-c", str(control_file),
            "-t", "4",
            "-p", "fasttree",
            "-o", "CDS",
            "--annotator", "viral",
            "--vadr-models", f"{vadr_models_dir}:calici",
            "--vigor-db-dir", str(vigor_db_dir),
            "--vigor-db", "auto",
        ]

        env = os.environ.copy()
        env["viral_env_run"] = viral_env_run_prefix

        print(f"\n{'='*80}\nRunning ReLeaf.sh --annotator viral\n"
              f"Command: {' '.join(cmd)}\n{'='*80}\n")

        result = subprocess.run(cmd, cwd=project_root, env=env, capture_output=True, text=True)

        if result.returncode != 0:
            pytest.fail(
                f"ReLeaf.sh --annotator viral failed with exit code {result.returncode}\n\n"
                f"STDOUT:\n{result.stdout[-8000:]}\n\nSTDERR:\n{result.stderr[-8000:]}"
            )

        return viral_orthophyl_run

    def test_releaf_genome_annotated_via_vadr(self, viral_releaf_run):
        rows = _read_annotation_sources(viral_releaf_run)
        assert _RELEAF_STEM in rows, f"Missing annotation_sources.tsv row for {_RELEAF_STEM}"
        assert rows[_RELEAF_STEM][1] == "vadr", (
            f"{_RELEAF_STEM}: expected VADR classification (real Norovirus genome), "
            f"got method={rows[_RELEAF_STEM][1]!r}"
        )

    def test_releaf_tree_includes_new_genome(self, viral_releaf_run):
        releaf_dir = viral_releaf_run / "ReLeaf_dir"
        assert releaf_dir.exists(), f"ReLeaf_dir not created: {releaf_dir}"

        tree_candidates = [
            releaf_dir / "phylogeny_with_new_genomes.nwk",
            releaf_dir / "new_trees" / "phylogeny.nwk",
            releaf_dir / "new_trees" / "phylogeny_with_new_genomes.nwk",
        ]
        new_trees_dir = releaf_dir / "new_trees"
        if new_trees_dir.exists():
            tree_candidates.extend(new_trees_dir.glob("*.tree"))
            tree_candidates.extend(new_trees_dir.glob("*.treefile"))

        tree_path = next((p for p in tree_candidates if p.exists()), None)
        assert tree_path, f"No ReLeaf output tree found. Checked: {[str(p) for p in tree_candidates]}"

        tree_text = tree_path.read_text()
        assert _RELEAF_STEM in tree_text, f"{_RELEAF_STEM} not found as a leaf in {tree_path}"
