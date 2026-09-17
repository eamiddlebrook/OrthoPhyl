"""Tests for local genome-ingest mode in orthophyl_pipeline_wrapper.py (--genome-dir).

This third mode builds a tree + database from genomes already on disk, under a
user-supplied --clade-name that may or may not resolve against the local NCBI
taxdump. Covers:
  - mode detection and pairwise mode conflicts (incl. a regression pin on the
    original --input/--taxon message)
  - taxonomy resolution precedence (--clade-taxonomy > resolved --clade-name >
    unresolvable fallback), and that resolvability actually determines
    is_within_clade routability (cross-module, via router_module)
  - staging: extension normalization, originals untouched, duplicate-stem guard
  - --skip-qc vs default QC, --use-bbmap/--must-keep passthrough
  - subsample-over-ceiling uses a distinct qc/ dir
  - genome-count guards (<4 before/after QC)
  - DB-collision pre-flight short-circuits before staging
  - --dry-run writes nothing
  - config provenance fields written by _create_local_database
  - _phase_initialization does not raise for local mode with no existing DB dir

All external process execution (gather script, OrthoPhyl, DB creator) is mocked.
"""

import json
from pathlib import Path

import pytest


@pytest.fixture
def Wrapper(wrapper_module):
    return wrapper_module.PipelineWrapper


@pytest.fixture
def recording_run(monkeypatch, wrapper_module):
    calls = []

    class FakeCompleted:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, *args, **kwargs):
        calls.append([str(x) for x in cmd])
        return FakeCompleted()

    monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)
    return calls


def _make_genome_dir(tmp_path, stems=("g0", "g1", "g2", "g3"), ext=".fna"):
    d = tmp_path / "my_genomes"
    d.mkdir(exist_ok=True)
    for s in stems:
        (d / f"{s}{ext}").write_text(">c\nACGT\n")
    return d


def _make_local_wrapper(Wrapper, tmp_path, genome_dir=None, clade_name="Blorptaxon",
                         **overrides):
    db_dir = tmp_path / "databases"
    db_dir.mkdir(exist_ok=True)
    fake_gather = tmp_path / "fake_gather.sh"
    fake_gather.write_text("#!/bin/bash\necho fake\n")
    fake_gather.chmod(0o755)
    if genome_dir is None:
        genome_dir = _make_genome_dir(tmp_path)
    kwargs = dict(
        database_dir=db_dir,
        output_dir=tmp_path / "out",
        threads=4,
        gather_script=fake_gather,
        genome_dir=genome_dir,
        clade_name=clade_name,
    )
    kwargs.update(overrides)
    w = Wrapper(**kwargs)
    w.logs_dir.mkdir(parents=True, exist_ok=True)
    w.checkpoint_dir.mkdir(parents=True, exist_ok=True)
    return w


class TestModeDetectionAndConflicts:
    def test_genome_dir_enables_local_mode(self, Wrapper, tmp_path):
        w = _make_local_wrapper(Wrapper, tmp_path)
        assert w.local_mode is True
        assert w.clade_name == "Blorptaxon"
        assert w.taxon_mode is False

    def test_input_and_taxon_regression_message(self, Wrapper, tmp_path):
        """Regression pin: --input + --taxon must still raise the ORIGINAL
        exact message (test_wrapper_taxon.py's contract), unaffected by the
        new three-way guard."""
        with pytest.raises(ValueError) as exc:
            Wrapper(
                input_file=tmp_path / "assemblies.tsv",
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                taxon="Methylorubrum",
            )
        assert str(exc.value) == (
            "Cannot specify both --input and --taxon. Use one or the other.")

    def test_input_and_genome_dir_conflict(self, Wrapper, tmp_path):
        genome_dir = _make_genome_dir(tmp_path)
        with pytest.raises(ValueError):
            Wrapper(
                input_file=tmp_path / "assemblies.tsv",
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                genome_dir=genome_dir,
                clade_name="Blorptaxon",
            )

    def test_taxon_and_genome_dir_conflict(self, Wrapper, tmp_path):
        genome_dir = _make_genome_dir(tmp_path)
        with pytest.raises(ValueError):
            Wrapper(
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                taxon="Methylorubrum",
                genome_dir=genome_dir,
                clade_name="Blorptaxon",
            )

    def test_all_three_modes_conflict(self, Wrapper, tmp_path):
        genome_dir = _make_genome_dir(tmp_path)
        with pytest.raises(ValueError):
            Wrapper(
                input_file=tmp_path / "assemblies.tsv",
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                taxon="Methylorubrum",
                genome_dir=genome_dir,
                clade_name="Blorptaxon",
            )

    def test_genome_dir_requires_clade_name(self, Wrapper, tmp_path):
        genome_dir = _make_genome_dir(tmp_path)
        with pytest.raises(ValueError, match="--clade-name is required"):
            Wrapper(
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                genome_dir=genome_dir,
            )

    def test_genome_dir_must_exist(self, Wrapper, tmp_path):
        with pytest.raises(ValueError, match="--genome-dir does not exist"):
            Wrapper(
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                genome_dir=tmp_path / "nope",
                clade_name="Blorptaxon",
            )

    def test_bad_clade_rank_rejected(self, Wrapper, tmp_path):
        genome_dir = _make_genome_dir(tmp_path)
        with pytest.raises(ValueError, match="--clade-rank must be one of"):
            Wrapper(
                database_dir=tmp_path / "db",
                output_dir=tmp_path / "out",
                genome_dir=genome_dir,
                clade_name="Blorptaxon",
                clade_rank="z",
            )


class TestPhaseInitializationLocalMode:
    def test_no_database_index_required(self, Wrapper, tmp_path):
        """Regression pin: local mode must not hit the FileNotFoundError hazard
        that non-taxon/non-local modes trigger when database_index.json is
        absent."""
        w = _make_local_wrapper(Wrapper, tmp_path)
        # database_dir exists but is otherwise empty -- no database_index.json.
        w._phase_initialization()  # must not raise
        assert w.database_dir.exists()


class TestTaxonomyResolution:
    def test_explicit_clade_taxonomy_used_verbatim(self, Wrapper, tmp_path):
        w = _make_local_wrapper(
            Wrapper, tmp_path, clade_name="Blorptaxon",
            clade_taxonomy="d__Bacteria;p__Foo;g__Blorptaxon")
        taxonomy, routable = w._resolve_local_taxonomy()
        assert taxonomy == "d__Bacteria;p__Foo;g__Blorptaxon"
        assert routable is True

    def test_clade_name_resolves_against_taxdump(self, Wrapper, tmp_path, monkeypatch):
        w = _make_local_wrapper(Wrapper, tmp_path, clade_name="Pseudomonas")

        class FakeStub:
            def _ensure_taxonomy_database(self):
                pass

        class FakeNCBITaxonomy:
            def __init__(self, taxdump_dir):
                pass

            def resolve_taxon(self, name):
                return "286" if name == "Pseudomonas" else None

            def get_lineage(self, taxid):
                return {
                    "domain": "Bacteria",
                    "phylum": "Pseudomonadota",
                    "class": "Gammaproteobacteria",
                    "order": "Pseudomonadales",
                    "family": "Pseudomonadaceae",
                    "genus": "Pseudomonas",
                }

            def get_rank(self, taxid):
                return "genus"

        def fake_render(lineage):
            order = ["domain", "phylum", "class", "order", "family", "genus", "species"]
            prefix = {"domain": "d", "phylum": "p", "class": "c", "order": "o",
                      "family": "f", "genus": "g", "species": "s"}
            parts = []
            deepest = None
            for r in order:
                if r in lineage:
                    deepest = r
            for r in order:
                parts.append(f"{prefix[r]}__{lineage.get(r, '')}")
                if r == deepest:
                    break
            return ";".join(parts)

        import types
        import sys
        fake_mod = types.ModuleType("taxon_assembly_gatherer")
        fake_mod.NCBITaxonomy = FakeNCBITaxonomy
        fake_mod.render_gtdb_lineage = fake_render
        fake_mod.TaxonAssemblyGatherer = FakeStub
        monkeypatch.setitem(sys.modules, "taxon_assembly_gatherer", fake_mod)

        taxonomy, routable = w._resolve_local_taxonomy()
        assert routable is True
        assert taxonomy == (
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas")

    def test_unresolvable_clade_name_falls_back_with_hint(self, Wrapper, tmp_path, monkeypatch, caplog):
        w = _make_local_wrapper(Wrapper, tmp_path, clade_name="Blorptaxon")

        class FakeStub:
            def _ensure_taxonomy_database(self):
                pass

        class FakeNCBITaxonomy:
            def __init__(self, taxdump_dir):
                pass

            def resolve_taxon(self, name):
                return None

        import types
        import sys
        fake_mod = types.ModuleType("taxon_assembly_gatherer")
        fake_mod.NCBITaxonomy = FakeNCBITaxonomy
        fake_mod.render_gtdb_lineage = lambda lineage: ""
        fake_mod.TaxonAssemblyGatherer = FakeStub
        monkeypatch.setitem(sys.modules, "taxon_assembly_gatherer", fake_mod)

        import logging
        with caplog.at_level(logging.WARNING):
            taxonomy, routable = w._resolve_local_taxonomy()
        assert routable is False
        assert taxonomy == "g__Blorptaxon"
        assert any("--clade-taxonomy" in rec.message for rec in caplog.records)

    def test_taxdump_failure_falls_back_rather_than_raising(self, Wrapper, tmp_path, monkeypatch):
        """Best-effort: if taxdump download/parsing fails entirely, fall back to
        the name-only taxonomy instead of failing the run."""
        w = _make_local_wrapper(Wrapper, tmp_path, clade_name="Blorptaxon")

        import types
        import sys
        fake_mod = types.ModuleType("taxon_assembly_gatherer")

        class BoomStub:
            def _ensure_taxonomy_database(self):
                raise RuntimeError("network unavailable")

        fake_mod.TaxonAssemblyGatherer = BoomStub
        fake_mod.NCBITaxonomy = object
        fake_mod.render_gtdb_lineage = lambda lineage: ""
        monkeypatch.setitem(sys.modules, "taxon_assembly_gatherer", fake_mod)

        taxonomy, routable = w._resolve_local_taxonomy()
        assert routable is False
        assert taxonomy == "g__Blorptaxon"


class TestTaxonomyRoutabilityCrossModule:
    """Pins finding 1 from the plan: a name-only taxonomy is unroutable against
    a fully-specified query, while a resolved full lineage matches correctly."""

    def test_resolved_lineage_matches_full_query(self, router_module):
        GTDBTaxonomy = router_module.GTDBTaxonomy
        resolved = GTDBTaxonomy(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas")
        query = GTDBTaxonomy(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas putida")
        assert query.is_within_clade(resolved, "g") is True

    def test_unresolvable_name_only_does_not_match_full_query(self, router_module):
        GTDBTaxonomy = router_module.GTDBTaxonomy
        unresolvable = GTDBTaxonomy("g__Blorptaxon")
        query = GTDBTaxonomy(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Pseudomonadales;f__Pseudomonadaceae;g__Blorptaxon;s__Blorptaxon sp.")
        assert query.is_within_clade(unresolvable, "g") is False


class TestStaging:
    def test_mixed_extensions_normalize_to_fna(self, Wrapper, tmp_path):
        genome_dir = tmp_path / "src"
        genome_dir.mkdir()
        (genome_dir / "a.fna").write_text(">a\nAC\n")
        (genome_dir / "b.fa").write_text(">b\nAC\n")
        (genome_dir / "c.fasta").write_text(">c\nAC\n")

        w = _make_local_wrapper(Wrapper, tmp_path, genome_dir=genome_dir)
        dest = tmp_path / "staged"
        result = w._stage_local_genomes(dest)

        assert sorted(p.name for p in result.glob("*.fna")) == [
            "a.fna", "b.fna", "c.fna"]
        # Originals untouched.
        assert (genome_dir / "a.fna").exists()
        assert (genome_dir / "b.fa").exists()
        assert (genome_dir / "c.fasta").exists()

    def test_gzip_variant_decompressed(self, Wrapper, tmp_path):
        import gzip
        genome_dir = tmp_path / "src"
        genome_dir.mkdir()
        with gzip.open(genome_dir / "d.fna.gz", "wb") as f:
            f.write(b">d\nACGT\n")

        w = _make_local_wrapper(Wrapper, tmp_path, genome_dir=genome_dir)
        dest = tmp_path / "staged"
        result = w._stage_local_genomes(dest)

        out = result / "d.fna"
        assert out.exists()
        assert out.read_bytes() == b">d\nACGT\n"
        # Original gz untouched.
        assert (genome_dir / "d.fna.gz").exists()

    def test_duplicate_stem_raises(self, Wrapper, tmp_path):
        genome_dir = tmp_path / "src"
        genome_dir.mkdir()
        (genome_dir / "foo.fna").write_text(">a\nAC\n")
        (genome_dir / "foo.fasta").write_text(">a\nAC\n")

        w = _make_local_wrapper(Wrapper, tmp_path, genome_dir=genome_dir)
        with pytest.raises(ValueError, match="normalize to the same stem"):
            w._stage_local_genomes(tmp_path / "staged")


class TestRunLocalGenomesMode:
    def _stub_common(self, w, monkeypatch, n_kept=4):
        """Stub taxonomy resolution, QC, OrthoPhyl, DB creation, tree publish."""
        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))

        def fake_qc(subclade_dir, raw, taxon_label, query_assemblies=None):
            gtk = subclade_dir / "genomes_to_keep"
            gtk.mkdir(parents=True, exist_ok=True)
            for i in range(n_kept):
                (gtk / f"k{i}.fna").write_text(">c\nAC\n")
            return gtk

        monkeypatch.setattr(w, "_qc_subclade", fake_qc)
        monkeypatch.setattr(w, "_run_orthophyl", lambda **k: None)
        db_calls = []
        monkeypatch.setattr(w, "_create_local_database",
                             lambda **k: db_calls.append(k))
        monkeypatch.setattr(w, "_locate_species_tree",
                             lambda out: out / "nonexistent.nwk")
        monkeypatch.setattr(w, "_save_final_status", lambda: None)
        return db_calls

    def test_skip_qc_issues_no_gather_subprocess(self, Wrapper, tmp_path, monkeypatch, recording_run):
        w = _make_local_wrapper(Wrapper, tmp_path, skip_qc=True)
        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))
        monkeypatch.setattr(w, "_run_orthophyl", lambda **k: None)
        db_calls = []
        monkeypatch.setattr(w, "_create_local_database",
                             lambda **k: db_calls.append(k))
        monkeypatch.setattr(w, "_locate_species_tree",
                             lambda out: out / "nonexistent.nwk")
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        rc = w._run_local_genomes_mode()
        assert rc == 0
        assert recording_run == []
        assert db_calls[0]["qc_applied"] is False

    def test_default_qc_runs_and_passes_flags(self, Wrapper, tmp_path, monkeypatch, recording_run):
        """Let the real _qc_subclade -> _download_genomes(qc_only=True) path run
        (subprocess.run itself is faked via `recording_run`), so the --qc-only
        argv actually reaches the gather script invocation."""
        w = _make_local_wrapper(Wrapper, tmp_path, use_bbmap=True, must_keep="k0")
        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))

        safe = w._default_run_name(w.clade_name)
        local_dir = w.orthophyl_dir / "local_input" / safe
        # _qc_subclade/_download_genomes validate these post-QC artifacts;
        # since subprocess.run is faked, pre-create them as the real gather
        # script would have.
        gtk = local_dir / "genomes_to_keep"
        gtk.mkdir(parents=True, exist_ok=True)
        for i in range(4):
            (gtk / f"k{i}.fna").write_text(">c\nAC\n")
        (local_dir / "assemblies_all.stats.txt").write_text("h\nk0\tx\n")

        monkeypatch.setattr(w, "_run_orthophyl", lambda **k: None)
        db_calls = []
        monkeypatch.setattr(w, "_create_local_database",
                             lambda **k: db_calls.append(k))
        monkeypatch.setattr(w, "_locate_species_tree",
                             lambda out: out / "nonexistent.nwk")
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        rc = w._run_local_genomes_mode()
        assert rc == 0
        assert db_calls[0]["qc_applied"] is True
        qc_cmd = [c for c in recording_run if "fake_gather.sh" in " ".join(c)]
        assert qc_cmd, "expected a --qc-only gather invocation"
        assert "--qc-only" in qc_cmd[0]
        assert "--use-bbmap" in qc_cmd[0]
        assert "--must-keep" in qc_cmd[0]
        assert "k0" in qc_cmd[0]

    def test_lt4_genomes_before_qc_returns_1(self, Wrapper, tmp_path, monkeypatch):
        genome_dir = _make_genome_dir(tmp_path, stems=("g0", "g1"))
        w = _make_local_wrapper(Wrapper, tmp_path, genome_dir=genome_dir)
        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))
        called = []
        monkeypatch.setattr(w, "_run_orthophyl",
                             lambda **k: called.append(k) or pytest.fail(
                                 "OrthoPhyl must not run with <4 genomes"))

        rc = w._run_local_genomes_mode()
        assert rc == 1
        assert called == []

    def test_lt4_genomes_after_qc_returns_1(self, Wrapper, tmp_path, monkeypatch):
        w = _make_local_wrapper(Wrapper, tmp_path)
        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))

        def fake_qc(subclade_dir, raw, taxon_label, query_assemblies=None):
            gtk = subclade_dir / "genomes_to_keep"
            gtk.mkdir(parents=True, exist_ok=True)
            (gtk / "k0.fna").write_text(">c\nAC\n")  # only 1 survives QC
            return gtk

        monkeypatch.setattr(w, "_qc_subclade", fake_qc)
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        called = []
        monkeypatch.setattr(w, "_run_orthophyl",
                             lambda **k: called.append(k) or pytest.fail(
                                 "OrthoPhyl must not run with <4 post-QC genomes"))

        rc = w._run_local_genomes_mode()
        assert rc == 1
        assert called == []

    def test_db_collision_returns_1_before_staging(self, Wrapper, tmp_path, monkeypatch):
        w = _make_local_wrapper(Wrapper, tmp_path, clade_name="Existing")
        existing_db = w.database_dir / "Existing_db"
        existing_db.mkdir(parents=True)
        (existing_db / "database_config.json").write_text("{}")

        stage_calls = []
        monkeypatch.setattr(w, "_stage_local_genomes",
                             lambda dest: stage_calls.append(dest) or pytest.fail(
                                 "must not stage after a DB collision"))

        rc = w._run_local_genomes_mode()
        assert rc == 1
        assert stage_calls == []

    def test_dry_run_writes_nothing(self, Wrapper, tmp_path, monkeypatch):
        w = _make_local_wrapper(Wrapper, tmp_path, dry_run=True)
        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))
        saved = []
        monkeypatch.setattr(w, "_save_final_status", lambda: saved.append(1))
        stage_calls = []
        monkeypatch.setattr(w, "_stage_local_genomes",
                             lambda dest: stage_calls.append(dest) or pytest.fail(
                                 "dry-run must not stage"))

        rc = w._run_local_genomes_mode()
        assert rc == 0
        assert stage_calls == []
        assert saved == [1]
        assert not (w.orthophyl_dir / "local_input").exists()

    def test_subsample_over_ceiling_uses_distinct_qc_dir(self, Wrapper, tmp_path, monkeypatch):
        genome_dir = _make_genome_dir(
            tmp_path, stems=[f"g{i}" for i in range(10)])
        w = _make_local_wrapper(
            Wrapper, tmp_path, genome_dir=genome_dir,
            max_tree_genomes=5, subsample_size=4)
        monkeypatch.setattr(w, "_resolve_local_taxonomy",
                             lambda: ("g__" + w.clade_name, False))

        sub_calls = []

        def fake_subsample(taxon_name, raw_dir, target, must_keep_stems=None):
            sub_calls.append({"raw_dir": raw_dir, "target": target})
            sel = w.orthophyl_dir / "subsample" / taxon_name / "selected"
            sel.mkdir(parents=True, exist_ok=True)
            for i in range(target):
                (sel / f"s{i}.fna").write_text(">c\nAC\n")
            return sel

        monkeypatch.setattr(w, "_subsample_genomes", fake_subsample)

        qc_calls = []

        def fake_qc(subclade_dir, raw, taxon_label, query_assemblies=None):
            qc_calls.append(subclade_dir)
            gtk = subclade_dir / "genomes_to_keep"
            gtk.mkdir(parents=True, exist_ok=True)
            for i in range(4):
                (gtk / f"k{i}.fna").write_text(">c\nAC\n")
            return gtk

        monkeypatch.setattr(w, "_qc_subclade", fake_qc)
        monkeypatch.setattr(w, "_run_orthophyl", lambda **k: None)
        monkeypatch.setattr(w, "_create_local_database", lambda **k: None)
        monkeypatch.setattr(w, "_locate_species_tree",
                             lambda out: out / "nonexistent.nwk")
        monkeypatch.setattr(w, "_save_final_status", lambda: None)

        rc = w._run_local_genomes_mode()
        assert rc == 0
        assert len(sub_calls) == 1
        assert sub_calls[0]["target"] == 4
        # QC ran against a distinct qc/ dir, not the full staged dir.
        safe = w._default_run_name(w.clade_name)
        local_dir = w.orthophyl_dir / "local_input" / safe
        assert qc_calls[0] == local_dir / "qc"
        assert qc_calls[0] != local_dir


class TestCreateLocalDatabase:
    def test_config_provenance_fields(self, Wrapper, tmp_path, monkeypatch, wrapper_module):
        w = _make_local_wrapper(Wrapper, tmp_path, clade_name="Blorptaxon")

        genomes_to_keep = tmp_path / "genomes_to_keep"
        genomes_to_keep.mkdir()
        for i in range(4):
            (genomes_to_keep / f"k{i}.fna").write_text(">c\nAC\n")

        orthophyl_output = tmp_path / "orthophyl_run"
        orthophyl_output.mkdir()

        # Fake the DB creator subprocess call to materialize a config file the
        # way create_hierarchical_database.py's --update path would.
        db_dir = w.database_dir / "Blorptaxon_db"

        def fake_run(cmd, *args, **kwargs):
            # Mirror create_hierarchical_database.py: taxonomy_source/qc_applied
            # in the written config reflect the --taxonomy-source/--qc-not-applied
            # argv the wrapper passed, not a fixed default.
            cmd = [str(x) for x in cmd]
            source = cmd[cmd.index("--taxonomy-source") + 1] if "--taxonomy-source" in cmd else "ncbi"
            qc_applied = "--qc-not-applied" not in cmd
            db_dir.mkdir(parents=True, exist_ok=True)
            (db_dir / "database_config.json").write_text(json.dumps({
                "clade_name": "Blorptaxon",
                "taxonomy_source": source,
                "qc_applied": qc_applied,
            }))

            class R:
                returncode = 0
            return R()

        monkeypatch.setattr(wrapper_module.subprocess, "run", fake_run)

        w._create_local_database(
            clade_name="Blorptaxon",
            orthophyl_output=orthophyl_output,
            taxonomy="g__Blorptaxon",
            taxonomy_source="user_supplied",
            qc_applied=False,
            genomes_to_keep=genomes_to_keep,
            source_genome_dir=w.genome_dir,
        )

        config = json.loads((db_dir / "database_config.json").read_text())
        assert config["taxonomy_source"] == "user_supplied"
        assert config["qc_applied"] is False
        assert config["source_taxon_name"] == "Blorptaxon"
        assert config["source_taxid"] is None
        assert config["source_rank"] is None
        assert config["source_genome_dir"] == str(w.genome_dir.resolve())
        assert sorted(config["assembly_accessions"]) == ["k0", "k1", "k2", "k3"]
        assert config["n_assemblies_at_creation"] == 4
