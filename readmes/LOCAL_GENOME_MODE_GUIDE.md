# OrthoPhyl Local Genome-Ingest Mode - Complete Guide

## Overview

**Local genome-ingest mode** (`--genome-dir`) builds a routable OrthoPhyl database from
genomes you already have on disk — unpublished isolates, a curated set, or output from
another pipeline — instead of downloading from NCBI. It reuses the same tree-building and
database-creation machinery as taxon mode (`TAXON_MODE_GUIDE.md`); the only difference is
where the genomes and the clade name come from.

1. **QC runs by default** (CheckM2, or bbmap stats with `--use-bbmap`), skippable with
   `--skip-qc` for genomes you've already quality-checked.
2. **You name the clade** with `--clade-name`. If the name resolves against NCBI's
   taxonomy, the database gets a full, routable lineage automatically. If it doesn't
   (e.g. an informal or novel name), the database still builds, but the run logs a clear
   note that **this taxonomy was not assigned by NCBI**, and the database is not matched by
   fully-specified taxonomy queries unless you also supply `--clade-taxonomy`.

## Quick Start

### Build a database from a directory of genomes

```bash
python orthophyl_pipeline_wrapper.py \
    --genome-dir /data/my_isolates/ \
    --clade-name Pseudomonas \
    --database-dir databases/ \
    --output-dir pseudomonas_local_run/ \
    --gather-script utils/gather_filter_asms.sh \
    --threads 32
```

This will:
1. Stage every `.fna`/`.fa`/`.fasta` (optionally `.gz`) file in `/data/my_isolates/` into
   a wrapper-owned working directory, normalized to `<stem>.fna`. Your original files are
   never modified or renamed.
2. Try to resolve `Pseudomonas` against the local NCBI taxdump. Since it's a real genus,
   the database gets the full lineage
   (`d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas`)
   and is routable by later queries.
3. QC the staged genomes with CheckM2 (via `gather_filter_asms.sh --qc-only`).
4. Run OrthoPhyl on the QC'd genomes to build a phylogeny.
5. Create `Pseudomonas_db` with full metadata, including `taxonomy_source: "ncbi"` and
   `qc_applied: true`.

### Skip QC for already-checked genomes, with a name NCBI won't resolve

```bash
python orthophyl_pipeline_wrapper.py \
    --genome-dir /data/qcd_isolates/ \
    --clade-name MyLabCollection \
    --skip-qc \
    --database-dir databases/ \
    --output-dir mylab_run/ \
    --threads 32
```

Since `MyLabCollection` isn't a real NCBI taxon, the run logs a warning and falls back to a
name-only taxonomy (`g__MyLabCollection` by default rank). The database still builds and is
usable directly, but won't be matched by a fully-specified `--input` taxonomy query later.
The warning includes a ready-to-paste `--clade-taxonomy` suggestion if you know the real
lineage and want the database to be routable.

### Supply the real lineage explicitly for a novel/informal name

```bash
python orthophyl_pipeline_wrapper.py \
    --genome-dir /data/novel_clade/ \
    --clade-name Blorptaxon \
    --clade-taxonomy "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Blorptaxon" \
    --database-dir databases/ \
    --output-dir blorptaxon_run/ \
    --threads 32
```

`--clade-taxonomy` is used verbatim and makes the database fully routable, even though
`Blorptaxon` itself isn't an NCBI-recognized name. `database_config.json` still records
`taxonomy_source: "user_supplied"` so the provenance is never lost, and the router's
startup log annotates the database with `[user-supplied taxonomy: not NCBI-assigned]`.

## Command-Line Arguments

#### Local genome-ingest mode arguments

- `--genome-dir DIR` - Directory of genome FASTAs already on disk (`.fna`/`.fa`/`.fasta`,
  optionally `.gz`). Mutually exclusive with `--input` and `--taxon`.
- `--clade-name NAME` - **Required** with `--genome-dir`. Names the clade and the
  resulting `<name>_db`. Auto-resolved against the local NCBI taxdump when it's a real
  taxon name; otherwise used as-is with a "not NCBI-assigned" note.
- `--clade-taxonomy STR` - Optional escape hatch: a full GTDB taxonomy string
  (`d__...;p__...;...;g__MyClade`), used verbatim. Needed only when `--clade-name` doesn't
  resolve to a known NCBI taxon and you know the real lineage.
- `--clade-rank {d,p,c,o,f,g,s}` - Rank letter at which an unresolvable `--clade-name` is
  attached. Default `g` (genus).
- `--skip-qc` - Skip CheckM2/bbmap QC on `--genome-dir` genomes (default: QC runs). Use
  only for genomes you've already quality-checked.

#### Standard arguments that also apply

- `--database-dir DIR`, `--output-dir DIR`, `--threads N`, `--gather-script`,
  `--use-bbmap`, `--low-ram`, `--must-keep`, `--keep-failing-query`, `--max-tree-genomes`,
  `--subsample-size`, `--dry-run`, `--resume`, `-v`/`-vv` — same meaning as in batch and
  taxon mode.
- `--megatree` is **not supported** with `--genome-dir` (rejected at argument-parsing
  time). Use `--subsample-size`/`--max-tree-genomes` for oversized local sets instead.

## Database Metadata

`database_config.json` for a local-genome database includes the same provenance fields as
taxon mode, plus:

```json
{
  "clade_name": "Blorptaxon",
  "clade_taxonomy": "d__Bacteria;p__Pseudomonadota;...;g__Blorptaxon",
  "taxonomy_source": "user_supplied",
  "qc_applied": true,
  "source_taxon_name": "Blorptaxon",
  "source_taxid": null,
  "source_rank": null,
  "source_genome_dir": "/data/novel_clade",
  "assembly_accessions": ["isolate_01", "isolate_02", ...],
  "n_assemblies_at_creation": 42
}
```

`taxonomy_source` is `"ncbi"` when `--clade-name` resolved against the taxdump, or
`"user_supplied"` when it fell back to a name-only taxonomy or `--clade-taxonomy` was given
explicitly. `qc_applied` is `false` only when `--skip-qc` was used. Both fields default to
`"ncbi"`/`true` when absent, so pre-existing databases created before this mode existed
remain valid.

**Note:** this is provenance information only — it is surfaced in logs and metadata, but
does not gate routing. A `user_supplied` database with a fully-specified taxonomy still
routes normally; the only thing that affects routing is whether the taxonomy itself is
fully specified (see the routability note above).
