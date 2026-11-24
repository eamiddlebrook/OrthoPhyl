#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
viral_annotate.py

Per-genome annotator for eukaryotic viral assemblies, used by functions.sh's
VIRAL_PREDICT (the --annotator viral path) in place of Prodigal.

For one input genome FASTA (which may hold several sequences/segments):
  1. Try NCBI VADR (v-annotate.pl) against an ordered list of model libraries.
     The first library that assigns a model to at least one sequence wins.
     CDS feature nucleotide sequences are pulled from VADR's --out_allfasta
     output and translated in-house (VADR only emits nucleotide features).
     Fatal alerts do NOT trigger a VIGOR4 fallback -- they are only logged to
     the annotation_sources.tsv row (see CLAUDE.md's "log alerts" decision).
  2. If no library classifies the genome at all, fall back to VIGOR4. The
     reference database is picked automatically (diamond blastx against a
     pooled reference-protein DB, highest summed bitscore wins) unless
     --vigor-db names one explicitly. VIGOR4's own CDS/PEP translations are
     used as-is (not re-translated).
  3. Output is normalized to <genome_stem>_<n> record IDs, written in matching
     order to the CDS (--out-trans) and protein (--out-prots) FASTAs, mirroring
     Prodigal's own "<genome>_<n>" gene numbering so downstream DEDUP_annot_trans/
     FIX_TRANS_NAMES/FIX_PROTS_NAMES need no changes.

Parsing/translation/picking logic is pure (operates on strings/dicts already
read from disk) so it is unit-testable without running VADR/VIGOR4/diamond;
see tests/unit/test_viral_annotate.py and the real captured fixtures under
tests/fixtures/viral_annotate/.
"""
import argparse
import glob
import os
import shlex
import subprocess
import sys

STANDARD_CODON_TABLE = {}
_BASES = "TCAG"
_AAS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
for _i, _b1 in enumerate(_BASES):
    for _j, _b2 in enumerate(_BASES):
        for _k, _b3 in enumerate(_BASES):
            STANDARD_CODON_TABLE[_b1 + _b2 + _b3] = _AAS[_i * 16 + _j * 4 + _k]


def translate(nt_seq, table=1):
    """Translate a nucleotide CDS to protein using the standard genetic code.

    Tables 1 and 11 (the only ones relevant to viral/bacterial CDS translation
    here) share an identical sense/stop codon mapping -- they differ only in
    which codons may additionally serve as alternative starts, which does not
    affect in-frame translation of an already-delimited CDS. `table` is kept
    for interface clarity/future codes but is not currently used to vary the
    lookup.
    """
    nt_seq = nt_seq.upper().replace("U", "T")
    protein = []
    for pos in range(0, len(nt_seq) - len(nt_seq) % 3, 3):
        codon = nt_seq[pos:pos + 3]
        protein.append(STANDARD_CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def parse_fasta(text):
    """Minimal multi-FASTA parser. Returns an ordered list of (header, seq)
    where header is the full header line (no leading '>'), seq has no
    whitespace/newlines."""
    records = []
    header = None
    chunks = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


# ---------------------------------------------------------------------------
# VADR
# ---------------------------------------------------------------------------

def parse_sqc(text):
    """Parse a VADR .vadr.sqc table into a list of row dicts.

    Column layout is fixed by the installed VADR version (confirmed against
    a real 1.6.4 run): idx name len pf ant model1 grp1 subgrp1 score scnt cov
    subcov bias hits strand model2 grp2 subgrp2 diff diffnt alerts
    """
    cols = [
        "idx", "name", "len", "pf", "ant", "model1", "grp1", "subgrp1",
        "score", "scnt", "cov", "subcov", "bias", "hits", "strand",
        "model2", "grp2", "subgrp2", "diff", "diffnt", "alerts",
    ]
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < len(cols):
            continue
        row = dict(zip(cols, fields))
        rows.append(row)
    return rows


def parse_alt(text):
    """Parse a VADR .vadr.alt alerts table into a list of row dicts.

    Column layout (confirmed against a real 1.6.4 run's header, though the
    fixture run itself had zero alert rows): idx seqname model ftrtype
    ftrname ftridx code fail description seqcoords seqlen mdlcoords mdllen
    detail
    """
    cols = [
        "idx", "seqname", "model", "ftrtype", "ftrname", "ftridx", "code",
        "fail", "description", "seqcoords", "seqlen", "mdlcoords", "mdllen",
        "detail",
    ]
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(None, len(cols) - 1)
        if len(fields) < 9:
            # not every alert row necessarily fills every trailing column
            continue
        row = dict(zip(cols, fields))
        rows.append(row)
    return rows


def is_classified(sqc_row):
    """A sequence is "classified" by VADR if a model was assigned, regardless
    of pass/fail -- fatal alerts are logged, not treated as non-classification
    (see CLAUDE.md: "Use VADR CDS, log alerts")."""
    model = sqc_row.get("model1", "-")
    return bool(model) and model != "-"


def choose_vadr_library(library_results):
    """Pick the first VADR model library (in caller-supplied order) that
    classifies at least one sequence.

    library_results: ordered list of (mdir, mkey, sqc_rows).
    Returns the winning (mdir, mkey, sqc_rows) tuple, or None if no library
    classified anything.
    """
    for mdir, mkey, sqc_rows in library_results:
        if any(is_classified(r) for r in sqc_rows):
            return (mdir, mkey, sqc_rows)
    return None


def parse_vadr_feature_fasta(text):
    """Parse a VADR --out_allfasta feature FASTA (e.g. *.vadr.<model>.CDS.*.fa).

    Real header format (confirmed): >gi|770146038|gb|KM031079.1|/CDS.1/1..515:+
    The sequence ID is everything before the first '/'; it may itself contain
    '|' characters (NCBI gi-style headers), which is fine since we discard it
    after use (replaced by a <genome_stem>_<n> name in normalize_records).

    Returns an ordered list of (seq_id, feature_label, nt_seq).
    """
    out = []
    for header, seq in parse_fasta(text):
        seq_id, _, rest = header.partition("/")
        feature_label = rest.split("/", 1)[0] if rest else ""
        out.append((seq_id, feature_label, seq.replace("-", "")))
    return out


def collect_vadr_cds(vadr_out_dir, out_prefix):
    """Glob every CDS feature FASTA VADR wrote for one run and parse them.

    VADR writes one file per model actually used, e.g.
    <prefix>.vadr.NC_001474.CDS.1.fa and <prefix>.vadr.NC_001477.CDS.1.fa when
    two sequences in the input were assigned different models within the same
    library. Sorted by (seq_id, feature_label) for a deterministic gene order
    independent of filesystem glob ordering.
    """
    pattern = os.path.join(vadr_out_dir, f"{out_prefix}.vadr.*.CDS.*.fa")
    records = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as fh:
            records.extend(parse_vadr_feature_fasta(fh.read()))
    records.sort(key=lambda r: (r[0], r[1]))
    return records


# ---------------------------------------------------------------------------
# VIGOR4
# ---------------------------------------------------------------------------

def parse_vigor_fasta(text):
    """Parse a VIGOR4 .cds/.pep FASTA into {record_id: seq}, keyed by the
    first whitespace-delimited token of the header (the accession, or
    <accession>.<n> for mat_peptide records)."""
    out = {}
    for header, seq in parse_fasta(text):
        rec_id = header.split(None, 1)[0]
        out[rec_id] = seq
    return out


def filter_vigor_pep(cds_records, pep_records):
    """Keep only .pep entries whose ID exactly matches a .cds ID.

    Real VIGOR4 output (confirmed): .pep contains the full top-level CDS
    translation under the SAME id as .cds (e.g. "NC_001477.1"), plus extra
    "mat_peptide"-tagged records under suffixed ids (e.g. "NC_001477.1.1",
    ".2", ...) for individual cleaved peptides of a polyprotein. We only want
    the full CDS-level translations here, matched 1:1 against .cds.
    """
    return {rid: pep_records[rid] for rid in cds_records if rid in pep_records}


def pick_vigor_db(db_scores):
    """Given {db_name: summed_bitscore}, pick the best-scoring DB.

    Ties broken by ascending db_name for determinism. Returns None if
    db_scores is empty or every score is <= 0 (no usable hits)."""
    usable = {db: score for db, score in db_scores.items() if score and score > 0}
    if not usable:
        return None
    best_score = max(usable.values())
    best = sorted(db for db, score in usable.items() if score == best_score)
    return best[0]


def parse_diamond_outfmt6(text):
    """Aggregate a diamond blastx outfmt6 table into {ref_db_name: summed_bitscore}.

    Expects the pooled reference DB's sequence ids to be tagged
    "<db_name>|<original_id>" (see build_vigor_diamond_db), so the ref db a
    hit came from can be recovered from column 2 (sseqid) alone.
    """
    scores = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        sseqid = fields[1]
        bitscore = float(fields[11])
        db_name = sseqid.split("|", 1)[0]
        scores[db_name] = scores.get(db_name, 0.0) + bitscore
    return scores


# ---------------------------------------------------------------------------
# Shared normalization
# ---------------------------------------------------------------------------

def normalize_records(genome_stem, nt_records, pep_records):
    """Rename matched CDS/protein pairs to <genome_stem>_<n> (1-indexed, in
    nt_records' iteration order), dropping any pair whose protein contains an
    internal stop codon (premature stop before the final residue) -- a sign
    of a frameshift/misassembly that would otherwise corrupt downstream SCO
    membership.

    nt_records / pep_records: ordered lists of (id, seq) with matching ids in
    matching order (callers are responsible for pairing/ordering upstream).

    Returns (trans_fasta_records, prots_fasta_records, n_dropped), each a list
    of (name, seq) ready to write out, and n_dropped the count of CDS dropped
    for an internal stop.
    """
    trans_out = []
    prots_out = []
    n_dropped = 0
    n = 0
    for (nt_id, nt_seq), (pep_id, pep_seq) in zip(nt_records, pep_records):
        pep_seq = pep_seq.rstrip("*")
        if "*" in pep_seq:
            n_dropped += 1
            continue
        n += 1
        name = f"{genome_stem}_{n}"
        trans_out.append((name, nt_seq))
        prots_out.append((name, pep_seq))
    return trans_out, prots_out, n_dropped


def write_fasta(records, path):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


# ---------------------------------------------------------------------------
# Subprocess-isolated I/O (not unit tested directly -- thin wrappers only)
# ---------------------------------------------------------------------------

def _run(cmd, env_run_prefix=""):
    full_cmd = (shlex.split(env_run_prefix) if env_run_prefix else []) + cmd
    return subprocess.run(full_cmd, capture_output=True, text=True)


def run_vadr(genome_fa, mdir, mkey, workdir, out_prefix, vadr_options="", env_run_prefix=""):
    """Run v-annotate.pl for one model library. Returns the output directory
    (removed/recreated first since v-annotate.pl refuses to run into an
    existing directory).

    v-annotate.pl names every output file after the basename of its output
    directory argument (e.g. "<out_dir_basename>.vadr.sqc") -- out_dir must
    therefore be workdir/out_prefix exactly, not some other name, or callers
    looking for "<out_prefix>.vadr.sqc" won't find it."""
    out_dir = os.path.join(workdir, out_prefix)
    if os.path.isdir(out_dir):
        import shutil
        shutil.rmtree(out_dir)
    cmd = ["v-annotate.pl", "--mdir", mdir, "--mkey", mkey, "--out_allfasta"]
    if vadr_options:
        cmd += shlex.split(vadr_options)
    cmd += [genome_fa, out_dir]
    result = _run(cmd, env_run_prefix)
    return out_dir, result


def build_vigor_diamond_db(vigor_db_dir, out_dmnd, env_run_prefix=""):
    """Pool every VIGOR_DB Reference_DBs/*_db protein reference FASTA into one
    diamond DB, tagging each sequence "<db_name>|<original_id>" so a blastx
    hit can be traced back to its source DB (see parse_diamond_outfmt6).

    A handful of real VIGOR_DB entries have their header and sequence
    squashed onto one physical line with no newline between them (confirmed
    on dengue_txid12637_db's NP_073286.1 and a few others) -- parse_fasta
    reads these as a header with a zero-length sequence. `diamond makedb`
    hard-errors ("sequence of length 0") on even one such record, which would
    otherwise take down pooled-DB auto-detection for every genome over a
    handful of bad upstream records. Skip them (logged) rather than guess
    where metadata ends and sequence begins in the malformed line."""
    n_skipped = 0
    pooled_fa = out_dmnd + ".pooled.fa"
    with open(pooled_fa, "w") as out_fh:
        for ini_path in sorted(glob.glob(os.path.join(vigor_db_dir, "*_db.ini"))):
            db_name = os.path.basename(ini_path)[:-len("_db.ini")] + "_db"
            fa_path = os.path.join(vigor_db_dir, db_name)
            if not os.path.isfile(fa_path):
                continue
            with open(fa_path) as fh:
                for header, seq in parse_fasta(fh.read()):
                    if not seq:
                        n_skipped += 1
                        continue
                    rec_id = header.split(None, 1)[0]
                    out_fh.write(f">{db_name}|{rec_id}\n{seq}\n")
    if n_skipped:
        sys.stderr.write(f"build_vigor_diamond_db: skipped {n_skipped} zero-length record(s) "
                          f"(malformed header/sequence line in a VIGOR_DB reference file)\n")
    cmd = ["diamond", "makedb", "--in", pooled_fa, "-d", out_dmnd]
    return _run(cmd, env_run_prefix)


def run_diamond_blastx(genome_fa, dmnd_db, out_path, env_run_prefix=""):
    cmd = ["diamond", "blastx", "-q", genome_fa, "-d", dmnd_db, "-o", out_path,
           "--outfmt", "6", "--quiet"]
    return _run(cmd, env_run_prefix)


def run_vigor4(genome_fa, db, out_prefix, config_file, env_run_prefix=""):
    cmd = ["vigor4", "-i", genome_fa, "-o", out_prefix, "-d", db,
           "--config-file", config_file, "--overwrite-output"]
    return _run(cmd, env_run_prefix)


def resolve_exonerate_path(env_run_prefix=""):
    """VIGOR4 needs exonerate_path to be an absolute filesystem path -- it does
    NOT fall back to a bare command name via its own PATH lookup (confirmed:
    "exonerate_path exonerate does not exist"). Resolve it with `which` run
    through the same env_run_prefix VIGOR4 itself will use, so this finds the
    exonerate binary inside the pinned viral-tools env rather than whatever is
    on the calling (orthophyl) process's PATH."""
    result = _run(["which", "exonerate"], env_run_prefix)
    path = result.stdout.strip()
    return path if result.returncode == 0 and path else None


def build_vigor_config(vigor_db_dir, workdir, genome_stem, user_config_file="", env_run_prefix=""):
    """VIGOR4 crashes outright on an empty/missing --config-file (Java
    NullPointerException: "virus specific config path not specified"), and
    even a real config needs reference_database_path (must match
    --vigor-db-dir) and an absolute exonerate_path set explicitly. Rather than
    requiring the caller to hand-author and keep a config file in sync with
    --vigor-db-dir/--env-run, generate the effective config here every run:
    reference_database_path/temporary_directory are always ours (so they can
    never drift from the current --vigor-db-dir/--workdir); any other lines
    from a user-supplied --vigor-config (e.g. a non-default exonerate_path)
    are preserved as-is."""
    lines = []
    have_exonerate_path = False
    if user_config_file and os.path.isfile(user_config_file):
        with open(user_config_file) as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                key = stripped.split("=", 1)[0].strip()
                if key in ("reference_database_path", "temporary_directory"):
                    continue
                if key == "exonerate_path":
                    have_exonerate_path = True
                lines.append(stripped)
    if not have_exonerate_path:
        exonerate_path = resolve_exonerate_path(env_run_prefix)
        if exonerate_path:
            lines.append(f"exonerate_path={exonerate_path}")
    lines.append(f"reference_database_path={vigor_db_dir}")
    lines.append(f"temporary_directory={workdir}")
    out_path = os.path.join(workdir, f"{genome_stem}.vigor.ini")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return out_path


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def annotate_via_vadr(genome_fa, genome_stem, vadr_libraries, workdir, vadr_options="", env_run_prefix=""):
    """Try each VADR library in order. Returns a dict with keys:
    classified (bool), mdir, mkey, sqc_rows, alt_rows, cds_records (list of
    (seq_id, feature_label, nt_seq)), or classified=False if none matched."""
    library_results = []
    all_outputs = {}
    for lib in vadr_libraries:
        mdir, _, mkey = lib.partition(":")
        mkey = mkey or "vadr"
        out_prefix = f"{genome_stem}.{mkey}"
        out_dir, result = run_vadr(genome_fa, mdir, mkey, workdir, out_prefix,
                                    vadr_options, env_run_prefix)
        sqc_path = os.path.join(out_dir, f"{out_prefix}.vadr.sqc")
        sqc_rows = []
        if os.path.isfile(sqc_path):
            with open(sqc_path) as fh:
                sqc_rows = parse_sqc(fh.read())
        library_results.append((mdir, mkey, sqc_rows))
        all_outputs[mkey] = (out_dir, out_prefix, result)

    winner = choose_vadr_library(library_results)
    if winner is None:
        return {"classified": False}

    mdir, mkey, sqc_rows = winner
    out_dir, out_prefix, result = all_outputs[mkey]
    alt_path = os.path.join(out_dir, f"{out_prefix}.vadr.alt")
    alt_rows = []
    if os.path.isfile(alt_path):
        with open(alt_path) as fh:
            alt_rows = parse_alt(fh.read())
    cds_records = collect_vadr_cds(out_dir, out_prefix)
    return {
        "classified": True,
        "mdir": mdir,
        "mkey": mkey,
        "sqc_rows": sqc_rows,
        "alt_rows": alt_rows,
        "cds_records": cds_records,
    }


def annotate_via_vigor4(genome_fa, genome_stem, vigor_db, vigor_db_dir, config_file,
                         workdir, dmnd_db=None, env_run_prefix=""):
    """Auto-detect (if vigor_db=="auto") or use the given VIGOR_DB reference
    database, run VIGOR4, and return a dict with keys: db, cds_records,
    pep_records (both {id: seq}, already filtered to CDS-level ids), or
    failed=True/reason if nothing usable came back."""
    chosen_db = vigor_db
    if vigor_db == "auto":
        if dmnd_db is None:
            return {"failed": True, "reason": "no diamond DB built for auto-detection"}
        hits_path = os.path.join(workdir, f"{genome_stem}.vigor_refs.outfmt6")
        run_diamond_blastx(genome_fa, dmnd_db, hits_path, env_run_prefix)
        db_scores = {}
        if os.path.isfile(hits_path):
            with open(hits_path) as fh:
                db_scores = parse_diamond_outfmt6(fh.read())
        chosen_db = pick_vigor_db(db_scores)
        if chosen_db is None:
            return {"failed": True, "reason": "no diamond hits against any VIGOR_DB reference"}

    out_prefix = os.path.join(workdir, f"{genome_stem}.vigor4")
    effective_config = build_vigor_config(vigor_db_dir, workdir, genome_stem, config_file, env_run_prefix)
    run_vigor4(genome_fa, chosen_db, out_prefix, effective_config, env_run_prefix)
    cds_path = out_prefix + ".cds"
    pep_path = out_prefix + ".pep"
    if not os.path.isfile(cds_path) or not os.path.isfile(pep_path):
        return {"failed": True, "reason": f"vigor4 produced no .cds/.pep for db {chosen_db}"}

    with open(cds_path) as fh:
        cds_records = parse_vigor_fasta(fh.read())
    with open(pep_path) as fh:
        pep_records = parse_vigor_fasta(fh.read())
    pep_records = filter_vigor_pep(cds_records, pep_records)
    if not pep_records:
        return {"failed": True, "reason": f"vigor4 db {chosen_db} produced zero matched CDS/PEP pairs"}
    return {"failed": False, "db": chosen_db, "cds_records": cds_records, "pep_records": pep_records}


def annotate_genome(args):
    genome_fa = args.genome
    genome_stem = args.genome_stem
    workdir = args.workdir
    os.makedirs(workdir, exist_ok=True)

    vadr_libraries = shlex.split(args.vadr_models) if args.vadr_models else []
    vadr_result = {"classified": False}
    if vadr_libraries:
        vadr_result = annotate_via_vadr(
            genome_fa, genome_stem, vadr_libraries, workdir,
            vadr_options=args.vadr_options, env_run_prefix=args.env_run,
        )

    method = "failed"
    library_or_db = "-"
    model = "-"
    n_cds = 0
    vadr_status = "NA"
    fatal_alerts = "-"

    if vadr_result.get("classified"):
        method = "vadr"
        library_or_db = vadr_result["mkey"]
        models_hit = sorted({r["model1"] for r in vadr_result["sqc_rows"] if is_classified(r)})
        model = ",".join(models_hit)
        vadr_status = "fail" if any(r.get("pf") == "FAIL" for r in vadr_result["sqc_rows"]) else "pass"
        alerts = sorted({r["code"] for r in vadr_result.get("alt_rows", [])})
        fatal_alerts = ",".join(alerts) if alerts else "-"

        nt_records = [(f"{seq_id}/{feat}", seq) for seq_id, feat, seq in vadr_result["cds_records"]]
        pep_records = [(name, translate(seq)) for name, seq in nt_records]
        trans_out, prots_out, n_dropped = normalize_records(genome_stem, nt_records, pep_records)
        n_cds = len(trans_out)
        if n_dropped:
            sys.stderr.write(f"{genome_stem}: dropped {n_dropped} VADR CDS with internal stop codons\n")
    else:
        vigor_result = annotate_via_vigor4(
            genome_fa, genome_stem, args.vigor_db, args.vigor_db_dir, args.vigor_config,
            workdir, dmnd_db=args.vigor_refs_dmnd, env_run_prefix=args.env_run,
        )
        if not vigor_result.get("failed"):
            method = "vigor4"
            library_or_db = vigor_result["db"]
            ids = sorted(vigor_result["pep_records"])
            model = ",".join(ids)
            nt_records = [(i, vigor_result["cds_records"][i]) for i in ids]
            pep_records = [(i, vigor_result["pep_records"][i]) for i in ids]
            trans_out, prots_out, n_dropped = normalize_records(genome_stem, nt_records, pep_records)
            n_cds = len(trans_out)
            if n_dropped:
                sys.stderr.write(f"{genome_stem}: dropped {n_dropped} VIGOR4 CDS with internal stop codons\n")
        else:
            trans_out, prots_out = [], []
            sys.stderr.write(f"{genome_stem}: VIGOR4 fallback failed ({vigor_result.get('reason')})\n")

    if trans_out:
        write_fasta(trans_out, args.out_trans)
        write_fasta(prots_out, args.out_prots)

    row = "\t".join([genome_stem, method, library_or_db, model, str(n_cds), vadr_status, fatal_alerts])
    if args.annot_row:
        with open(args.annot_row, "w") as fh:
            fh.write(row + "\n")
    print(row)
    return 0 if trans_out else 1


def build_arg_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome", help="path to one genome FASTA (may hold multiple segments)")
    p.add_argument("--genome-stem", help="basename used for output record naming, e.g. <genome_stem>_1")
    p.add_argument("--out-trans", help="output CDS nucleotide FASTA path")
    p.add_argument("--out-prots", help="output protein FASTA path")
    p.add_argument("--annot-row", help="path to write this genome's annotation_sources.tsv row")
    p.add_argument("--workdir", help="scratch directory for VADR/VIGOR4 intermediate output")
    p.add_argument("--vadr-models", default="", help="space-separated ordered 'DIR[:MKEY] ...' list")
    p.add_argument("--vadr-options", default="", help="extra v-annotate.pl options")
    p.add_argument("--vigor-db", default="auto", help="VIGOR_DB reference database name, or 'auto'")
    p.add_argument("--vigor-db-dir", default="", help="path to VIGOR_DB's Reference_DBs")
    p.add_argument("--vigor-config", default="", help="VIGOR4 --config-file path")
    p.add_argument("--vigor-refs-dmnd", default=None, help="prebuilt pooled diamond DB (see build_vigor_diamond_db), required for --vigor-db auto")
    p.add_argument("--env-run", default="", help="command prefix to run viral tools (e.g. 'micromamba run -n viral_annot'), empty = use PATH")
    p.add_argument("--build-vigor-refs-db", action="store_true",
                    help="instead of annotating a genome, pool --vigor-db-dir's reference proteins into one "
                         "diamond DB at --out (for later --vigor-db auto calls) and exit")
    p.add_argument("--out", help="output path for --build-vigor-refs-db")
    return p


def main(argv=None):
    p = build_arg_parser()
    args = p.parse_args(argv)

    if args.build_vigor_refs_db:
        if not args.vigor_db_dir or not args.out:
            p.error("--build-vigor-refs-db requires --vigor-db-dir and --out")
        result = build_vigor_diamond_db(args.vigor_db_dir, args.out, env_run_prefix=args.env_run)
        if result.returncode != 0:
            sys.stderr.write(result.stderr)
            return 1
        return 0

    missing = [name for name, val in (
        ("--genome", args.genome), ("--genome-stem", args.genome_stem),
        ("--out-trans", args.out_trans), ("--out-prots", args.out_prots),
        ("--workdir", args.workdir),
    ) if not val]
    if missing:
        p.error(f"missing required argument(s): {', '.join(missing)}")

    return annotate_genome(args)


if __name__ == "__main__":
    sys.exit(main())
