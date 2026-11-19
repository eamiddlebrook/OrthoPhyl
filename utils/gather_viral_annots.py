#!/usr/bin/env python3
"""Fetch NCBI's existing CDS+protein annotations for a viral taxon.

Unlike `gather_filter_asms.sh` (built around GC[AF] assembly accessions,
rehydration, and CheckM2), viral taxa are fetched through NCBI's separate
"virus genome" `datasets` endpoint, which returns pooled `cds.fna` /
`protein.faa` / `annotation_report.jsonl` / `data_report.jsonl` files -- not
per-genome assemblies. This script splits that pooled package back into
per-genome `<id>.fna` (CDS nucleotide) / `<id>.faa` (protein) pairs suitable
for OrthoPhyl.sh's/ReLeaf.sh's `-a trans_dir,prots_dir` input, so the
expensive VADR/VIGOR4 annotation path (see `python_scripts/viral_annotate.py`)
only ever runs on user-supplied query genomes, never on all of NCBI.

Key format facts this module's logic depends on, confirmed against real
`datasets download virus genome taxon ... --include cds,protein,annotation`
output (not guessed -- see the session notes for the exact commands run):

- `cds.fna` / `protein.faa` are keyed by a `seqId` token
  (e.g. ``PV275052.1:5-5104``) that appears as the first whitespace-delimited
  word after '>'. These files pool EVERY genome in the taxon together.
- `annotation_report.jsonl` has one line per nucleotide accession, with
  ``genes[].cds[]`` entries carrying the ``nucleotide.seqId``/``protein.seqId``
  pair that maps a CDS's feature range back into `cds.fna`/`protein.faa`.
  A CDS whose product is itself a polyprotein additionally lists its
  cleavage products under a nested ``maturePeptide`` key -- those are
  sub-ranges of the already-captured full CDS/protein pair, not separate
  genes, and are deliberately NOT extracted here.
- `data_report.jsonl` has one line per nucleotide accession too, with
  ``isolate.name`` and (for segmented viruses only) a top-level ``segment``
  field. Confirmed against a real segmented genome (Lassa mammarenavirus,
  segments "L"/"S"): both segment accessions shared one `isolate.name` and
  had NO common BioSample field (virus genome packages carry no BioSample at
  all) -- so `isolate.name` is the only cross-accession grouping key NCBI
  gives us for multi-segment genomes. Falls back to the record's own
  accession (always unique) when no isolate name is present.
- **`cds.fna` gives the literal genomic-strand span for a CDS feature and
  does NOT reverse-complement minus-strand ORFs.** Confirmed against a real
  ambisense arenavirus (Lassa virus S segment): GPC translates forward, but
  NP -- on the same accession, same coordinate convention (``begin<end``
  either way) -- only translates correctly on the reverse complement.
  Neither `annotation_report.jsonl` nor `data_report.jsonl` carries an
  explicit strand field for a CDS, so orientation can't be read off
  metadata; `resolve_cds_strand` instead translates both orientations and
  keeps whichever one actually matches the paired `protein.faa` sequence
  (which NCBI _does_ give correctly oriented).
"""

import argparse
import json
import re
import statistics
import subprocess
import sys
import zipfile
from pathlib import Path

STOP_CODONS = {"TAA", "TAG", "TGA"}
_UNSAFE_CHARS = re.compile(r"[^A-Za-z0-9_.-]")

# Standard genetic code (NCBI translation table 1) -- sufficient for the
# orientation-matching heuristic in resolve_cds_strand, which only needs an
# approximate match score against the real (correctly-oriented) protein.
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
for i, b1 in enumerate(_BASES):
    for j, b2 in enumerate(_BASES):
        for k, b3 in enumerate(_BASES):
            CODON_TABLE[b1 + b2 + b3] = _AAS[i * 16 + j * 4 + k]


# --------------------------------------------------------------------------- #
# Pure parsing / grouping / QC (unit-tested directly, no subprocess/network). #
# --------------------------------------------------------------------------- #

def parse_jsonlines(text):
    """Parse newline-delimited JSON, skipping blank lines and unparseable
    noise (e.g. a stray `datasets` CLI update notice sharing the stream)."""
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return records


def load_fasta(path):
    """Parse a multi-FASTA into ``{seqid: (header, sequence)}``.

    ``seqid`` is the token up to the first whitespace after '>' -- matches
    how `datasets`' cds.fna/protein.faa key each record
    (e.g. ``"PV275052.1:5-5104"``).
    """
    records = {}
    seqid = None
    header = None
    chunks = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if seqid is not None:
                    records[seqid] = (header, "".join(chunks))
                header = line[1:]
                seqid = header.split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if seqid is not None:
            records[seqid] = (header, "".join(chunks))
    return records


def genome_group_key(data_report_rec):
    """Return the cross-accession grouping key for one data_report record.

    See the module docstring: isolate name when present, else the record's
    own (always-unique) accession.
    """
    isolate = (data_report_rec.get("isolate") or {}).get("name", "") or ""
    isolate = isolate.strip()
    if isolate:
        return isolate
    return data_report_rec["accession"]


def sanitize_genome_id(raw):
    """Turn a free-text isolate name/accession into a safe filename stem."""
    safe = _UNSAFE_CHARS.sub("_", raw.strip())
    safe = safe.strip("_.")
    return safe or "unnamed"


def dedupe_prefer_refseq(data_report_recs):
    """Collapse exact-duplicate sequences, keeping the RefSeq copy.

    Virus genome packages have no assembly-accession pairing like bacterial
    GCF_/GCA_, so duplicate detection instead uses NCBI's own
    ``nucleotide.sequenceHash`` (identical sequence => identical hash). When
    both a RefSeq and a GenBank record share a hash, the RefSeq one is kept;
    otherwise the first one seen wins. Preserves input order.
    """
    seen_hash_idx = {}
    kept = []
    for rec in data_report_recs:
        h = (rec.get("nucleotide") or {}).get("sequenceHash")
        if not h:
            kept.append(rec)
            continue
        if h not in seen_hash_idx:
            seen_hash_idx[h] = len(kept)
            kept.append(rec)
        else:
            idx = seen_hash_idx[h]
            if (rec.get("sourceDatabase") == "RefSeq"
                    and kept[idx].get("sourceDatabase") != "RefSeq"):
                kept[idx] = rec
    return kept


def extract_cds_entries(annotation_rec):
    """Yield one dict per top-level CDS in an annotation_report record.

    ``maturePeptide`` sub-entries are deliberately excluded -- see the module
    docstring.
    """
    accession = annotation_rec.get("accession")
    for gene in annotation_rec.get("genes", []):
        gene_name = gene.get("name", "")
        for cds in gene.get("cds", []):
            nt = cds.get("nucleotide") or {}
            prot = cds.get("protein") or {}
            nt_seqid = nt.get("seqId")
            prot_seqid = prot.get("seqId")
            if not nt_seqid or not prot_seqid:
                continue
            yield {
                "accession": accession,
                "gene_name": cds.get("name") or gene_name,
                "nt_seqid": nt_seqid,
                "prot_seqid": prot_seqid,
            }


def _unique_sanitized_id(raw_key, used_ids):
    base = sanitize_genome_id(raw_key)
    candidate = base
    i = 2
    while candidate in used_ids:
        candidate = f"{base}_{i}"
        i += 1
    used_ids.add(candidate)
    return candidate


def build_genomes(data_report_recs, annotation_recs):
    """Group per-accession CDS/protein pairs into per-genome gene lists.

    Returns ``{genome_id: {"accessions": [...], "representative_accession":
    acc, "entries": [...]}}``, keyed by a unique, filesystem-safe genome ID.
    Genomes with zero resolvable CDS entries (e.g. an accession present in
    data_report.jsonl but missing from annotation_report.jsonl) are omitted.
    """
    ann_by_accession = {r.get("accession"): r for r in annotation_recs}

    groups = {}
    for rec in data_report_recs:
        accession = rec.get("accession")
        if not accession:
            continue
        raw_key = genome_group_key(rec)
        groups.setdefault(raw_key, []).append(accession)

    genomes = {}
    used_ids = set()
    for raw_key, accessions in groups.items():
        entries = []
        for accession in accessions:
            ann = ann_by_accession.get(accession)
            if ann is None:
                continue
            entries.extend(extract_cds_entries(ann))
        if not entries:
            continue
        genome_id = _unique_sanitized_id(raw_key, used_ids)
        genomes[genome_id] = {
            "accessions": accessions,
            "representative_accession": accessions[0],
            "entries": entries,
        }
    return genomes


def filter_genomes_by_gene_count(genomes, min_frac=0.5):
    """Drop genomes with unusually few CDS relative to the taxon's norm.

    A crude per-genome QC stand-in for CheckM2 (which only works on cellular
    genomes, see ``viral_min_genomes``/the wrapper's auto-``skip_qc`` for
    viral `--genome-dir` builds): drop any genome whose CDS count is below
    ``min_frac`` times the across-genome median. Returns ``(kept, dropped_ids)``.
    """
    if not genomes:
        return {}, []
    counts = {gid: len(g["entries"]) for gid, g in genomes.items()}
    median = statistics.median(counts.values())
    threshold = median * min_frac
    kept = {gid: g for gid, g in genomes.items() if counts[gid] >= threshold}
    dropped = [gid for gid in genomes if gid not in kept]
    return kept, dropped


def has_internal_stop(nt_seq):
    """True if an in-frame stop codon appears before the final codon."""
    seq = nt_seq.upper().replace("\n", "")
    usable_len = len(seq) - (len(seq) % 3)
    codons = [seq[i:i + 3] for i in range(0, usable_len, 3)]
    if not codons:
        return False
    return any(c in STOP_CODONS for c in codons[:-1])


def ambiguous_base_fraction(nt_seq):
    """Fraction of non-ACGT characters in a nucleotide sequence."""
    seq = nt_seq.upper()
    if not seq:
        return 0.0
    n_ambig = sum(1 for c in seq if c not in "ACGT")
    return n_ambig / len(seq)


def reverse_complement(nt_seq):
    """Reverse-complement a nucleotide sequence (case-preserving, N-safe)."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return nt_seq.translate(comp)[::-1]


def translate(nt_seq, table=CODON_TABLE):
    """Translate a nucleotide sequence in frame 0. Unknown codons -> 'X'."""
    seq = nt_seq.upper().replace("\n", "")
    usable_len = len(seq) - (len(seq) % 3)
    return "".join(
        table.get(seq[i:i + 3], "X") for i in range(0, usable_len, 3)
    )


def _match_score(translated, aa_seq):
    """Fraction of identical residues over the longer of the two sequences
    (trailing stops stripped from both). Crude but enough to tell a correct
    orientation from a scrambled one -- see resolve_cds_strand."""
    t = translated.rstrip("*")
    a = aa_seq.rstrip("*")
    if not t or not a:
        return 0.0
    n = min(len(t), len(a))
    matches = sum(1 for x, y in zip(t[:n], a[:n]) if x == y)
    return matches / max(len(t), len(a))


def resolve_cds_strand(nt_seq, aa_seq, min_score=0.5):
    """Pick whichever orientation/frame of ``nt_seq`` translates to match
    ``aa_seq``, trimming leading bases as needed.

    Two independent, real-data-confirmed quirks make this necessary:
    - NCBI's `cds.fna` does not reverse-complement minus-strand/ambisense CDS
      features (see module docstring), so the forward, frame-0 extraction is
      only sometimes the correct orientation.
    - Some CDS ranges (confirmed on a Norovirus ORF1/ORF2 ribosomal-frameshift
      junction, e.g. `nucleotide.range` ``1-5099`` translating cleanly only
      2 bases into the record) start a base or two before the real
      translation start -- a side effect of how NCBI records the CDS span
      across a programmed ribosomal frameshift -- so frame 0 on the correct
      strand still is not guaranteed to be the correct frame.

    Tries all 3 frames on both the given sequence and its reverse
    complement (6 candidates), scores each against ``aa_seq``, and returns
    whichever trimmed sequence scores highest -- or ``None`` if none of the
    6 match ``aa_seq`` well enough to trust.
    """
    rev_seq = reverse_complement(nt_seq)
    best_seq, best_score = None, -1.0
    for seq in (nt_seq, rev_seq):
        for frame in range(3):
            trimmed = seq[frame:]
            score = _match_score(translate(trimmed), aa_seq)
            if score > best_score:
                best_seq, best_score = trimmed, score
    if best_score < min_score:
        return None
    return best_seq


def write_genome_fastas(genomes, cds_seqs, prot_seqs, nucl_dir, prot_dir,
                         max_ambiguous_frac=0.05):
    """Write one ``<genome_id>.fna``/``.faa`` pair per genome.

    Genes are renumbered ``<genome_id>_<n>`` (1-based, in entry order) so
    nucleotide and protein headers match exactly, mirroring the convention
    the rest of OrthoPhyl expects from an annotator
    (`python_scripts/viral_annotate.py`'s ``normalize_records`` does the same
    for the VADR/VIGOR4 path). A gene is skipped (and logged) when its
    sequence is missing from the pooled FASTAs, contains an in-frame stop
    before its end, or exceeds ``max_ambiguous_frac`` non-ACGT bases. A
    genome that loses every gene this way is dropped entirely.

    Returns ``{genome_id: n_genes_written}``.
    """
    nucl_dir = Path(nucl_dir)
    prot_dir = Path(prot_dir)
    nucl_dir.mkdir(parents=True, exist_ok=True)
    prot_dir.mkdir(parents=True, exist_ok=True)

    written = {}
    for genome_id, g in genomes.items():
        nt_lines = []
        aa_lines = []
        n = 0
        for entry in g["entries"]:
            nt_rec = cds_seqs.get(entry["nt_seqid"])
            aa_rec = prot_seqs.get(entry["prot_seqid"])
            if nt_rec is None or aa_rec is None:
                print(f"[gather_viral_annots] {genome_id}: missing sequence "
                      f"for {entry['nt_seqid']} / {entry['prot_seqid']}, "
                      f"skipping gene", file=sys.stderr)
                continue
            _, nt_seq = nt_rec
            _, aa_seq = aa_rec
            resolved = resolve_cds_strand(nt_seq, aa_seq)
            if resolved is None:
                print(f"[gather_viral_annots] {genome_id}: could not resolve "
                      f"strand orientation for {entry['nt_seqid']} (neither "
                      f"orientation matches {entry['prot_seqid']}), "
                      f"skipping gene", file=sys.stderr)
                continue
            nt_seq = resolved
            if has_internal_stop(nt_seq):
                print(f"[gather_viral_annots] {genome_id}: internal stop in "
                      f"{entry['nt_seqid']}, skipping gene", file=sys.stderr)
                continue
            if ambiguous_base_fraction(nt_seq) > max_ambiguous_frac:
                print(f"[gather_viral_annots] {genome_id}: too many "
                      f"ambiguous bases in {entry['nt_seqid']}, skipping "
                      f"gene", file=sys.stderr)
                continue
            n += 1
            gene_id = f"{genome_id}_{n}"
            nt_lines.append(f">{gene_id}\n{nt_seq}\n")
            aa_lines.append(f">{gene_id}\n{aa_seq.rstrip('*')}\n")
        if n == 0:
            print(f"[gather_viral_annots] {genome_id}: no usable genes, "
                  f"dropping", file=sys.stderr)
            continue
        (nucl_dir / f"{genome_id}.fna").write_text("".join(nt_lines))
        (prot_dir / f"{genome_id}.faa").write_text("".join(aa_lines))
        written[genome_id] = n
    return written


def write_metadata_tsv(written_genome_ids, genomes, data_report_recs, out_path):
    """Write ``all_asm_acc_metadata`` in the existing space-delimited format.

    Column order (see ``parse_asm_metadata.py``'s ``project_metadata``, which
    every downstream reader -- `annotate_tree_by_country.py`, the wrapper's
    ``--geo`` metadata copy -- keys off of by position, not header)::

        accession species species.strain biosample taxid taxcheck
        excl_from_refseq geo_loc collected_by sequenced_by biosample_owner
        asm_submitter lat_lon

    Virus genome packages carry no BioSample accession of their own, so
    ``accession`` and ``biosample`` are both set to the sanitized genome ID
    -- the same stem used for ``<genome_id>.fna``/``.faa`` -- so a later
    ``--geo`` lookup (keyed on column 1) resolves against the same ID used
    everywhere else in this pipeline. ``taxcheck``/``excl_from_refseq`` and
    the trailing submitter-provenance columns (9-13, see
    ``geo_provenance.py``) have no virus-package equivalent and are set to
    "NA" -- NOT "", which a whitespace-`.split()` reader (e.g.
    ``annotate_tree_by_country.py``'s ``parse_metadata``) collapses away,
    silently shifting every later column left.

    Every field is whitespace-collapsed-then-underscored, matching
    ``parse_asm_metadata.py``'s ``sanitize_field`` -- a species name
    ("Norovirus GII") or geo_loc value ("Japan: Tokyo") containing a raw
    space would otherwise split into extra columns and shift every field
    after it for whitespace-splitting readers like this file's own
    ``annotate_tree_by_country.py``.
    """
    by_accession = {r.get("accession"): r for r in data_report_recs}
    lines = []
    for genome_id in written_genome_ids:
        g = genomes[genome_id]
        rec = by_accession.get(g["representative_accession"], {})
        virus = rec.get("virus") or {}
        species = _sanitize_field(virus.get("organismName", ""))
        taxid = str(virus.get("taxId", ""))
        geo = _sanitize_field((rec.get("location") or {}).get("geographicLocation", ""))
        lines.append(" ".join([
            genome_id, species, f"{species}.{genome_id}", genome_id,
            taxid, "NA", "NA", geo, "NA", "NA", "NA", "NA", "NA",
        ]))
    Path(out_path).write_text("".join(f"{l}\n" for l in lines))


def _sanitize_field(value):
    """Collapse internal whitespace to underscores so a field can never
    split a space-delimited row (mirrors parse_asm_metadata.py's
    sanitize_field)."""
    return "_".join(str(value).split())


# --------------------------------------------------------------------------- #
# Thin subprocess/filesystem wrappers (not unit-tested against real network). #
# --------------------------------------------------------------------------- #

def fetch_virus_summary(taxon, datasets_bin="datasets", extra_filters=None):
    """Run ``datasets summary virus genome taxon`` and return parsed records.

    Used only to count matching genomes before committing to a full
    download (the `--max-total-genomes` guardrail).
    """
    cmd = [datasets_bin, "summary", "virus", "genome", "taxon", taxon,
           "--as-json-lines"]
    cmd.extend(extra_filters or [])
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return parse_jsonlines(proc.stdout)


def download_virus_package(taxon, dest_zip, datasets_bin="datasets",
                             extra_filters=None):
    """Run ``datasets download virus genome taxon ... --include
    cds,protein,annotation`` into ``dest_zip``."""
    cmd = [datasets_bin, "download", "virus", "genome", "taxon", taxon,
           "--include", "cds,protein,annotation",
           "--filename", str(dest_zip), "--no-progressbar"]
    cmd.extend(extra_filters or [])
    subprocess.run(cmd, check=True)


def unpack_package(zip_path, dest_dir):
    """Unzip a `datasets` virus package and return its ``data/`` dir."""
    dest_dir = Path(dest_dir)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest_dir)
    return dest_dir / "ncbi_dataset" / "data"


# --------------------------------------------------------------------------- #
# CLI entry point                                                             #
# --------------------------------------------------------------------------- #

def run(taxon, out_dir, datasets_bin="datasets", extra_filters=None,
        max_total_genomes=None, min_gene_frac=0.5, max_ambiguous_frac=0.05):
    """End-to-end: summarize -> guardrail -> download -> split -> QC -> write.

    Returns a dict summary (genome counts at each stage) for logging/tests.
    """
    extra_filters = extra_filters or []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if max_total_genomes is not None:
        summary = fetch_virus_summary(taxon, datasets_bin, extra_filters)
        if len(summary) > max_total_genomes:
            raise RuntimeError(
                f"{taxon}: {len(summary)} genomes exceed --max-total-genomes "
                f"{max_total_genomes}. Narrow the taxon with "
                f"--complete-only/--refseq/--host/--released-after, or raise "
                f"the ceiling.")

    zip_path = out_dir / "virus_package.zip"
    download_virus_package(taxon, zip_path, datasets_bin, extra_filters)
    data_dir = unpack_package(zip_path, out_dir / "unpacked")

    data_report_recs = parse_jsonlines((data_dir / "data_report.jsonl").read_text())
    annotation_recs = parse_jsonlines((data_dir / "annotation_report.jsonl").read_text())
    cds_seqs = load_fasta(data_dir / "cds.fna")
    prot_seqs = load_fasta(data_dir / "protein.faa")

    data_report_recs = dedupe_prefer_refseq(data_report_recs)
    genomes = build_genomes(data_report_recs, annotation_recs)
    n_before_qc = len(genomes)
    genomes, dropped = filter_genomes_by_gene_count(genomes, min_gene_frac)
    if dropped:
        print(f"[gather_viral_annots] dropped {len(dropped)} low-gene-count "
              f"genome(s): {', '.join(dropped)}", file=sys.stderr)

    nucl_dir = out_dir / "annots_nucls"
    prot_dir = out_dir / "annots_prots"
    written = write_genome_fastas(genomes, cds_seqs, prot_seqs, nucl_dir,
                                   prot_dir, max_ambiguous_frac)
    write_metadata_tsv(written, genomes, data_report_recs,
                        out_dir / "all_asm_acc_metadata")

    return {
        "n_accessions": len(data_report_recs),
        "n_genomes_grouped": n_before_qc,
        "n_genomes_after_qc": len(genomes),
        "n_genomes_written": len(written),
        "nucl_dir": str(nucl_dir),
        "prot_dir": str(prot_dir),
    }


def main():
    p = argparse.ArgumentParser(
        description="Download NCBI's existing CDS+protein annotations for a "
                     "viral taxon, split per genome.")
    p.add_argument("taxon", help="Taxon name or NCBI taxid, as accepted by "
                                  "`datasets ... virus genome taxon`.")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--datasets-bin", default="datasets")
    p.add_argument("--complete-only", action="store_true")
    p.add_argument("--refseq", action="store_true")
    p.add_argument("--host", default=None)
    p.add_argument("--released-after", default=None)
    p.add_argument("--geo-location", default=None)
    p.add_argument("--max-total-genomes", type=int, default=None)
    p.add_argument("--min-gene-frac", type=float, default=0.5)
    p.add_argument("--max-ambiguous-frac", type=float, default=0.05)
    args = p.parse_args()

    extra_filters = []
    if args.complete_only:
        extra_filters.append("--complete-only")
    if args.refseq:
        extra_filters.append("--refseq")
    if args.host:
        extra_filters.extend(["--host", args.host])
    if args.released_after:
        extra_filters.extend(["--released-after", args.released_after])
    if args.geo_location:
        extra_filters.extend(["--geo-location", args.geo_location])

    summary = run(args.taxon, args.out_dir, args.datasets_bin, extra_filters,
                  args.max_total_genomes, args.min_gene_frac,
                  args.max_ambiguous_frac)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
