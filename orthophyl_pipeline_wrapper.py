#!/usr/bin/env python3
"""
OrthoPhyl Pipeline Wrapper - Automated Assembly Routing and Phylogenetic Placement

REQUIRES: Python 3.7+ (uses subprocess.run with text=True)
RECOMMENDED: Run within OrthoPhyl conda environment

This wrapper orchestrates the complete pipeline:
1. Routes assemblies to appropriate pipeline (ReLeaf vs OrthoPhyl)
2. Executes ReLeaf for assemblies matching existing databases
3. For novel taxa: downloads genomes, adds queries, runs OrthoPhyl, creates databases
4. Aggregates all results

Usage:
    python orthophyl_pipeline_wrapper.py \\
        --input assemblies.tsv \\
        --database-dir databases/ \\
        --output-dir results/ \\
        --threads 32

Author: HGTool/OrthoPhyl Project
Date: 2025
"""

import os
import re
import sys
import json
import gzip
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from datetime import datetime
from collections import defaultdict

# Setup logging (will be reconfigured based on verbosity in main())
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class PipelineWrapper:
    """Main wrapper class for OrthoPhyl/ReLeaf pipeline."""

    # Valid GTDB rank letters, domain -> species.
    _GTDB_RANK_LETTERS = ('d', 'p', 'c', 'o', 'f', 'g', 's')

    # Accepted genome FASTA extensions for --genome-dir, gzip variants included.
    # Longest/most-specific first so a "*.fna.gz" file is never mistaken for "*.gz".
    _GENOME_EXTENSIONS = ('.fna.gz', '.fa.gz', '.fasta.gz', '.fna', '.fa', '.fasta')

    def __init__(
        self,
        input_file: Optional[Path] = None,
        database_dir: Path = None,
        output_dir: Optional[Path] = None,
        threads: int = 8,
        gather_script: Optional[Path] = None,
        orthophyl_runs_tsv: Optional[Path] = None,
        resume: bool = False,
        skip_download: bool = False,
        dry_run: bool = False,
        verbose: int = 0,
        low_ram: bool = False,
        use_bbmap: bool = False,
        must_keep: Optional[str] = None,
        keep_failing_query: bool = False,
        # NEW: Taxon mode parameters
        taxon: Optional[str] = None,
        taxon_rank: Optional[str] = None,
        update_existing: bool = False,
        # NEW: subclade partitioning
        max_tree_genomes: int = 2000,
        max_total_genomes: int = 5000,
        subsample_size: int = 500,
        # NEW: opt-in megatree (partition -> per-subclade trees -> backbone graft)
        megatree: bool = False,
        backbone_reps: int = 5,
        subclade_size: int = 150,
        conflict_min_support: int = 90,
        # NEW: local genome-ingest mode (build a DB from genomes already on disk,
        #   under a user-supplied clade name that was not assigned by NCBI)
        genome_dir: Optional[Path] = None,
        clade_name: Optional[str] = None,
        clade_rank: str = 'g',
        clade_taxonomy: Optional[str] = None,
        skip_qc: bool = False,
    ):
        # Validate mutually exclusive mode flags (--input / --taxon / --genome-dir)
        modes_given = sum(bool(x) for x in (input_file, taxon, genome_dir))
        if modes_given > 1:
            if input_file and taxon:
                raise ValueError("Cannot specify both --input and --taxon. Use one or the other.")
            raise ValueError(
                "Cannot specify more than one of --input, --taxon, --genome-dir. "
                "Use exactly one.")

        self.input_file = Path(input_file) if input_file else None
        self.database_dir = Path(database_dir) if database_dir else None

        # Default output_dir to database_dir/.pipeline_runs/<name>_<timestamp> if not provided
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            run_name = self._default_run_name(taxon or clade_name)
            self.output_dir = self.database_dir / '.pipeline_runs' / f'{run_name}_{timestamp}'
            logger.info(f"No --output-dir provided, using: {self.output_dir}")
        
        self.threads = threads
        self.gather_script = Path(gather_script) if gather_script else None
        self.orthophyl_runs_tsv = Path(orthophyl_runs_tsv) if orthophyl_runs_tsv else None
        self.resume = resume
        self.skip_download = skip_download
        self.dry_run = dry_run
        self.verbose = verbose
        self.low_ram = low_ram
        self.use_bbmap = use_bbmap
        # Genome-retention controls, passed through to gather_filter_asms.sh.
        #   must_keep : comma-sep accession list OR path to a file (one per line);
        #               these must survive QC or the download aborts.
        #   keep_failing_query : let query genomes that fail QC through with a
        #               warning instead of aborting.
        self.must_keep = must_keep
        self.keep_failing_query = keep_failing_query

        # NEW: Taxon mode
        self.taxon = taxon
        self.taxon_rank = taxon_rank
        self.update_existing = update_existing
        self.taxon_mode = taxon is not None

        # NEW: Local genome-ingest mode. Builds a tree + database from genomes
        #   already on disk under a user-supplied clade name/taxonomy that was
        #   NOT assigned by NCBI (unless --clade-name happens to resolve against
        #   the local taxdump -- see _resolve_local_taxonomy).
        self.local_mode = genome_dir is not None
        if self.local_mode:
            self.genome_dir = Path(genome_dir)
            if not self.genome_dir.exists():
                raise ValueError(f"--genome-dir does not exist: {self.genome_dir}")
            if not clade_name or not clade_name.strip():
                raise ValueError("--clade-name is required with --genome-dir.")
            if clade_rank not in self._GTDB_RANK_LETTERS:
                raise ValueError(
                    f"--clade-rank must be one of {self._GTDB_RANK_LETTERS}, "
                    f"got: {clade_rank!r}")
        else:
            self.genome_dir = None
        self.clade_name = clade_name
        self.clade_rank = clade_rank
        self.clade_taxonomy = clade_taxonomy
        self.skip_qc = skip_qc

        # NEW: large-taxon handling. When a taxon's RAW downloaded genome set
        #   exceeds max_tree_genomes, the DEFAULT behavior is to build ONE tree from
        #   a diverse MASH subsample of subsample_size genomes (greedy max-min, no
        #   O(n^2) matrix -- see python_scripts/subsample_genomes.py). The per-
        #   subclade partition/megatree path (subclade_partition.py) is opt-in and
        #   still guarded by max_total_genomes below.
        self.max_tree_genomes = max_tree_genomes
        self.subsample_size = subsample_size

        # Guardrail for the PARTITION/megatree path only: subclade_partition.py runs
        #   an all-vs-all `mash triangle` and builds a DENSE NxN distance matrix,
        #   which is O(n^2) in time and memory -- e.g. ~20 GB just for the matrix at
        #   n=50k. Above this ceiling we refuse to partition rather than OOM-kill the
        #   node. The default subsample path does NOT hit this (it never builds the
        #   matrix); this only bounds the opt-in partition/megatree route.
        self.max_total_genomes = max_total_genomes

        # Opt-in MEGATREE path: instead of subsampling an oversized taxon down to
        #   one tree, partition the raw set into size-bounded subclades
        #   (subclade_size each), build a full tree per subclade, build a small
        #   BACKBONE tree from backbone_reps diverse reps per subclade, and GRAFT
        #   each subclade tree onto its reps in the backbone -> one merged tree with
        #   every genome. High-support bipartition disagreements between a subclade
        #   tree and the backbone are FLAGGED (>= conflict_min_support), not
        #   resolved. This path builds the dense matrix and so IS subject to
        #   max_total_genomes above.
        self.megatree = megatree
        self.backbone_reps = backbone_reps
        self.subclade_size = subclade_size
        self.conflict_min_support = conflict_min_support

        # Script paths (relative to this wrapper)
        self.script_dir = Path(__file__).parent
        self.assembly_router = self.script_dir / "assembly_router" / "assembly_router.py"
        self.database_creator = self.script_dir / "assembly_router" / "create_hierarchical_database.py"
        self.releaf_versioner = self.script_dir / "assembly_router" / "add_releaf_version.py"
        self.subclade_partitioner = self.script_dir / "python_scripts" / "subclade_partition.py"
        self.subsampler = self.script_dir / "python_scripts" / "subsample_genomes.py"
        self.megatree_grafter = self.script_dir / "python_scripts" / "megatree_graft.py"
        self.orthophyl_script = self.script_dir / "OrthoPhyl.sh"
        self.releaf_script = self.script_dir / "ReLeaf.sh"
        
        # Output subdirectories
        self.routing_dir = self.output_dir / "00_routing"
        self.releaf_dir = self.output_dir / "01_releaf_only"
        self.orthophyl_dir = self.output_dir / "02_orthophyl_novel"
        self.results_dir = self.output_dir / "03_results"
        self.logs_dir = self.output_dir / "logs"
        
        # Checkpoint tracking
        self.checkpoint_dir = self.output_dir / "checkpoints"
        
        # Status tracking with success/failure lists
        self.pipeline_status = {
            'start_time': datetime.now().isoformat(),
            'phases': {},
            'summary': {},
            'successes': [],
            'failures': []
        }

    @staticmethod
    def _default_run_name(taxon: Optional[str]) -> str:
        """Build a filesystem-safe run name from the taxon.

        Numeric taxa (NCBI TaxIDs) become ``TaxID<num>``; named taxa are
        sanitized to alphanumerics/underscores. Batch mode (no taxon) falls
        back to ``run``.
        """
        if not taxon:
            return 'run'
        taxon = taxon.strip()
        if taxon.isdigit():
            return f'TaxID{taxon}'
        # Sanitize the name for use as a directory component
        safe = re.sub(r'[^A-Za-z0-9._-]+', '_', taxon).strip('_')
        return safe or 'run'

    @staticmethod
    def _validate_clade_taxonomy(taxonomy: str, clade_name: str) -> str:
        """Validate a verbatim --clade-taxonomy string.

        Must be semicolon-separated rank__name tokens (e.g. "d__Bacteria;p__...")
        using only the GTDB rank letters, and must actually mention clade_name
        somewhere in it (a copy-paste of the wrong taxon's string is otherwise a
        silent, hard-to-notice mistake). Raises ValueError on malformed input.
        """
        taxonomy = taxonomy.strip()
        if not taxonomy:
            raise ValueError("--clade-taxonomy must not be empty.")
        tokens = taxonomy.split(';')
        for tok in tokens:
            if '__' not in tok:
                raise ValueError(
                    f"--clade-taxonomy token {tok!r} is not in 'rank__name' form. "
                    f"Expected e.g. 'd__Bacteria;p__Pseudomonadota;...;g__{clade_name}'.")
            prefix, _, _name = tok.partition('__')
            if prefix not in PipelineWrapper._GTDB_RANK_LETTERS:
                raise ValueError(
                    f"--clade-taxonomy token {tok!r} has an invalid rank letter "
                    f"{prefix!r}; must be one of {PipelineWrapper._GTDB_RANK_LETTERS}.")
        if clade_name.lower() not in taxonomy.lower():
            raise ValueError(
                f"--clade-taxonomy {taxonomy!r} does not mention --clade-name "
                f"{clade_name!r}. Double-check it's the right string.")
        return taxonomy

    def _resolve_local_taxonomy(self) -> Tuple[str, bool]:
        """Resolve the taxonomy string for local genome-ingest mode.

        Returns (taxonomy, is_routable) by precedence:
          1. --clade-taxonomy given: validated and used verbatim. Routable.
          2. --clade-name resolves against the local taxdump (a real NCBI taxon):
             render its full lineage. Routable.
          3. Unresolvable: fall back to "<clade-rank>__<clade-name>" and warn
             that the database will build but won't be matched by fully
             specified queries.

        Best-effort: taxdump absence/download failure falls back to step 3
        rather than failing the run (an offline --clade-taxonomy user needs no
        taxdump at all).
        """
        if self.clade_taxonomy:
            taxonomy = self._validate_clade_taxonomy(self.clade_taxonomy, self.clade_name)
            logger.info(f"  Using --clade-taxonomy verbatim: {taxonomy}")
            return taxonomy, True

        # Try resolving --clade-name against the local taxdump.
        try:
            sys.path.insert(0, str(self.script_dir / "utils"))
            from taxon_assembly_gatherer import (
                NCBITaxonomy, render_gtdb_lineage, TaxonAssemblyGatherer)

            taxdump_dir = self.output_dir / "taxon_query" / "taxdump"
            # Reuse the gatherer's download-if-missing logic without invoking its
            # __init__ (which raises on an unresolvable taxon -- not an error here).
            stub = TaxonAssemblyGatherer.__new__(TaxonAssemblyGatherer)
            stub.output_dir = self.output_dir / "taxon_query"
            stub.output_dir.mkdir(parents=True, exist_ok=True)
            stub.taxdump_dir = taxdump_dir
            stub._ensure_taxonomy_database()

            taxonomy_db = NCBITaxonomy(taxdump_dir)
            taxid = taxonomy_db.resolve_taxon(self.clade_name)
            if taxid:
                lineage = taxonomy_db.get_lineage(taxid)
                rendered = render_gtdb_lineage(lineage)
                if rendered:
                    rank = taxonomy_db.get_rank(taxid)
                    logger.info(
                        f"  ✓ --clade-name '{self.clade_name}' resolved against NCBI "
                        f"taxdump (TaxID {taxid}, rank {rank}): {rendered}")
                    return rendered, True
        except Exception as e:
            logger.warning(f"  ⚠ Could not resolve --clade-name against local taxdump: {e}")

        # Unresolvable: fall back to a name-only taxonomy at --clade-rank.
        taxonomy = f"{self.clade_rank}__{self.clade_name}"
        template = (
            f'd__Bacteria;p__...;c__...;o__...;f__...;g__{self.clade_name}'
            if self.clade_rank == 'g' else
            f'd__...;...;{self.clade_rank}__{self.clade_name}')
        logger.warning(
            f"  ⚠ '{self.clade_name}' did not resolve to a known NCBI taxon. "
            f"Falling back to name-only taxonomy: {taxonomy!r}")
        logger.warning(
            "  ⚠ This database WILL build, but will NOT be matched by "
            "fully-specified taxonomy queries (routing compares every rank "
            "from domain down). If you know the real lineage, pass it "
            f'explicitly, e.g.: --clade-taxonomy "{template}"')
        return taxonomy, False

    def run(self):
        """Main execution pipeline."""
        try:
            logger.info("=" * 70)
            logger.info("ORTHOPHYL PIPELINE WRAPPER")
            if self.taxon_mode:
                logger.info(f"*** TAXON MODE: {self.taxon} ***")
                if self.update_existing:
                    logger.info("*** UPDATE MODE: Checking for new assemblies ***")
            if self.local_mode:
                logger.info(f"*** LOCAL GENOME-INGEST MODE: {self.clade_name} ***")
                logger.info("*** NOTE: clade name/taxonomy is user-supplied, not "
                            "assigned by NCBI ***")
            if self.dry_run:
                logger.info("*** DRY RUN MODE - No commands will be executed ***")
            if self.verbose == 1:
                logger.info("*** VERBOSE MODE (Level 1) - Showing stdout ***")
            elif self.verbose >= 2:
                logger.info("*** VERBOSE MODE (Level 2) - Showing stdout and stderr ***")
            logger.info("=" * 70)
            
            if self.verbose:
                logger.info(f"Configuration:")
                if self.taxon_mode:
                    logger.info(f"  Taxon: {self.taxon}")
                    logger.info(f"  Taxon rank: {self.taxon_rank or 'auto-detect'}")
                    logger.info(f"  Update existing: {self.update_existing}")
                else:
                    logger.info(f"  Input file: {self.input_file}")
                logger.info(f"  Database dir: {self.database_dir}")
                logger.info(f"  Output dir: {self.output_dir}")
                logger.info(f"  Threads: {self.threads}")
                logger.info(f"  Gather script: {self.gather_script}")
                logger.info(f"  Resume: {self.resume}")
                logger.info(f"  Skip download: {self.skip_download}")
                logger.info(f"  Low RAM mode: {self.low_ram}")
                logger.info(f"  Use bbmap stats: {self.use_bbmap}")
            
            # Phase 1: Initialization
            self._phase_initialization()
            
            # Branch based on mode
            if self.taxon_mode:
                # NEW: Taxon mode workflow
                return self._run_taxon_mode()
            elif self.local_mode:
                # NEW: Local genome-ingest mode workflow
                return self._run_local_genomes_mode()
            else:
                # Original: Batch mode workflow
                return self._run_batch_mode()
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}", exc_info=True)
            self.pipeline_status['status'] = 'failed'
            self.pipeline_status['error'] = str(e)
            self._save_final_status()
            return 1
    
    def _run_batch_mode(self) -> int:
        """Run original batch mode workflow."""
        # Phase 2: Routing
        routing_results = self._phase_routing()
        
        # Phase 3a: ReLeaf route
        if routing_results['releaf_batch']:
            self._phase_releaf(routing_results['releaf_batch'])
        
        # Phase 3b: OrthoPhyl route
        if routing_results['orthophyl_batch']:
            self._phase_orthophyl(routing_results['orthophyl_batch'])

        # Phase 4: Results aggregation
        self._phase_aggregation()
        
        # Check for failures
        failures = len(self.pipeline_status.get('failures', []))
        
        logger.info("=" * 70)
        if failures == 0:
            logger.info("PIPELINE COMPLETE!")
        else:
            logger.error(f"PIPELINE COMPLETED WITH {failures} FAILURE(S)")
            logger.error("See summary report for details")
        logger.info("=" * 70)
        
        self._save_final_status()
        return 1 if failures > 0 else 0
    
    def _run_taxon_mode(self) -> int:
        """Run taxon mode workflow."""
        # Check if database exists for this taxon
        existing_db = self._check_existing_taxon_database()
        
        if existing_db and self.update_existing:
            # Update mode: add new assemblies to existing database
            logger.info(f"\n✓ Found existing database: {existing_db['db_dir']}")
            logger.info(f"  Current assemblies: {existing_db['n_assemblies']}")
            return self._run_taxon_update_mode(existing_db)
        elif existing_db and not self.update_existing:
            # Database exists but not in update mode
            logger.error(f"\n✗ Database already exists for taxon '{self.taxon}': {existing_db['db_dir']}")
            logger.error(f"  Use --update-existing to add new assemblies to this database")
            logger.error(f"  Or use a different --output-dir to create a new run")
            return 1
        else:
            # Create new database from taxon
            logger.info(f"\n→ No existing database found for taxon '{self.taxon}'")
            logger.info(f"  Creating new database...")
            return self._run_taxon_create_mode()
    
    def _phase_initialization(self):
        """Phase 1: Initialize directory structure and validate dependencies."""
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 1: INITIALIZATION")
        logger.info("=" * 70)
        
        if self._check_checkpoint('initialization') and self.resume:
            logger.info("✓ Initialization already complete (resuming)")
            return
        
        # Create directory structure
        for dir_path in [self.routing_dir, self.releaf_dir, self.orthophyl_dir,
                        self.results_dir, self.logs_dir, self.checkpoint_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        logger.info("✓ Created directory structure")
        
        # Validate dependencies
        self._validate_dependencies()
        logger.info("✓ Dependencies validated")
        
        # Initialize or validate databases (skip in taxon create mode / local mode)
        if self.local_mode:
            # Local genome-ingest mode: database directory will be created during
            # workflow, same as taxon create mode. No database_index.json required.
            logger.info(f"✓ Local genome-ingest mode: database will be created for "
                        f"'{self.clade_name}'")
            if self.database_dir:
                self.database_dir.mkdir(parents=True, exist_ok=True)
        elif self.taxon_mode and not self.update_existing:
            # Taxon create mode: database directory will be created during workflow
            logger.info(f"✓ Taxon create mode: database will be created for '{self.taxon}'")
            # Ensure database_dir exists as a directory (but may be empty)
            if self.database_dir:
                self.database_dir.mkdir(parents=True, exist_ok=True)
        elif self.taxon_mode and self.update_existing:
            # Taxon update mode: database directory must exist with databases
            logger.info(f"✓ Taxon update mode: checking for existing database for '{self.taxon}'")
            if not self.database_dir.exists():
                raise FileNotFoundError(
                    f"Database directory not found: {self.database_dir}\n"
                    f"Cannot update non-existent database. Use create mode first."
                )
            # Database index is optional in taxon mode - we search for matching databases directly
        elif not self.database_dir.exists() or not (self.database_dir / "database_index.json").exists():
            if self.orthophyl_runs_tsv and self.orthophyl_runs_tsv.exists():
                logger.info("Creating initial databases...")
                self._create_initial_databases()
            else:
                raise FileNotFoundError(
                    f"Database directory not found: {self.database_dir}\n"
                    f"Please provide --orthophyl-runs to create initial databases"
                )
        else:
            logger.info(f"✓ Found existing database directory: {self.database_dir}")
            with open(self.database_dir / "database_index.json", 'r') as f:
                db_index = json.load(f)
            logger.info(f"  Loaded {db_index['n_databases']} databases")
        
        self._write_checkpoint('initialization')
        self.pipeline_status['phases']['initialization'] = {'status': 'complete'}
    
    def _validate_dependencies(self):
        """Check that all required scripts and tools are available."""
        required_scripts = {
            'Assembly Router': self.assembly_router,
            'Database Creator': self.database_creator,
            'OrthoPhyl': self.orthophyl_script,
            'ReLeaf': self.releaf_script
        }
        
        missing = []
        for name, script_path in required_scripts.items():
            if not script_path.exists():
                missing.append(f"{name}: {script_path}")
        
        if missing:
            raise FileNotFoundError(
                "Missing required scripts:\n" + "\n".join(f"  - {m}" for m in missing)
            )
        
        # Check if gather script is provided and exists
        if self.gather_script and not self.gather_script.exists():
            logger.warning(f"Genome download script not found: {self.gather_script}")
            logger.warning("  Will generate manual download instructions instead")
            self.gather_script = None

    def _log_command(self, cmd: List) -> None:
        """Print a constructed command to real STDOUT, unconditionally.

        Every subprocess the wrapper runs (or skips under --dry-run) is
        announced here rather than via logger.info, because logger's default
        StreamHandler writes to stderr and most call sites previously gated
        this line behind --verbose. Call this right after `cmd` is fully
        built, before the --dry-run branch, so dry runs show the same line
        real runs do.
        """
        tag = "[DRY RUN] " if self.dry_run else ""
        print(f"{tag}COMMAND: {' '.join(str(c) for c in cmd)}", flush=True)

    def _create_initial_databases(self):
        """Create initial databases from orthophyl_runs.tsv."""
        cmd = [
            'python', str(self.database_creator),
            '--input', str(self.orthophyl_runs_tsv),
            '--output-dir', str(self.database_dir)
        ]

        self._log_command(cmd)

        if self.dry_run:
            logger.info("  [DRY RUN] Would create initial databases")
            return
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Database creation failed:\n{result.stderr}")
        
        logger.info("✓ Initial databases created")
    
    def _phase_routing(self) -> Dict:
        """Phase 2: Route all assemblies."""
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 2: ASSEMBLY ROUTING")
        logger.info("=" * 70)
        
        if self._check_checkpoint('routing') and self.resume:
            logger.info("✓ Routing already complete (resuming)")
            return self._load_routing_results()
        
        # Build command
        cmd = [
            'python', str(self.assembly_router),
            '--batch', str(self.input_file),
            '--database-dir', str(self.database_dir),
            '--output-dir', str(self.routing_dir),
            '--threads', str(self.threads)
        ]
        
        if self.gather_script:
            cmd.extend(['--gather-filter-script', str(self.gather_script)])

        logger.info(f"Running assembly router...")
        self._log_command(cmd)
        logger.info(f"  Input: {self.input_file}")
        logger.info(f"  Database: {self.database_dir}")

        if self.dry_run:
            logger.info("  [DRY RUN] Would run assembly routing")
            # In dry run, create mock routing results for preview
            return self._create_mock_routing_results()
        
        # Run routing
        log_file = self.logs_dir / "routing.log"
        with open(log_file, 'w') as f:
            if self.verbose == 1:
                result = subprocess.run(cmd, stderr=f, text=True)
            elif self.verbose >= 2:
                result = subprocess.run(cmd, text=True)
            else:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Routing failed. Check log: {log_file}")
        
        logger.info("✓ Routing complete")
        
        # Parse routing results
        routing_results = self._parse_routing_results()

        logger.info(f"\nRouting Summary:")
        logger.info(f"  ReLeaf route: {len(routing_results['releaf_batch'])} assemblies")
        logger.info(f"  OrthoPhyl route: {sum(len(v) for v in routing_results['orthophyl_batch'].values())} assemblies")
        logger.info(f"    ({len(routing_results['orthophyl_batch'])} unique taxa)")

        self._write_checkpoint('routing')
        self.pipeline_status['phases']['routing'] = {
            'status': 'complete',
            'releaf_count': len(routing_results['releaf_batch']),
            'orthophyl_count': sum(len(v) for v in routing_results['orthophyl_batch'].values()),
        }

        return routing_results

    def _parse_routing_results(self) -> Dict:
        """Parse routing decision JSON files."""
        releaf_batch = []
        orthophyl_batch = defaultdict(list)

        for json_file in self.routing_dir.glob("routing_decision_*.json"):
            with open(json_file, 'r') as f:
                decision = json.load(f)

            pipeline = decision['pipeline']
            if pipeline == 'ReLeaf':
                releaf_batch.append({
                    'assembly_id': decision['assembly_id'],
                    'assembly_path': decision['assembly'],
                    'database': decision['matched_database'],
                    'database_dir': decision['database_dir'],
                    'tree_method': decision.get('tree_method', 'iqtree'),
                    'tree_data': decision.get('tree_data', 'CDS')
                })
            else:  # OrthoPhyl
                taxon = decision['download_value']
                orthophyl_batch[taxon].append({
                    'assembly_id': decision['assembly_id'],
                    'assembly_path': decision['assembly'],
                    'download_rank': decision['download_rank'],
                    'taxonomy': decision['query_taxonomy'],
                    'download_taxonomy': decision['download_taxonomy']
                })

        return {
            'releaf_batch': releaf_batch,
            'orthophyl_batch': dict(orthophyl_batch),
        }
    
    def _phase_releaf(self, releaf_batch: List[Dict]):
        """Phase 3a: Execute ReLeaf for matched assemblies."""
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 3A: RELEAF ROUTE (Matched Databases)")
        logger.info("=" * 70)
        
        # Group by database
        by_database = defaultdict(list)
        for item in releaf_batch:
            by_database[item['database']].append(item)
        
        logger.info(f"Processing {len(releaf_batch)} assemblies across {len(by_database)} databases")
        
        for database_name, assemblies in by_database.items():
            checkpoint_name = f"releaf_{database_name}"
            
            if self._check_checkpoint(checkpoint_name) and self.resume:
                logger.info(f"\n✓ ReLeaf for {database_name} already complete (resuming)")
                continue
            
            logger.info(f"\n--- Processing database: {database_name} ({len(assemblies)} assemblies) ---")
            
            # Prepare input directory
            input_dir = self.releaf_dir / database_name / "input_genomes"
            input_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy assemblies
            for asm in assemblies:
                src = Path(asm['assembly_path'])
                dst = input_dir / f"{asm['assembly_id']}.fna"
                if not dst.exists():
                    shutil.copy(src, dst)
                logger.info(f"  Prepared: {asm['assembly_id']}")
            
            # Get database path and parameters
            db_dir = Path(assemblies[0]['database_dir'])
            tree_method = assemblies[0]['tree_method']
            tree_data = assemblies[0]['tree_data']
            
            # Run ReLeaf
            output_dir = self.releaf_dir / database_name
            try:
                self._run_releaf(
                    database_dir=db_dir,
                    input_genomes=input_dir,
                    output_dir=output_dir,
                    tree_method=tree_method,
                    tree_data=tree_data,
                    database_name=database_name,
                    n_assemblies=len(assemblies)
                )
                self._write_checkpoint(checkpoint_name)
            except Exception as e:
                logger.error(f"✗ ReLeaf failed for {database_name}: {e}")
                # Continue with other databases rather than failing entire pipeline
                continue
        
        self.pipeline_status['phases']['releaf'] = {
            'status': 'complete',
            'databases_processed': len(by_database)
        }
    
    def _run_releaf(
        self,
        database_dir: Path,
        input_genomes: Path,
        output_dir: Path,
        tree_method: str,
        tree_data: str,
        database_name: str,
        n_assemblies: int = 0
    ):
        """Execute ReLeaf for a single database."""
        # Handle stale ReLeaf_dir from previous run
        releaf_dir = database_dir / 'orthophyl_run' / 'ReLeaf_dir'
        if releaf_dir.exists() and not self.resume:
            logger.warning(f"  ⚠ Removing stale ReLeaf output: {releaf_dir}")
            shutil.rmtree(releaf_dir)
        
        # Fix: Use correct ReLeaf.sh flags (see script_lib/arg_parse_addem.sh)
        # -s/--storage_dir, -g/--genome_dir, -t/--threads, -p/--phylo_tool, -o/--omics
        cmd = [
            str(self.releaf_script),
            '-s', str(database_dir / 'orthophyl_run'),  # --storage_dir (not --store)
            '-g', str(input_genomes),                    # --genome_dir (not --input_genomes)
            '-t', str(self.threads),                     # --threads
            '-p', tree_method,                           # --phylo_tool (not --tree_method)
            '-o', tree_data                              # --omics (not --TREE_DATA, expects CDS|PROT|BOTH)
        ]
        
        logger.info(f"  Running ReLeaf...")
        self._log_command(cmd)
        logger.info(f"    Database: {database_dir}")
        logger.info(f"    Method: {tree_method}, Data: {tree_data}")

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would run ReLeaf for {database_name}")
            return

        log_file = self.logs_dir / f"releaf_{database_name}.log"
        with open(log_file, 'w') as f:
            if self.verbose == 1:
                result = subprocess.run(
                    cmd,
                    stderr=f,
                    text=True,
                    cwd=str(output_dir)
                )
            elif self.verbose >= 2:
                result = subprocess.run(
                    cmd,
                    text=True,
                    cwd=str(output_dir)
                )
            else:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    text=True,
                    cwd=str(output_dir)
                )
        
        if result.returncode != 0:
            error_msg = f"ReLeaf failed for {database_name}. Check log: {log_file}"
            self.pipeline_status['failures'].append({
                'type': 'releaf',
                'database': database_name,
                'error': error_msg,
                'log': str(log_file)
            })
            raise RuntimeError(error_msg)
        
        # Verify expected outputs exist (ReLeaf writes to $store/ReLeaf_dir)
        releaf_output = database_dir / 'orthophyl_run' / 'ReLeaf_dir'
        expected_files = [
            releaf_output / 'new_prot_alignments.trm.nm',
            releaf_output / 'new_CDS_alignments.trm.nm',
            releaf_output / 'new_trees'
        ]
        missing = [f for f in expected_files if not f.exists()]
        if missing:
            error_msg = (
                f"ReLeaf completed but missing expected outputs:\n" +
                "\n".join(f"  - {f}" for f in missing) +
                f"\nCheck log: {log_file}"
            )
            self.pipeline_status['failures'].append({
                'type': 'releaf',
                'database': database_name,
                'error': error_msg,
                'log': str(log_file)
            })
            raise RuntimeError(error_msg)
        
        logger.info(f"  ✓ ReLeaf complete for {database_name}")
        
        # Track success
        self.pipeline_status['successes'].append({
            'type': 'releaf',
            'database': database_name,
            'assemblies': n_assemblies
        })
            
        # Create new database version from ReLeaf output
        if self.releaf_versioner.exists():
            self._create_releaf_version(database_name, database_dir)
        else:
            logger.warning(f"  ⚠ ReLeaf versioner not found, skipping version creation")
    
    def _create_releaf_version(self, database_name: str, database_dir: Path):
        """Create a new database version from ReLeaf output."""
        logger.info(f"\n  Creating new database version from ReLeaf output...")
        
        # ReLeaf writes to $store/ReLeaf_dir (inside database's orthophyl_run)
        releaf_output_dir = database_dir / 'orthophyl_run' / 'ReLeaf_dir'
        
        if not releaf_output_dir.exists():
            logger.warning(f"  ⚠ ReLeaf output not found: {releaf_output_dir}")
            return
        
        # Find the database directory
        db_dir = database_dir

        cmd = [
            'python', str(self.releaf_versioner),
            '--database-dir', str(db_dir),
            '--releaf-output', str(releaf_output_dir)
        ]
        
        self._log_command(cmd)

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would create new database version")
            return

        log_file = self.logs_dir / f"releaf_version_{database_name}.log"
        with open(log_file, 'w') as f:
            if self.verbose == 1:
                result = subprocess.run(cmd, stderr=f, text=True)
            elif self.verbose >= 2:
                result = subprocess.run(cmd, text=True)
            else:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
        
        if result.returncode != 0:
            logger.error(f"  ✗ Failed to create database version. Check log: {log_file}")
        else:
            logger.info(f"  ✓ Created new database version for {database_name}")
    
    def _phase_orthophyl(self, orthophyl_batch: Dict[str, List[Dict]]):
        """Phase 3b: Execute OrthoPhyl route for novel taxa."""
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 3B: ORTHOPHYL ROUTE (Novel Taxa)")
        logger.info("=" * 70)
        
        logger.info(f"Processing {len(orthophyl_batch)} novel taxa")
        
        for taxon_name, assemblies in orthophyl_batch.items():
            logger.info(f"\n{'=' * 60}")
            logger.info(f"Processing taxon: {taxon_name} ({len(assemblies)} assemblies)")
            logger.info(f"{'=' * 60}")
            try:
                self._process_orthophyl_taxon(taxon_name, assemblies)
            except Exception as e:
                logger.error(f"✗ OrthoPhyl failed for {taxon_name}: {e}")
                # Continue with other taxa rather than failing entire pipeline
                continue

        self.pipeline_status['phases']['orthophyl'] = {
            'status': 'complete',
            'taxa_processed': len(orthophyl_batch)
        }

    def _process_orthophyl_taxon(self, taxon_name: str, assemblies: List[Dict]):
        """Handle one novel taxon: download RAW -> partition -> QC+build needed subclades.

        Partitioning happens BEFORE the expensive CheckM2 QC: only the subclade(s)
        that will actually become trees get QC'd. When the raw set fits under
        max_tree_genomes, this collapses to the classic single-tree flow.
        """
        download_dir = self.orthophyl_dir / "downloads" / taxon_name

        # ---- Stage 1: download the RAW genome set (no QC yet) ----
        raw_dir = download_dir / "assemblies_all.TMP"
        if self._check_checkpoint(f"download_{taxon_name}") and self.resume and raw_dir.exists():
            logger.info(f"  ✓ Raw download already complete (resuming)")
        elif not self.skip_download:
            raw_dir = self._download_raw(taxon_name, download_dir, query_assemblies=assemblies)
            self._write_checkpoint(f"download_{taxon_name}")
        else:
            logger.info(f"  Skipping download (--skip-download)")
            if not raw_dir.exists():
                raise FileNotFoundError(
                    f"--skip-download specified but no raw genomes at: {raw_dir}\n"
                    f"Provide raw FASTAs there or drop --skip-download.")

        raw_files = (list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))
                     if raw_dir.exists() else [])
        raw_count = len(raw_files)
        logger.info(f"  Raw genomes downloaded: {raw_count}")

        # ---- Opt-in megatree: partition -> per-subclade trees -> backbone graft ----
        # When --megatree is set, an oversized taxon is covered in full instead of
        # subsampled. Under the ceiling this collapses to the normal single tree.
        if self.megatree and raw_count > self.max_tree_genomes:
            self._run_megatree(
                taxon_name=taxon_name, raw_dir=raw_dir,
                query_assemblies=assemblies,
                taxonomy=assemblies[0]['download_taxonomy'])
            return

        # ---- Stage 2: cap tree size (default = diverse subsample) ----
        # Over the ceiling, the DEFAULT large-taxon behavior is to build ONE tree
        # from a MASH greedy max-min diverse subset (no O(n^2) matrix). Query
        # genomes seed the pick so they are always retained. The subsampled dir
        # then flows through the single-tree path below exactly like a raw dir.
        if raw_count > self.max_tree_genomes:
            logger.info(f"  Raw count {raw_count} > max_tree_genomes "
                        f"{self.max_tree_genomes}: diverse-subsampling to "
                        f"{self.subsample_size} genomes")
            query_stems = [Path(a['assembly_path']).stem for a in assemblies]
            raw_dir = self._subsample_genomes(
                taxon_name, raw_dir, self.subsample_size,
                must_keep_stems=query_stems)
            raw_files = (list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))
                         if raw_dir.exists() else [])
            raw_count = len(raw_files)

        # ---- Stage 3: QC + build one tree over the (possibly subsampled) set ----
        entry = {'subclade_id': 1, 'name': taxon_name, 'n_genomes': raw_count,
                 'members_file': None, 'sketch_file': None}
        self._build_subclade(
            taxon_name=taxon_name, entry=entry, raw_dir=raw_dir,
            query_assemblies=assemblies,
            taxonomy=assemblies[0]['download_taxonomy'],
            is_subclade=False)

    def _run_megatree(self, taxon_name: str, raw_dir: Path,
                      query_assemblies: List[Dict], taxonomy: str) -> None:
        """Opt-in large-taxon strategy: partition -> per-subclade trees -> graft.

        For an oversized taxon (raw count > max_tree_genomes) build FULL coverage
        rather than a subsample:

          1. Enforce --max-total-genomes (the partitioner builds a dense O(n^2)
             MASH matrix; refuse rather than OOM).
          2. Partition the raw set into subclades of <= subclade_size genomes.
          3. Build a full OrthoPhyl tree for EVERY subclade (queries mapped to
             their subclade via the manifest's query_assignments).
          4. Pick backbone_reps diverse reps per subclade (MASH greedy max-min,
             seeded by that subclade's queries), pool them, and build one BACKBONE
             OrthoPhyl tree.
          5. Graft each subclade tree onto its reps in the backbone -> one merged
             megatree, flagging (not resolving) high-support bipartition conflicts.
          6. Publish the merged tree + conflict report and create the taxon DB from
             the backbone run so ReLeaf has a coherent HMM set.

        Checkpointed per stage; dry-run short-circuits.
        """
        logger.info("\n" + "=" * 70)
        logger.info(f"MEGATREE: full-coverage build for {taxon_name}")
        logger.info("=" * 70)

        raw_files = (list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))
                     if raw_dir.exists() else [])
        raw_count = len(raw_files)

        # (1) Guard the dense matrix.
        self._enforce_total_genome_ceiling(taxon_name, raw_count)

        # (2) Partition into size-bounded subclades.
        manifest = self._partition_genomes(
            taxon_name, raw_dir, query_assemblies, max_size=self.subclade_size)
        subclades = manifest['subclades']
        assignments = manifest.get('query_assignments', {})
        logger.info(f"  Partitioned into {len(subclades)} subclade(s)")

        # Map queries to their subclade.
        query_by_subclade: Dict[str, List[Dict]] = {}
        for asm in query_assemblies:
            name = Path(asm['assembly_path']).name
            sc = assignments.get(name)
            query_by_subclade.setdefault(sc, []).append(asm)

        # (3) Build a full tree for every subclade + (4) collect backbone reps.
        backbone_dir = self.orthophyl_dir / "megatree" / taxon_name / "backbone_genomes"
        backbone_dir.mkdir(parents=True, exist_ok=True)
        # subclade name -> {tree, reps} accumulated for the graft.
        subclade_specs: Dict[str, Dict] = {}

        for entry in subclades:
            sc_name = entry['name']
            queries_here = query_by_subclade.get(sc_name, [])
            logger.info(f"\n  Building subclade {sc_name} "
                        f"({entry['n_genomes']} raw genomes, "
                        f"{len(queries_here)} query)")
            self._build_subclade(
                taxon_name=taxon_name, entry=entry, raw_dir=raw_dir,
                query_assemblies=queries_here, taxonomy=taxonomy,
                is_subclade=True)

            # Backbone reps: diverse pick over this subclade's QC-kept genomes,
            # seeded by its queries so they anchor the backbone.
            sc_genomes = self.orthophyl_dir / "downloads" / sc_name / "genomes_to_keep"
            seed_stems = [Path(a['assembly_path']).stem for a in queries_here]
            reps_dir = self._subsample_genomes(
                f"{sc_name}_backbone", sc_genomes, self.backbone_reps,
                must_keep_stems=seed_stems)
            rep_stems = []
            if not self.dry_run:
                for p in (list(reps_dir.glob("*.fna")) +
                          list(reps_dir.glob("*.fasta"))):
                    rep_stems.append(p.stem)
                    dst = backbone_dir / p.name
                    if not dst.exists() and not dst.is_symlink():
                        try:
                            os.symlink(os.path.abspath(p), dst)
                        except OSError:
                            shutil.copy(p, dst)
            sc_tree = self._locate_species_tree(
                self.orthophyl_dir / "orthophyl_runs" / sc_name)
            subclade_specs[sc_name] = {'tree': sc_tree, 'reps': rep_stems}

        # (5) Build the backbone tree over pooled reps.
        backbone_out = self.orthophyl_dir / "megatree" / taxon_name / "backbone_run"
        if self._check_checkpoint(f"megatree_backbone_{taxon_name}") and self.resume:
            logger.info(f"  ✓ Backbone tree already built for {taxon_name} (resuming)")
        else:
            logger.info(f"\n  Building backbone tree from "
                        f"{len(subclade_specs)} subclade rep sets")
            self._run_orthophyl(
                input_dir=backbone_dir, output_dir=backbone_out,
                taxon_name=f"{taxon_name}_backbone", assemblies=[])
            self._write_checkpoint(f"megatree_backbone_{taxon_name}")
        backbone_tree = self._locate_species_tree(backbone_out)

        # (6) Graft subclade trees onto the backbone.
        megatree_dir = self.results_dir / "trees" / "orthophyl"
        megatree_dir.mkdir(parents=True, exist_ok=True)
        merged_tree = megatree_dir / f"{taxon_name}_megatree.nwk"
        conflict_report = megatree_dir / f"{taxon_name}_megatree_conflicts.json"

        cmd = ['python', str(self.megatree_grafter),
               '--backbone', str(backbone_tree),
               '--out-tree', str(merged_tree),
               '--out-report', str(conflict_report),
               '--min-support', str(self.conflict_min_support)]
        for sc_name, spec in sorted(subclade_specs.items()):
            if not spec['reps']:
                logger.warning(f"    Subclade {sc_name} contributed no backbone "
                               f"reps; skipping its graft.")
                continue
            cmd.extend(['--subclade',
                        f"{sc_name}:{spec['tree']}:{','.join(spec['reps'])}"])
        self._log_command(cmd)

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would graft {len(subclade_specs)} subclades "
                        f"onto backbone -> {merged_tree}")
        else:
            log_file = self.logs_dir / f"megatree_{taxon_name}.log"
            with open(log_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                                        text=True)
            if result.returncode != 0:
                raise RuntimeError(
                    f"Megatree graft failed for {taxon_name}. Check log: {log_file}")
            logger.info(f"  ✓ Megatree written: {merged_tree}")
            logger.info(f"    Conflict report: {conflict_report}")
            self._write_checkpoint(f"megatree_graft_{taxon_name}")

        # Create the taxon DB from the backbone run (coherent HMM set for ReLeaf).
        self._create_database_entry(
            taxon_name=taxon_name, orthophyl_output=backbone_out,
            taxonomy=taxonomy)

    @staticmethod
    def _locate_species_tree(orthophyl_output: Path) -> Path:
        """Find an OrthoPhyl run's species tree by pattern, tolerating the known
        filename variants.

        script_lib/functions.sh writes FINAL_SPECIES_TREES/iqtree.SCO_strict.CDS.tree
        while other code globs SCO_strict.CDS.iqtree.treefile; we match by pattern
        (mirroring create_hierarchical_database.py) so grafting is not broken by the
        ordering. Prefers IQ-TREE trees, then any .treefile/.tree/.nwk.
        """
        tree_dirs = [
            orthophyl_output / "FINAL_SPECIES_TREES",
            orthophyl_output / "phylo_current" / "SpeciesTree",
        ]
        for tree_dir in tree_dirs:
            if not tree_dir.exists():
                continue
            candidates = (list(tree_dir.glob("*iqtree*.tree*")) +
                          list(tree_dir.glob("*.treefile")) +
                          list(tree_dir.glob("*.tree")) +
                          list(tree_dir.glob("*.nwk")))
            # De-dup preserving order.
            seen = set()
            for c in candidates:
                if c not in seen and c.exists():
                    return c
                seen.add(c)
        # Fall back to the canonical path so the error message is actionable.
        return orthophyl_output / "FINAL_SPECIES_TREES" / "SCO_strict.CDS.iqtree.treefile"

    def _build_subclade(self, taxon_name: str, entry: Dict, raw_dir: Path,
                        query_assemblies: List[Dict], taxonomy: str,
                        is_subclade: bool, force: bool = False):
        """QC one subclade's raw members, run OrthoPhyl, and create its DB entry.

        For the unpartitioned case (is_subclade=False) sc_name == taxon_name and
        the whole raw set is the member list. CheckM2 QC runs HERE (deferred from
        partition time), only on this subclade's members.

        force=True overwrites an existing DB dir for this subclade (e.g. when
        rebuilding), rather than aborting on FileExistsError.
        """
        sc_name = entry['name']
        sc_id = entry.get('subclade_id')
        subclade_dir = self.orthophyl_dir / "downloads" / sc_name

        # Determine this subclade's raw member paths.
        member_names = self._read_members_file(entry.get('members_file'))
        if member_names:
            raw_members = [raw_dir / m for m in member_names]
        else:
            # Unpartitioned: everything in raw_dir.
            raw_members = list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))

        ckey = sc_name  # checkpoint key component (already suffixed for subclades)

        # ---- QC ----
        if self._check_checkpoint(f"qc_{ckey}") and self.resume:
            logger.info(f"  ✓ QC already complete for {sc_name} (resuming)")
            genomes_to_keep = subclade_dir / "genomes_to_keep"
        else:
            genomes_to_keep = self._qc_subclade(
                subclade_dir, raw_members, taxon_label=sc_name,
                query_assemblies=query_assemblies)
            self._write_checkpoint(f"qc_{ckey}")

        # Ensure query genomes present (fallback, e.g. --skip-download / dry-run).
        if not self.dry_run:
            if not genomes_to_keep.exists():
                genomes_to_keep.mkdir(parents=True, exist_ok=True)
            for asm in query_assemblies:
                src = Path(asm['assembly_path'])
                normalized = genomes_to_keep / (src.stem + ".fna")
                dst = genomes_to_keep / src.name
                if normalized.exists() or dst.exists():
                    continue
                shutil.copy(src, normalized)
                logger.info(f"    Added (no QC): {asm['assembly_id']}")

            # Post-QC floor: OrthoPhyl needs >=4 genomes for SCO analysis.
            n_kept = len(list(genomes_to_keep.glob("*.fna")) +
                         list(genomes_to_keep.glob("*.fasta")))
            logger.info(f"  Genomes passing QC for {sc_name}: {n_kept}")
            if n_kept < 4:
                raise RuntimeError(
                    f"Subclade {sc_name} has only {n_kept} genomes after QC "
                    f"(< 4 required for OrthoPhyl). Raw member count was "
                    f"{len(raw_members)}; QC dropped too many. Consider raising "
                    f"--max-tree-genomes or relaxing QC.")

        # ---- OrthoPhyl ----
        orthophyl_output = self.orthophyl_dir / "orthophyl_runs" / sc_name
        if self._check_checkpoint(f"orthophyl_{ckey}") and self.resume:
            logger.info(f"  ✓ OrthoPhyl already complete for {sc_name} (resuming)")
        else:
            self._run_orthophyl(
                input_dir=genomes_to_keep,
                output_dir=orthophyl_output,
                taxon_name=sc_name,
                assemblies=query_assemblies)
            self._write_checkpoint(f"orthophyl_{ckey}")

        # ---- Database entry ----
        subclade_meta = None
        if is_subclade:
            subclade_meta = {
                'is_subclade': True,
                'parent_taxon': taxon_name,
                'subclade_id': sc_id,
                'sketch_file': entry.get('sketch_file'),
                'members_file': entry.get('members_file'),
                'source_genome_dir': str(raw_dir),
                'built': True,
            }
        if self._check_checkpoint(f"database_{ckey}") and self.resume:
            logger.info(f"  ✓ Database already created for {sc_name} (resuming)")
        else:
            self._create_database_entry(
                taxon_name=sc_name,
                orthophyl_output=orthophyl_output,
                taxonomy=taxonomy,
                subclade_meta=subclade_meta,
                force=force)
            self._write_checkpoint(f"database_{ckey}")

    def _download_genomes(self, taxon_name: str, output_dir: Path,
                          query_assemblies: Optional[List[Dict]] = None,
                          download_only: bool = False, qc_only: bool = False):
        """Run gather_filter_asms.sh to download and/or QC-filter genomes.

        query_assemblies, if given, are the user's input assemblies; their FASTA
        paths are passed via --query-genomes so they are QC-filtered alongside the
        downloads instead of bypassing QC.

        Modes (mutually exclusive; default runs the full download+QC pipeline as
        before, byte-identical to prior behaviour):
          download_only : run only the NCBI download + staging steps, stopping
                          before CheckM2 QC. Raw genomes land at
                          <output_dir>/assemblies_all.TMP/*.fna. No genomes_to_keep/
                          or stats table is produced, so those checks are skipped.
          qc_only       : skip the download; assume <output_dir>/assemblies_all.TMP/
                          already holds the raw member FASTAs (the wrapper staged
                          them). Runs only CheckM2 QC + filtering -> genomes_to_keep/.
        """
        if download_only and qc_only:
            raise ValueError("download_only and qc_only are mutually exclusive")

        if not self.gather_script:
            logger.warning(f"  No gather script provided, skipping download for {taxon_name}")
            logger.warning(f"  Please manually download genomes to: {output_dir}/genomes_to_keep/")
            return

        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            str(self.gather_script),
            taxon_name,
            str(output_dir),
            str(self.threads)
        ]

        if download_only:
            cmd.append('--download-only')
        elif qc_only:
            cmd.append('--qc-only')

        # QC-time flags only matter when QC runs (full or qc_only). In
        # download_only mode CheckM2 never runs, so skip them.
        qc_runs = not download_only
        if qc_runs:
            if self.use_bbmap:
                cmd.append('--use-bbmap')
                logger.info(f"  Using bbmap statswrapper instead of CheckM2")
            elif self.low_ram:
                cmd.append('--lowmem')
                logger.info(f"  Using CheckM2 --lowmem option (low RAM mode)")

            # Genome-retention controls. Pass the raw --must-keep value straight
            #   through -- the gather script resolves list-vs-file itself.
            if self.must_keep:
                cmd.extend(['--must-keep', self.must_keep])
                logger.info(f"  Enforcing must-keep genomes: {self.must_keep}")
            if self.keep_failing_query:
                cmd.append('--keep-failing-query')
                logger.info(f"  Query genomes failing QC will be kept with a warning")

        # --query-genomes is needed by BOTH phases: download_only stages the query
        # as a raw leaf (so partitioning sees it); qc_only QCs it. In full mode it
        # does both. Only omit for a qc_only call with no queries in this subclade.
        if query_assemblies:
            query_paths = ','.join(str(Path(a['assembly_path'])) for a in query_assemblies)
            cmd.extend(['--query-genomes', query_paths])
            logger.info(f"  Including {len(query_assemblies)} query genome(s)")

        mode_label = ("download-only" if download_only else
                      "qc-only" if qc_only else "download+QC")
        logger.info(f"  Running gather ({mode_label}) for {taxon_name}...")
        self._log_command(cmd)
        logger.info(f"    Output: {output_dir}")

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would run gather ({mode_label}) for {taxon_name}")
            return

        log_file = self.logs_dir / f"download_{taxon_name}.log"
        with open(log_file, 'w') as f:
            if self.verbose == 1:
                result = subprocess.run(
                    cmd,
                    stderr=f,
                    text=True
                )
            elif self.verbose >= 2:
                result = subprocess.run(
                    cmd,
                    text=True
                )
            else:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    text=True
                )

        if result.returncode != 0:
            raise RuntimeError(f"Genome download failed for {taxon_name}. Check log: {log_file}")

        if download_only:
            # No QC ran: verify only that raw genomes were materialized.
            raw_dir = output_dir / "assemblies_all.TMP"
            raw_files = list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))
            if not raw_files:
                raise RuntimeError(
                    f"Download-only completed but no raw genomes found in {raw_dir}\n"
                    f"This could mean no genomes are available for taxon '{taxon_name}',\n"
                    f"or an incorrect taxon name. Check log: {log_file}"
                )
            logger.info(f"  ✓ Downloaded {len(raw_files)} raw genomes (pre-QC)")
            return

        # QC ran (full or qc_only): validate the stats table and genomes_to_keep/.
        # Defense-in-depth: verify the QC stats table was actually populated.
        # If CheckM2 crashes (e.g. OOM), gather_filter_asms.sh can leave
        # assemblies_all.stats.txt as a header-only file, which causes the stats
        # filter to match nothing and pass EVERY raw assembly through unfiltered.
        # The bash script now aborts in that case, but we double-check here so a
        # regression can never silently feed unfiltered genomes into OrthoPhyl.
        if not self.use_bbmap:
            stats_file = output_dir / "assemblies_all.stats.txt"
            if not stats_file.exists():
                raise RuntimeError(
                    f"QC stats file not found after download: {stats_file}\n"
                    f"CheckM2 may have failed (possibly out of memory).\n"
                    f"Check log: {log_file}"
                )
            with open(stats_file) as sf:
                n_stats_rows = sum(1 for _ in sf) - 1  # subtract header
            if n_stats_rows <= 0:
                raise RuntimeError(
                    f"QC stats file contains no per-assembly rows: {stats_file}\n"
                    f"This usually means CheckM2 crashed (e.g. out of memory) and no\n"
                    f"quality filtering was applied. Refusing to proceed with unfiltered\n"
                    f"assemblies.\n"
                    f"Fix: re-run with --low-ram (CheckM2 --lowmem) or --use-bbmap,\n"
                    f"or allocate more memory.\n"
                    f"Check log: {log_file}"
                )

        # Verify download success
        genomes_to_keep = output_dir / "genomes_to_keep"
        if not genomes_to_keep.exists():
            raise FileNotFoundError(
                f"Download appeared to succeed but expected directory not found: {genomes_to_keep}\n"
                f"Check log: {log_file}"
            )

        # Count downloaded genomes
        genome_files = list(genomes_to_keep.glob("*.fna")) + list(genomes_to_keep.glob("*.fasta"))
        n_genomes = len(genome_files)

        if n_genomes == 0:
            raise RuntimeError(
                f"Download completed but no genomes found in {genomes_to_keep}\n"
                f"This could mean:\n"
                f"  - No genomes available for taxon '{taxon_name}' in NCBI\n"
                f"  - All genomes filtered out due to quality thresholds\n"
                f"  - Incorrect taxon name\n"
                f"Check log: {log_file}"
            )

        logger.info(f"  ✓ Downloaded and filtered {n_genomes} genomes")

        # Write success marker for this download
        success_file = output_dir / ".download_complete"
        with open(success_file, 'w') as f:
            f.write(f"Download completed: {datetime.now().isoformat()}\n")
            f.write(f"Taxon: {taxon_name}\n")
            f.write(f"Genomes: {n_genomes}\n")

    # ------------------------------------------------------------------ #
    # Subclade partitioning helpers (pre-QC MASH partition, used by --megatree)
    # ------------------------------------------------------------------ #

    def _download_raw(self, taxon_name: str, output_dir: Path,
                      query_assemblies: Optional[List[Dict]] = None) -> Path:
        """Download the RAW genome set (no QC) via gather --download-only.

        Returns the directory holding the raw FASTAs (<output_dir>/assemblies_all.TMP).
        Query genomes are staged into it as leaves so partitioning can see them.
        """
        self._download_genomes(taxon_name, output_dir,
                               query_assemblies=query_assemblies, download_only=True)
        return output_dir / "assemblies_all.TMP"

    def _qc_subclade(self, subclade_dir: Path, raw_member_paths: List[Path],
                     taxon_label: str,
                     query_assemblies: Optional[List[Dict]] = None) -> Path:
        """QC one subclade's raw members via gather --qc-only.

        Stages the given raw member FASTAs into <subclade_dir>/assemblies_all.TMP/
        (symlinks, falling back to copies), then runs CheckM2 QC + filtering.
        Returns <subclade_dir>/genomes_to_keep.
        """
        raw_tmp = subclade_dir / "assemblies_all.TMP"
        raw_tmp.mkdir(parents=True, exist_ok=True)

        if not self.dry_run:
            for src in raw_member_paths:
                src = Path(src)
                dst = raw_tmp / src.name
                if dst.exists() or dst.is_symlink():
                    continue
                try:
                    os.symlink(os.path.abspath(src), dst)
                except OSError:
                    shutil.copy(src, dst)

        self._download_genomes(taxon_label, subclade_dir,
                               query_assemblies=query_assemblies, qc_only=True)
        return subclade_dir / "genomes_to_keep"

    def _must_keep_stems(self) -> List[str]:
        """Parse self.must_keep (comma-sep accession list OR file path) to stems.

        Returns [] when unset. Accessions are returned as-is; the subsampler
        matches them against member basenames by stem, so extensions don't matter.
        """
        if not self.must_keep:
            return []
        p = Path(self.must_keep)
        if p.exists() and p.is_file():
            with open(p) as f:
                return [ln.strip() for ln in f if ln.strip()]
        return [s.strip() for s in self.must_keep.split(',') if s.strip()]

    def _stage_local_genomes(self, dest: Path) -> Path:
        """Normalize every FASTA in self.genome_dir into <dest>/<stem>.fna.

        Mirrors gather_filter_asms.sh's stage_query_genomes normalization (the
        QC/OrthoPhyl halves only ever glob "*.fna"), so any mix of
        .fna/.fa/.fasta (optionally .gz) is accepted without silently dropping
        files. Originals in self.genome_dir are never modified: gzip files are
        decompressed into dest, plain files are symlinked (falling back to a
        copy) so OrthoPhyl.sh's in-place contig-name sed never touches the
        user's source directory.

        Raises ValueError if two source files normalize to the same stem
        (e.g. "foo.fna" and "foo.fasta" both present) -- silently picking one
        would drop a genome without any indication to the user.
        """
        dest.mkdir(parents=True, exist_ok=True)
        by_stem: Dict[str, Path] = {}
        for ext in self._GENOME_EXTENSIONS:
            for src in sorted(self.genome_dir.glob(f"*{ext}")):
                stem = src.name[:-len(ext)]
                if stem in by_stem:
                    raise ValueError(
                        f"--genome-dir has two files that normalize to the same "
                        f"stem '{stem}': {by_stem[stem].name} and {src.name}. "
                        f"Rename one so each genome has a unique stem.")
                by_stem[stem] = src

        for stem, src in by_stem.items():
            dst = dest / f"{stem}.fna"
            if dst.exists() or dst.is_symlink():
                continue
            if src.name.endswith('.gz'):
                with gzip.open(src, 'rb') as fin, open(dst, 'wb') as fout:
                    shutil.copyfileobj(fin, fout)
            else:
                try:
                    os.symlink(os.path.abspath(src), dst)
                except OSError:
                    shutil.copy(src, dst)

        logger.info(f"  Staged {len(by_stem)} genome(s) from {self.genome_dir} -> {dest}")
        return dest

    def _subsample_genomes(self, taxon_name: str, raw_dir: Path,
                           target: int,
                           must_keep_stems: Optional[List[str]] = None) -> Path:
        """Diverse-subsample an oversized raw set to `target` genomes; return a dir.

        Runs python_scripts/subsample_genomes.py (MASH greedy max-min -- linear
        sketch, no O(n^2) matrix) on the raw FASTAs, then stages the selected
        members into <orthophyl_dir>/subsample/<taxon>/selected/ (symlinks, falling
        back to copies). Returns that directory so the caller's single-tree
        _build_subclade path can consume it exactly like a raw download dir.

        must_keep_stems seed the greedy pick so query/must-keep genomes are always
        retained. Guarded by a `subsample_<taxon>` checkpoint for --resume.
        """
        sub_dir = self.orthophyl_dir / "subsample" / taxon_name
        sub_dir.mkdir(parents=True, exist_ok=True)
        selected_dir = sub_dir / "selected"
        manifest_path = sub_dir / "subsample_manifest.json"

        # Resume: re-mashing risks a different pick, so reuse the manifest.
        if self.resume and manifest_path.exists() and selected_dir.exists():
            logger.info(f"  ✓ Subsample already complete (resuming): {manifest_path}")
            return selected_dir

        cmd = [
            'python', str(self.subsampler),
            '--genome-dir', str(raw_dir),
            '--out-dir', str(sub_dir),
            '--n', str(target),
            '--threads', str(self.threads),
        ]
        for stem in (must_keep_stems or []):
            cmd.extend(['--must-keep', stem])

        self._log_command(cmd)
        if self.dry_run:
            logger.info(f"  [DRY RUN] Would subsample {taxon_name} to {target} genomes")
            selected_dir.mkdir(parents=True, exist_ok=True)
            return selected_dir

        log_file = self.logs_dir / f"subsample_{taxon_name}.log"
        with open(log_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"Subsampling failed for {taxon_name}. Check log: {log_file}")

        with open(manifest_path) as f:
            manifest = json.load(f)

        # Stage the selected members into selected/ (symlink, copy fallback).
        selected_dir.mkdir(parents=True, exist_ok=True)
        for name in manifest.get('members', []):
            src = raw_dir / name
            dst = selected_dir / name
            if dst.exists() or dst.is_symlink():
                continue
            if not src.exists():
                logger.warning(f"    Subsample member missing from raw dir: {name}")
                continue
            try:
                os.symlink(os.path.abspath(src), dst)
            except OSError:
                shutil.copy(src, dst)

        logger.info(f"  Subsampled {taxon_name}: {manifest.get('n_selected')} of "
                    f"{manifest.get('n_total')} genomes selected (target {target})")
        return selected_dir

    def _enforce_total_genome_ceiling(self, taxon_name: str, raw_count: int):
        """Refuse to partition an over-large raw set (guardrail).

        The partitioner's dense NxN MASH-distance matrix is O(n^2) in time and
        memory; past a few thousand genomes it OOM-kills the node. Rather than
        fail opaquely mid-`mash triangle`, stop here with actionable guidance.

        This is a hard stop until one of the large-taxon strategies (--subsample
        to build one tree from a diverse subset, or --megatree to build per-
        subclade trees and merge) is selected. Both will route around this check
        with their own bounded handling.
        """
        if raw_count <= self.max_total_genomes:
            return
        raise RuntimeError(
            f"Taxon '{taxon_name}' has {raw_count} raw genomes, exceeding the "
            f"--max-total-genomes ceiling of {self.max_total_genomes}.\n"
            f"Partitioning builds an all-vs-all MASH distance matrix that grows "
            f"as O(n^2) in memory (~{(raw_count ** 2 * 8) / 1e9:.1f} GB at this "
            f"size) and would likely exhaust RAM.\n"
            f"Options:\n"
            f"  - Choose a more specific taxon/rank so fewer genomes are pulled.\n"
            f"  - Raise --max-total-genomes if you have the memory for an "
            f"{raw_count}x{raw_count} matrix.\n"
            f"  - Use a large-taxon strategy (--subsample / --megatree) once "
            f"available.")

    def _partition_genomes(self, taxon_name: str, raw_genome_dir: Path,
                           query_assemblies: List[Dict],
                           max_size: Optional[int] = None) -> Dict:
        """Run subclade_partition.py on the RAW genome set; return the manifest dict.

        Builds <orthophyl_dir>/partitions/<taxon>/ holding MASH_out, per-subclade
        .msh/.members.txt, and partition_manifest.json. Every query stem must land
        in exactly one subclade (it is a clustering leaf), which we assert.

        max_size is the per-subclade genome ceiling (`--max-size`); it defaults to
        self.max_tree_genomes. The megatree path passes self.subclade_size instead.
        """
        if max_size is None:
            max_size = self.max_tree_genomes
        part_dir = self.orthophyl_dir / "partitions" / taxon_name
        part_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = part_dir / "partition_manifest.json"

        # Resume: never re-mash (would risk renumbering) -- re-read the manifest.
        if self.resume and manifest_path.exists():
            logger.info(f"  ✓ Partition manifest exists (resuming): {manifest_path}")
            with open(manifest_path) as f:
                return json.load(f)

        cmd = [
            'python', str(self.subclade_partitioner),
            '--genome-dir', str(raw_genome_dir),
            '--taxon', taxon_name,
            '--out-dir', str(part_dir),
            '--max-size', str(max_size),
            '--threads', str(self.threads),
        ]
        for asm in query_assemblies:
            cmd.extend(['--query', Path(asm['assembly_path']).name])

        self._log_command(cmd)
        if self.dry_run:
            logger.info(f"  [DRY RUN] Would partition {taxon_name} into subclades")
            # Synthesize a trivial single-subclade manifest for dry-run flow.
            return {
                'partitioned': False, 'parent_taxon': taxon_name,
                'max_size': max_size, 'n_subclades': 1,
                'subclades': [{'subclade_id': 1, 'name': taxon_name, 'n_genomes': 0,
                               'members_file': None, 'sketch_file': None}],
                'query_assignments': {Path(a['assembly_path']).name: taxon_name
                                      for a in query_assemblies},
            }

        log_file = self.logs_dir / f"partition_{taxon_name}.log"
        with open(log_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Partitioning failed for {taxon_name}. Check log: {log_file}")

        with open(manifest_path) as f:
            manifest = json.load(f)

        # Invariant: every query stem lands in exactly one subclade.
        assignments = manifest.get('query_assignments', {})
        for asm in query_assemblies:
            name = Path(asm['assembly_path']).name
            if name not in assignments:
                raise RuntimeError(
                    f"Query {name} not assigned to any subclade in {manifest_path}; "
                    f"partitioning invariant violated.")
        return manifest

    @staticmethod
    def _read_members_file(members_file: Optional[str]) -> List[str]:
        """Read a subclade .members.txt (one basename per line)."""
        if not members_file:
            return []
        p = Path(members_file)
        if not p.exists():
            return []
        with open(p) as f:
            return [ln.strip() for ln in f if ln.strip()]
    
    def _run_orthophyl(
        self,
        input_dir: Path,
        output_dir: Path,
        taxon_name: str,
        assemblies: List[Dict]
    ):
        """Run OrthoPhyl on combined genome set."""
        cmd = [
            str(self.orthophyl_script),
            '-g', str(input_dir),
            '-s', str(output_dir),
            '-t', str(self.threads),
            '-p', 'iqtree',
            '-o', 'CDS'
        ]
        
        logger.info(f"\n  Running OrthoPhyl for {taxon_name}...")
        self._log_command(cmd)
        logger.info(f"    Input: {input_dir}")
        logger.info(f"    Output: {output_dir}")
        logger.info(f"    INCLUDES {len(assemblies)} query genomes!")

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would run OrthoPhyl for {taxon_name}")
            return
        
        log_file = self.logs_dir / f"orthophyl_{taxon_name}.log"
        with open(log_file, 'w') as f:
            if self.verbose == 1:
                result = subprocess.run(
                    cmd,
                    stderr=f,
                    text=True
                )
            elif self.verbose >= 2:
                result = subprocess.run(
                    cmd,
                    text=True
                )
            else:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    text=True
                )
        
        if result.returncode != 0:
            error_msg = f"OrthoPhyl failed for {taxon_name}. Check log: {log_file}"
            self.pipeline_status['failures'].append({
                'type': 'orthophyl',
                'taxon': taxon_name,
                'error': error_msg,
                'log': str(log_file)
            })
            raise RuntimeError(error_msg)
        
        logger.info(f"  ✓ OrthoPhyl complete for {taxon_name}")
        
        # Verify queries are in tree
        tree_file = output_dir / "FINAL_SPECIES_TREES" / "SCO_strict.CDS.iqtree.treefile"
        if tree_file.exists():
            self._verify_queries_in_tree(tree_file, assemblies)
        else:
            logger.warning(f"  ⚠ Tree file not found: {tree_file}")
        
        # Track success
        self.pipeline_status['successes'].append({
            'type': 'orthophyl',
            'taxon': taxon_name,
            'assemblies': len(assemblies)
        })
    
    def _verify_queries_in_tree(self, tree_file: Path, assemblies: List[Dict]):
        """Verify that query assemblies appear in the tree."""
        with open(tree_file, 'r') as f:
            tree_content = f.read()
        
        all_found = True
        for asm in assemblies:
            if asm['assembly_id'] in tree_content:
                logger.info(f"    ✓ Query {asm['assembly_id']} found in tree")
            else:
                logger.warning(f"    ⚠ Query {asm['assembly_id']} NOT found in tree")
                all_found = False
        
        if all_found:
            logger.info(f"  ✓ All {len(assemblies)} queries verified in tree")
        else:
            logger.warning(f"  ⚠ Some queries missing from tree")
    
    def _create_database_entry(
        self,
        taxon_name: str,
        orthophyl_output: Path,
        taxonomy: str,
        subclade_meta: Optional[Dict] = None,
        force: bool = False,
        taxonomy_source: str = 'ncbi',
        qc_applied: bool = True,
    ):
        """Create new database entry from OrthoPhyl run.

        When subclade_meta is given, the entry is written as a built subclade
        (is_subclade + parent_taxon + sketch/members recorded) via the DB creator's
        --single-clade path; otherwise the classic TSV --update path is used.

        force=True passes --force so an existing DB dir is overwritten (e.g.
        rebuilding a subclade); without it the DB creator refuses (FileExistsError).

        taxonomy_source/qc_applied are provenance flags (default 'ncbi'/True keep
        every existing caller byte-identical); the local-genome-ingest mode passes
        'user_supplied' and the actual QC status.
        """
        logger.info(f"\n  Creating database entry for {taxon_name}...")

        if subclade_meta is not None:
            # Built subclade: register a single clade carrying subclade metadata.
            cmd = [
                'python', str(self.database_creator),
                '--single-clade', taxon_name, taxonomy, str(orthophyl_output),
                '--output-dir', str(self.database_dir),
                '--is-subclade',
                '--parent-taxon', str(subclade_meta.get('parent_taxon')),
                '--taxonomy-source', taxonomy_source,
            ]
            if not qc_applied:
                cmd.append('--qc-not-applied')
            if force:
                cmd.append('--force')
            if subclade_meta.get('subclade_id') is not None:
                cmd.extend(['--subclade-id', str(subclade_meta['subclade_id'])])
            if subclade_meta.get('sketch_file'):
                cmd.extend(['--sketch-file', str(subclade_meta['sketch_file'])])
            if subclade_meta.get('members_file'):
                cmd.extend(['--members-file', str(subclade_meta['members_file'])])
            if subclade_meta.get('source_genome_dir'):
                cmd.extend(['--source-genome-dir', str(subclade_meta['source_genome_dir'])])

            self._log_command(cmd)
            if self.dry_run:
                logger.info(f"  [DRY RUN] Would create subclade database for {taxon_name}")
                return

            log_file = self.logs_dir / f"database_{taxon_name}.log"
            with open(log_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                raise RuntimeError(
                    f"Database creation failed for {taxon_name}. Check log: {log_file}")
            logger.info(f"  ✓ Database created: {taxon_name}_db (subclade of "
                        f"{subclade_meta.get('parent_taxon')})")
            return

        # Update orthophyl_runs.tsv
        tsv_file = self.database_dir / "orthophyl_runs.tsv"

        # Append new entry
        with open(tsv_file, 'a') as f:
            f.write(f"{taxon_name}\t{orthophyl_output}\t{taxonomy}\n")

        if self.verbose:
            logger.info(f"    Added to orthophyl_runs.tsv: {taxon_name}")

        # Run database creator in update mode
        cmd = [
            'python', str(self.database_creator),
            '--input', str(tsv_file),
            '--output-dir', str(self.database_dir),
            '--update',
            '--taxonomy-source', taxonomy_source,
        ]
        if not qc_applied:
            cmd.append('--qc-not-applied')

        self._log_command(cmd)

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would create database for {taxon_name}")
            return
        
        log_file = self.logs_dir / f"database_{taxon_name}.log"
        with open(log_file, 'w') as f:
            if self.verbose == 1:
                result = subprocess.run(cmd, stderr=f, text=True)
            elif self.verbose >= 2:
                result = subprocess.run(cmd, text=True)
            else:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Database creation failed for {taxon_name}. Check log: {log_file}")
        
        logger.info(f"  ✓ Database created: {taxon_name}_db")
    
    def _phase_aggregation(self):
        """Phase 4: Aggregate all results."""
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 4: RESULTS AGGREGATION")
        logger.info("=" * 70)
        
        # Create results directories
        trees_dir = self.results_dir / "trees"
        (trees_dir / "releaf").mkdir(parents=True, exist_ok=True)
        (trees_dir / "orthophyl").mkdir(parents=True, exist_ok=True)
        
        # Collect ReLeaf trees - look in database's ReLeaf_dir (where ReLeaf actually writes)
        releaf_trees = []
        for db_dir in self.database_dir.glob("*_db"):
            releaf_output = db_dir / 'orthophyl_run' / 'ReLeaf_dir'
            # Try multiple possible tree file names
            tree_candidates = [
                releaf_output / "new_trees" / "SCO_strict.CDS.iqtree.treefile.addasm",
                releaf_output / "new_trees" / "SCO_strict.CDS.iqtree.treefile",
                releaf_output / "new_trees" / "SCO_strict.PROT.iqtree.treefile.addasm",
                releaf_output / "new_trees" / "SCO_strict.PROT.iqtree.treefile"
            ]
            for tree_file in tree_candidates:
                if tree_file.exists():
                    dst = trees_dir / "releaf" / f"{db_dir.name}_phylogeny.nwk"
                    shutil.copy(tree_file, dst)
                    releaf_trees.append(db_dir.name)
                    logger.info(f"  ✓ Collected ReLeaf tree: {db_dir.name} ({tree_file.name})")
                    break
        
        # Collect OrthoPhyl trees
        orthophyl_trees = []
        for taxon_dir in (self.orthophyl_dir / "orthophyl_runs").glob("*/"):
            tree_file = taxon_dir / "FINAL_SPECIES_TREES" / "SCO_strict.CDS.iqtree.treefile"
            if tree_file.exists():
                dst = trees_dir / "orthophyl" / f"{taxon_dir.name}_phylogeny.nwk"
                shutil.copy(tree_file, dst)
                orthophyl_trees.append(taxon_dir.name)
                logger.info(f"  ✓ Collected OrthoPhyl tree: {taxon_dir.name}")
        
        # Generate summary report
        self._generate_summary_report(releaf_trees, orthophyl_trees)
        
        logger.info(f"\n✓ Results aggregated in: {self.results_dir}")
        
        self.pipeline_status['phases']['aggregation'] = {
            'status': 'complete',
            'releaf_trees': len(releaf_trees),
            'orthophyl_trees': len(orthophyl_trees)
        }
    
    def _generate_summary_report(self, releaf_trees: List[str], orthophyl_trees: List[str]):
        """Generate human-readable summary report with actual success/failure counts."""
        report_file = self.results_dir / "pipeline_summary.txt"
        
        successes = self.pipeline_status.get('successes', [])
        failures = self.pipeline_status.get('failures', [])
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("ORTHOPHYL PIPELINE - SUMMARY REPORT\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Input File: {self.input_file}\n")
            f.write(f"Database Directory: {self.database_dir}\n")
            f.write(f"Output Directory: {self.output_dir}\n\n")
            
            f.write("RESULTS:\n")
            f.write("-" * 70 + "\n\n")
            
            # Success counts
            f.write(f"✓ Successful placements: {len(successes)}\n")
            if successes:
                for s in successes:
                    if s['type'] == 'releaf':
                        f.write(f"  - ReLeaf: {s['database']} ({s.get('assemblies', '?')} assemblies)\n")
                    else:
                        f.write(f"  - OrthoPhyl: {s['taxon']} ({s.get('assemblies', '?')} assemblies)\n")
                f.write("\n")
            
            # Failure counts
            if failures:
                f.write(f"✗ Failed placements: {len(failures)}\n")
                for fail in failures:
                    if fail['type'] == 'releaf':
                        f.write(f"  - ReLeaf: {fail['database']}\n")
                    else:
                        f.write(f"  - OrthoPhyl: {fail['taxon']}\n")
                    f.write(f"    Error: {fail['error']}\n")
                    f.write(f"    Log: {fail['log']}\n")
                f.write("\n")
            
            f.write("OUTPUT LOCATIONS:\n")
            f.write(f"  Trees: {self.results_dir}/trees/\n")
            f.write(f"  Logs: {self.logs_dir}/\n")
            f.write(f"  Routing decisions: {self.routing_dir}/\n")
            f.write("\n")
            
            # Only show success banner if no failures
            f.write("=" * 70 + "\n")
            if not failures:
                f.write("All assemblies have been placed in phylogenetic trees!\n")
            else:
                f.write(f"PIPELINE COMPLETED WITH {len(failures)} FAILURE(S)\n")
                f.write("Please review the errors above and check the log files.\n")
            f.write("=" * 70 + "\n")
        
        logger.info(f"  Summary report: {report_file}")
        
        # Return failure count for exit code
        return len(failures)
    
    def _check_checkpoint(self, name: str) -> bool:
        """Check if a checkpoint exists."""
        checkpoint_file = self.checkpoint_dir / f"{name}.flag"
        return checkpoint_file.exists()
    
    def _write_checkpoint(self, name: str):
        """Write a checkpoint flag."""
        checkpoint_file = self.checkpoint_dir / f"{name}.flag"
        checkpoint_file.write_text(datetime.now().isoformat())
    
    def _verify_download_complete(self, download_dir: Path, taxon_name: str) -> bool:
        """Verify that genome download completed successfully.
        
        Checks for:
        1. Checkpoint flag exists
        2. Download success marker exists
        3. genomes_to_keep directory exists and contains genomes
        
        Returns:
            bool: True if download is complete and valid
        """
        # Check checkpoint
        if not self._check_checkpoint(f"download_{taxon_name}"):
            return False
        
        # Check success marker
        success_file = download_dir / ".download_complete"
        if not success_file.exists():
            logger.warning(f"  ⚠ Checkpoint exists but no success marker found")
            logger.warning(f"    Download may have been interrupted")
            return False
        
        # Check genomes_to_keep directory
        genomes_to_keep = download_dir / "genomes_to_keep"
        if not genomes_to_keep.exists():
            logger.warning(f"  ⚠ Success marker exists but genomes_to_keep directory not found")
            return False
        
        # Check for genome files
        genome_files = list(genomes_to_keep.glob("*.fna")) + list(genomes_to_keep.glob("*.fasta"))
        if len(genome_files) == 0:
            logger.warning(f"  ⚠ genomes_to_keep directory exists but contains no genomes")
            return False
        
        # All checks passed
        return True
    
    def _load_routing_results(self) -> Dict:
        """Load previously computed routing results."""
        return self._parse_routing_results()
    
    def _create_mock_routing_results(self) -> Dict:
        """Create mock routing results for dry run preview."""
        logger.info("\n  [DRY RUN] Parsing input file to preview routing...")
        
        # This is a simplified preview - in real mode, assembly_router does this
        releaf_batch = []
        orthophyl_batch = defaultdict(list)
        
        # Parse input file
        with open(self.input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue
                
                assembly_path = fields[0]
                taxonomy = fields[1]
                assembly_id = fields[2] if len(fields) > 2 else Path(assembly_path).stem
                
                # Simple heuristic: if taxonomy is very specific, might match a database
                # In reality, assembly_router checks actual databases
                logger.info(f"    Preview: {assembly_id} - {taxonomy[:50]}...")
        
        logger.info("\n  [DRY RUN] In real mode, assembly_router would determine:")
        logger.info("    - Which assemblies match existing databases (→ ReLeaf)")
        logger.info("    - Which assemblies need new databases (→ OrthoPhyl)")
        logger.info("    Run without --dry-run to see actual routing decisions\n")

        return {'releaf_batch': releaf_batch, 'orthophyl_batch': dict(orthophyl_batch)}
    
    def _check_existing_taxon_database(self) -> Optional[Dict]:
        """Check if a database exists for the specified taxon.
        
        Returns:
            Dict with database info if found, None otherwise
        """
        logger.info(f"\nChecking for existing database for taxon: {self.taxon}")
        
        if not self.database_dir or not self.database_dir.exists():
            return None
        
        # Search all databases for matching taxon
        for db_dir in self.database_dir.glob("*_db"):
            config_file = db_dir / "database_config.json"
            if not config_file.exists():
                continue
            
            try:
                with open(config_file, 'r') as f:
                    config = json.load(f)
                
                # Check if source_taxon_name matches
                if config.get('source_taxon_name') == self.taxon:
                    logger.info(f"  ✓ Found matching database by taxon name")
                    return {
                        'db_dir': db_dir,
                        'config': config,
                        'n_assemblies': len(config.get('assembly_accessions', [])),
                        'clade_name': config.get('clade_name')
                    }
                
                # Also check clade_name for fuzzy match
                if config.get('clade_name', '').lower() == self.taxon.lower():
                    logger.info(f"  ✓ Found matching database by clade name")
                    return {
                        'db_dir': db_dir,
                        'config': config,
                        'n_assemblies': config.get('n_genomes', 0),
                        'clade_name': config.get('clade_name')
                    }
            except Exception as e:
                logger.warning(f"  ⚠ Error reading {config_file}: {e}")
                continue
        
        return None
    
    def _run_taxon_create_mode(self) -> int:
        """Create new database from taxon query."""
        logger.info("\n" + "=" * 70)
        logger.info("TAXON MODE: CREATE NEW DATABASE")
        logger.info("=" * 70)
        
        # Import taxon gatherer
        sys.path.insert(0, str(self.script_dir / "utils"))
        try:
            from taxon_assembly_gatherer import TaxonAssemblyGatherer
        except ImportError as e:
            raise ImportError(f"Failed to import TaxonAssemblyGatherer: {e}")
        
        # Resolve taxon identity (taxid/rank/lineage) for database metadata.
        # This constructor does NOT download the assembly_summary tables; those
        # are only fetched by query_ncbi(), which create mode no longer calls.
        logger.info(f"\nResolving taxonomy for {self.taxon}...")
        gatherer = TaxonAssemblyGatherer(
            taxon=self.taxon,
            rank=self.taxon_rank,
            output_dir=self.output_dir / "taxon_query"
        )

        if self.dry_run:
            logger.info("  [DRY RUN] Would download and QC-filter assemblies")
            logger.info("=" * 70)
            logger.info("PIPELINE COMPLETE (DRY RUN)!")
            logger.info("=" * 70)
            self._save_final_status()
            return 0
        
        # NOTE: We intentionally do NOT call gatherer.query_ncbi() here. In create
        # mode the genomes are downloaded and QC-filtered entirely by
        # gather_filter_asms.sh (below); query_ncbi() would re-download the full
        # assembly_summary tables just to build a metadata accession list. We
        # instead derive that list from the QC-filtered genomes_to_keep/ output,
        # which is both cheaper and more accurate (post-QC, not pre-QC candidates).
        # The gatherer object is still used for taxonomy metadata (taxid, rank,
        # lineage) resolved in its constructor. Empty-taxon detection is handled
        # by the genomes_to_keep QC check below.

        # Verify gather script is available (required for create mode)
        if not self.gather_script or not self.gather_script.exists():
            logger.error(f"\n❌ ERROR: Genome download script is required for taxon create mode")
            logger.error(f"  Please provide --gather-script utils/gather_filter_asms.sh")
            return 1
        
        # Download the RAW assemblies (no QC yet) so we can MASH-partition before
        # spending CheckM2 compute -- same pre-QC ordering as batch mode.
        logger.info(f"\nDownloading {self.taxon} assemblies (raw, pre-QC)...")
        logger.info(f"  Using: {self.gather_script}")
        download_dir = self.output_dir / "downloaded_assemblies"
        raw_dir = download_dir / "assemblies_all.TMP"

        if (self.skip_download or
                (self._check_checkpoint(f"download_{self.taxon}") and self.resume
                 and raw_dir.exists())):
            logger.info(f"  ✓ Raw download already complete, skipping")
        else:
            raw_dir = self._download_raw(self.taxon, download_dir)
            self._write_checkpoint(f"download_{self.taxon}")

        raw_files = (list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))
                     if raw_dir.exists() else [])
        raw_count = len(raw_files)
        if raw_count == 0:
            logger.error(f"\n❌ ERROR: No genomes downloaded for {self.taxon}")
            logger.error(f"  Check {download_dir} for logs")
            return 1
        logger.info(f"  ✓ {raw_count} raw genomes downloaded")

        taxonomy = gatherer.get_taxonomy_string()

        # Opt-in megatree: full-coverage partition -> per-subclade trees -> graft.
        # Create mode has no query, so every subclade is built. Under the ceiling
        # this falls through to the normal single-tree path.
        if self.megatree and raw_count > self.max_tree_genomes:
            self._run_megatree(
                taxon_name=self.taxon, raw_dir=raw_dir,
                query_assemblies=[], taxonomy=taxonomy)
            logger.info("=" * 70)
            logger.info("TAXON MODE COMPLETE (megatree)!")
            logger.info("=" * 70)
            self._save_final_status()
            return 0

        # Over the ceiling, the DEFAULT is to diverse-subsample to one tree (no
        # O(n^2) matrix). Create mode has no query; must-keep accessions (if any)
        # seed the pick so they are retained. The subsampled dir replaces raw_dir/
        # raw_files, then flows through the single-tree QC+build path below.
        if raw_count > self.max_tree_genomes:
            logger.info(f"  Raw count {raw_count} > max_tree_genomes "
                        f"{self.max_tree_genomes}: diverse-subsampling to "
                        f"{self.subsample_size} genomes")
            raw_dir = self._subsample_genomes(
                self.taxon, raw_dir, self.subsample_size,
                must_keep_stems=self._must_keep_stems())
            raw_files = (list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))
                         if raw_dir.exists() else [])
            raw_count = len(raw_files)

        # Single tree: QC the whole (possibly subsampled) raw set, build, and
        # create the taxon-flavored DB.
        genomes_to_keep = self._qc_subclade(
            download_dir, raw_files, taxon_label=self.taxon, query_assemblies=[])
        kept = (list(genomes_to_keep.glob("*.fna")) +
                list(genomes_to_keep.glob("*.fasta")))
        if len(kept) < 4:
            logger.error(f"\n❌ ERROR: only {len(kept)} genomes passed QC (<4)")
            return 1
        orthophyl_output = self.output_dir / "orthophyl_run"
        self._run_orthophyl(
            input_dir=genomes_to_keep, output_dir=orthophyl_output,
            taxon_name=self.taxon, assemblies=[])
        logger.info(f"\nCreating database for {self.taxon}...")
        self._create_taxon_database(
            taxon_name=self.taxon, orthophyl_output=orthophyl_output,
            gatherer=gatherer, genomes_to_keep=genomes_to_keep)
        logger.info("=" * 70)
        logger.info("TAXON MODE COMPLETE!")
        logger.info("=" * 70)
        logger.info(f"  Database created: {self.taxon}_db")
        logger.info(f"  Assemblies: {len(kept)}")
        self._save_final_status()
        return 0
    
    def _run_taxon_update_mode(self, existing_db: Dict) -> int:
        """Update existing database with new assemblies."""
        logger.info("\n" + "=" * 70)
        logger.info("TAXON MODE: UPDATE EXISTING DATABASE")
        logger.info("=" * 70)
        
        # Import taxon gatherer
        sys.path.insert(0, str(self.script_dir / "utils"))
        try:
            from taxon_assembly_gatherer import TaxonAssemblyGatherer
        except ImportError as e:
            raise ImportError(f"Failed to import TaxonAssemblyGatherer: {e}")
        
        # Query NCBI for assemblies
        logger.info(f"\nQuerying NCBI for {self.taxon} assemblies...")
        gatherer = TaxonAssemblyGatherer(
            taxon=self.taxon,
            rank=self.taxon_rank,
            output_dir=self.output_dir / "taxon_query"
        )

        if self.dry_run:
            logger.info("  [DRY RUN] Would check for new assemblies and update database")
            logger.info("=" * 70)
            logger.info("PIPELINE COMPLETE (DRY RUN)!")
            logger.info("=" * 70)
            self._save_final_status()
            return 0
        
        # Get all assemblies
        all_assemblies = gatherer.query_ncbi()
        
        if not all_assemblies:
            logger.error(f"✗ No assemblies found for taxon '{self.taxon}'")
            return 1
        
        logger.info(f"  ✓ Found {len(all_assemblies)} total assemblies in NCBI")
        
        # Compare with existing database
        existing_accessions = set(existing_db['config'].get('assembly_accessions', []))
        new_assemblies = [a for a in all_assemblies if a['accession'] not in existing_accessions]
        
        if not new_assemblies:
            logger.info(f"\n✓ Database is up to date! No new assemblies found.")
            logger.info(f"  Current assemblies: {len(existing_accessions)}")
            logger.info(f"  NCBI assemblies: {len(all_assemblies)}")
            self._save_final_status()
            return 0
        
        logger.info(f"\n→ Found {len(new_assemblies)} new assemblies to add")
        logger.info(f"  Existing: {len(existing_accessions)}")
        logger.info(f"  New: {len(new_assemblies)}")
        
        # Download new assemblies (raw, unfiltered) into assemblies_all.TMP so the
        # shared QC helper can stage + filter them exactly like create mode does.
        logger.info(f"\nDownloading {len(new_assemblies)} new assemblies...")
        download_dir = self.output_dir / "new_assemblies"
        raw_dir = download_dir / "assemblies_all.TMP"
        raw_dir.mkdir(parents=True, exist_ok=True)

        gatherer.download_assemblies(new_assemblies, raw_dir)

        raw_files = list(raw_dir.glob("*.fna")) + list(raw_dir.glob("*.fasta"))
        if not raw_files:
            raise RuntimeError(
                f"Downloaded {len(new_assemblies)} new assemblies but no FASTAs "
                f"landed in {raw_dir}. Check the gatherer's download step.")

        # QC-filter the new downloads (CheckM2 completeness/contamination/N50),
        # consistent with the create-mode filtered path. These are NCBI reference
        # genomes being added to grow the DB -- not user queries -- so dropping
        # low-quality ones is the intended behaviour. Falls back to unfiltered when
        # no gather script is configured (CheckM2 unavailable).
        if self.gather_script:
            logger.info(f"\nQC-filtering {len(raw_files)} new assemblies (CheckM2)...")
            genomes_to_keep = self._qc_subclade(
                download_dir, raw_files, taxon_label=self.taxon)
            kept = (list(genomes_to_keep.glob("*.fna")) +
                    list(genomes_to_keep.glob("*.fasta")))
            n_kept = len(kept)
            logger.info(f"  New assemblies passing QC: {n_kept} / {len(raw_files)}")
            if n_kept == 0:
                logger.warning(
                    f"\n⚠ All {len(raw_files)} new assemblies were dropped by QC; "
                    f"nothing to add. Database left unchanged.")
                self._save_final_status()
                return 0
            releaf_input = genomes_to_keep
        else:
            logger.warning(
                "  ⚠ No gather script configured: skipping CheckM2 QC on new "
                "assemblies (adding unfiltered).")
            releaf_input = raw_dir

        # Run ReLeaf to add to existing database
        logger.info(f"\nRunning ReLeaf to add new assemblies to database...")

        db_dir = existing_db['db_dir']
        config = existing_db['config']

        # Get tree method and data type from config
        tree_methods = config.get('available_tree_methods', ['iqtree'])
        tree_data_types = config.get('available_data_types', ['CDS'])
        tree_method = tree_methods[0] if tree_methods else 'iqtree'
        tree_data = tree_data_types[0] if tree_data_types else 'CDS'

        releaf_output = self.output_dir / "releaf_update"

        self._run_releaf(
            database_dir=db_dir,
            input_genomes=releaf_input,
            output_dir=releaf_output,
            tree_method=tree_method,
            tree_data=tree_data,
            database_name=existing_db['clade_name']
        )
        
        # Update database metadata. Record only the accessions that actually
        # survived QC and entered ReLeaf (filenames in releaf_input are
        # accessions), NOT the full pre-QC candidate list -- otherwise dropped
        # genomes would be marked present and never re-tried on a later update.
        added_accessions = sorted(
            p.stem for p in
            list(releaf_input.glob("*.fna")) + list(releaf_input.glob("*.fasta")))
        logger.info(f"\nUpdating database metadata...")
        self._update_taxon_database_metadata(
            db_dir=db_dir,
            added_accessions=added_accessions
        )

        logger.info("=" * 70)
        logger.info("TAXON UPDATE COMPLETE!")
        logger.info("=" * 70)
        logger.info(f"  Database: {existing_db['clade_name']}")
        logger.info(f"  Added assemblies: {len(added_accessions)}")
        logger.info(f"  Total assemblies: {len(existing_accessions) + len(added_accessions)}")

        self._save_final_status()
        return 0

    def _run_local_genomes_mode(self) -> int:
        """Build a tree + database from genomes already on disk (--genome-dir).

        Order is cheap-checks-first so nothing expensive precedes a knowable
        failure: DB-collision pre-flight, taxonomy resolution/validation,
        pre-QC genome count, then (skippable) QC, subsampling if oversized,
        OrthoPhyl, and database creation.
        """
        logger.info("\n" + "=" * 70)
        logger.info("LOCAL GENOME-INGEST MODE")
        logger.info("=" * 70)

        safe = self._default_run_name(self.clade_name)

        # ---- DB-collision pre-flight (before any compute) ----
        sys.path.insert(0, str(self.script_dir / "assembly_router"))
        try:
            from create_hierarchical_database import database_exists
        except ImportError as e:
            raise ImportError(f"Failed to import database_exists: {e}")

        existing = database_exists(self.clade_name, self.database_dir)
        if existing:
            logger.error(
                f"\n✗ Database already exists for clade '{self.clade_name}': {existing}")
            logger.error("  Use a different --clade-name, or remove the existing "
                          "database directory to rebuild.")
            return 1

        # ---- Taxonomy resolution ----
        taxonomy, is_routable = self._resolve_local_taxonomy()
        if not is_routable:
            logger.warning(f"  ⚠ Database '{self.clade_name}' will be created but is "
                            f"NOT NCBI-assigned and will not be matched by fully "
                            f"specified taxonomy queries.")

        # ---- Pre-QC genome count (fail in seconds, not after QC/OrthoPhyl) ----
        raw_candidates = []
        for ext in self._GENOME_EXTENSIONS:
            raw_candidates.extend(self.genome_dir.glob(f"*{ext}"))
        if len(raw_candidates) < 4:
            logger.error(
                f"\n❌ ERROR: --genome-dir has only {len(raw_candidates)} genome "
                f"file(s) (< 4 required for OrthoPhyl): {self.genome_dir}")
            return 1
        logger.info(f"  Found {len(raw_candidates)} genome file(s) in {self.genome_dir}")

        # ---- QC requires the gather script, unless explicitly skipped ----
        if not self.skip_qc and (not self.gather_script or not self.gather_script.exists()):
            logger.error(f"\n❌ ERROR: QC is enabled by default but no genome "
                          f"download/QC script is configured")
            logger.error(f"  Please provide --gather-script utils/gather_filter_asms.sh, "
                         f"or pass --skip-qc to skip QC")
            return 1

        local_dir = self.orthophyl_dir / "local_input" / safe

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would stage genomes from {self.genome_dir}")
            logger.info(f"  [DRY RUN] Would QC: {not self.skip_qc}")
            logger.info(f"  [DRY RUN] Would build tree and create database "
                        f"'{self.clade_name}_db' (taxonomy_source=user_supplied, "
                        f"taxonomy={taxonomy!r})")
            self._save_final_status()
            return 0

        # ---- Stage: normalize every input to <stem>.fna, originals untouched ----
        if self._check_checkpoint(f"stage_local_{safe}") and self.resume:
            logger.info(f"  ✓ Staging already complete for {safe} (resuming)")
            staged_dir = local_dir / "assemblies_all.TMP"
        else:
            staged_dir = self._stage_local_genomes(local_dir / "assemblies_all.TMP")
            self._write_checkpoint(f"stage_local_{safe}")

        staged_files = (list(staged_dir.glob("*.fna")) + list(staged_dir.glob("*.fasta")))

        # ---- Oversized: diverse-subsample before QC ----
        qc_source_dir = local_dir
        if len(staged_files) > self.max_tree_genomes:
            logger.info(f"  Genome count {len(staged_files)} > max_tree_genomes "
                        f"{self.max_tree_genomes}: diverse-subsampling to "
                        f"{self.subsample_size} genomes")
            staged_dir = self._subsample_genomes(
                safe, staged_dir, self.subsample_size,
                must_keep_stems=self._must_keep_stems())
            staged_files = (list(staged_dir.glob("*.fna")) +
                             list(staged_dir.glob("*.fasta")))
            # Use a distinct qc/ dir -- reusing local_dir would re-stage the FULL
            # (pre-subsample) set into assemblies_all.TMP and QC everything.
            qc_source_dir = local_dir / "qc"

        # ---- QC (default on; --skip-qc bypasses it entirely) ----
        qc_applied = not self.skip_qc
        if self.skip_qc:
            logger.warning("  ⚠ --skip-qc: genomes will NOT be quality-checked "
                            "(no CheckM2 completeness/contamination/N50 filtering)")
            genomes_dir = staged_dir
        else:
            if self._check_checkpoint(f"qc_{safe}") and self.resume:
                logger.info(f"  ✓ QC already complete for {safe} (resuming)")
                genomes_dir = qc_source_dir / "genomes_to_keep"
            else:
                genomes_dir = self._qc_subclade(
                    qc_source_dir, staged_files, taxon_label=safe)
                self._write_checkpoint(f"qc_{safe}")

        kept = list(genomes_dir.glob("*.fna")) + list(genomes_dir.glob("*.fasta"))
        if len(kept) < 4:
            logger.error(f"\n❌ ERROR: only {len(kept)} genomes remain "
                         f"(< 4 required for OrthoPhyl)")
            self._save_final_status()
            return 1
        logger.info(f"  Genomes proceeding to OrthoPhyl: {len(kept)}")

        # ---- OrthoPhyl ----
        orthophyl_output = self.output_dir / "orthophyl_run"
        if self._check_checkpoint(f"orthophyl_{safe}") and self.resume:
            logger.info(f"  ✓ OrthoPhyl already complete for {safe} (resuming)")
        else:
            self._run_orthophyl(
                input_dir=genomes_dir, output_dir=orthophyl_output,
                taxon_name=safe, assemblies=[])
            self._write_checkpoint(f"orthophyl_{safe}")

        # ---- Database ----
        logger.info(f"\nCreating database for {self.clade_name}...")
        if self._check_checkpoint(f"database_{safe}") and self.resume:
            logger.info(f"  ✓ Database already created for {safe} (resuming)")
        else:
            self._create_local_database(
                clade_name=self.clade_name,
                orthophyl_output=orthophyl_output,
                taxonomy=taxonomy,
                taxonomy_source='user_supplied',
                qc_applied=qc_applied,
                genomes_to_keep=genomes_dir,
                source_genome_dir=self.genome_dir)
            self._write_checkpoint(f"database_{safe}")

        # ---- Publish the tree into 03_results ----
        tree_dir = self.results_dir / "trees" / "orthophyl"
        tree_dir.mkdir(parents=True, exist_ok=True)
        tree_file = self._locate_species_tree(orthophyl_output)
        if tree_file.exists():
            shutil.copy(tree_file, tree_dir / f"{safe}_phylogeny.nwk")
            logger.info(f"  ✓ Published tree: {tree_dir / f'{safe}_phylogeny.nwk'}")
        else:
            logger.warning(f"  ⚠ Species tree not found at expected location: {tree_file}")

        logger.info("=" * 70)
        logger.info("LOCAL GENOME-INGEST MODE COMPLETE!")
        logger.info("*** NOTE: clade name/taxonomy is user-supplied, not assigned "
                     "by NCBI ***")
        logger.info("=" * 70)
        logger.info(f"  Database created: {self.clade_name}_db")
        logger.info(f"  Assemblies: {len(kept)}")
        logger.info(f"  QC applied: {qc_applied}")
        logger.info(f"  Taxonomy routable: {is_routable}")
        self._save_final_status()
        return 0

    @staticmethod
    def _accessions_from_dir(d: Path) -> List[str]:
        """Sorted list of genome stems (accessions) from a genomes_to_keep-style dir."""
        return sorted(
            p.stem for p in
            list(d.glob("*.fna")) + list(d.glob("*.fasta")))

    def _create_taxon_database(
        self,
        taxon_name: str,
        orthophyl_output: Path,
        gatherer,
        genomes_to_keep: Path
    ):
        """Create database with taxon metadata.

        Accession metadata is derived from the QC-filtered ``genomes_to_keep/``
        directory produced by gather_filter_asms.sh (each ``.fna``/``.fasta``
        filename is an assembly accession). This records exactly the genomes
        that went into the tree, rather than the pre-QC NCBI candidate list.
        """
        # Get taxonomy from gatherer
        taxonomy = gatherer.get_taxonomy_string()

        # Create database entry
        self._create_database_entry(
            taxon_name=taxon_name,
            orthophyl_output=orthophyl_output,
            taxonomy=taxonomy
        )

        # Update database config with taxon metadata
        db_dir = self.database_dir / f"{taxon_name}_db"
        if db_dir.exists():
            config_file = db_dir / "database_config.json"
            if config_file.exists():
                with open(config_file, 'r') as f:
                    config = json.load(f)

                # Derive the accession list from the QC-filtered genomes.
                accessions = self._accessions_from_dir(genomes_to_keep)

                # Add taxon metadata
                config['source_taxon_name'] = taxon_name
                config['source_taxid'] = gatherer.taxid
                config['source_rank'] = gatherer.taxon_rank
                config['assembly_accessions'] = accessions
                config['n_assemblies_at_creation'] = len(accessions)
                config['last_updated'] = datetime.now().isoformat()

                # Save updated config
                with open(config_file, 'w') as f:
                    json.dump(config, f, indent=2)

                logger.info(f"  ✓ Updated database metadata with taxon info")

    def _create_local_database(
        self,
        clade_name: str,
        orthophyl_output: Path,
        taxonomy: str,
        taxonomy_source: str,
        qc_applied: bool,
        genomes_to_keep: Path,
        source_genome_dir: Path,
    ):
        """Create database for local genome-ingest mode.

        Unlike _create_taxon_database, there is no gatherer/taxid -- the clade
        name/taxonomy is user-supplied (or resolved offline against the local
        taxdump, see _resolve_local_taxonomy). source_taxon_name is still set to
        clade_name so a later `--taxon <same name> --update-existing` run finds
        this database (_check_existing_taxon_database matches on it).
        """
        self._create_database_entry(
            taxon_name=clade_name,
            orthophyl_output=orthophyl_output,
            taxonomy=taxonomy,
            taxonomy_source=taxonomy_source,
            qc_applied=qc_applied,
        )

        db_dir = self.database_dir / f"{clade_name}_db"
        config_file = db_dir / "database_config.json"
        if not config_file.exists():
            return
        with open(config_file, 'r') as f:
            config = json.load(f)

        accessions = self._accessions_from_dir(genomes_to_keep)

        config['source_taxon_name'] = clade_name
        config['source_taxid'] = None
        config['source_rank'] = None
        config['source_genome_dir'] = str(source_genome_dir.resolve())
        config['assembly_accessions'] = accessions
        config['n_assemblies_at_creation'] = len(accessions)
        config['last_updated'] = datetime.now().isoformat()

        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)

        logger.info(f"  ✓ Updated database metadata with local genome-ingest info")
    
    def _update_taxon_database_metadata(
        self,
        db_dir: Path,
        added_accessions: List[str]
    ):
        """Update database metadata after adding new assemblies.

        added_accessions are the accessions that actually entered ReLeaf (i.e.
        survived QC), so the config reflects the true DB contents.
        """
        config_file = db_dir / "database_config.json"
        if not config_file.exists():
            logger.warning(f"  ⚠ Config file not found: {config_file}")
            return

        with open(config_file, 'r') as f:
            config = json.load(f)

        # Update metadata (dedupe in case an accession was somehow already listed)
        existing_accessions = config.get('assembly_accessions', [])
        merged = existing_accessions + [
            a for a in added_accessions if a not in set(existing_accessions)]
        config['assembly_accessions'] = merged
        config['last_updated'] = datetime.now().isoformat()
        config['n_genomes'] = len(config['assembly_accessions'])
        
        # Save updated config
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        logger.info(f"  ✓ Updated database metadata")
        logger.info(f"    Total assemblies: {len(config['assembly_accessions'])}")
    
    def _save_final_status(self):
        """Save final pipeline status."""
        self.pipeline_status['end_time'] = datetime.now().isoformat()
        self.pipeline_status['dry_run'] = self.dry_run
        status_file = self.output_dir / "pipeline_status.json"
        
        with open(status_file, 'w') as f:
            json.dump(self.pipeline_status, f, indent=2)
        
        if self.dry_run:
            logger.info(f"\n[DRY RUN] Pipeline status would be saved to: {status_file}")
        else:
            logger.info(f"\nPipeline status saved: {status_file}")


def main():
    parser = argparse.ArgumentParser(
        description="OrthoPhyl Pipeline Wrapper - Automated routing and phylogenetic placement",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python orthophyl_pipeline_wrapper.py \\
      --input assemblies.tsv \\
      --database-dir databases/ \\
      --output-dir results/ \\
      --threads 32
  
  # With genome downloading
  python orthophyl_pipeline_wrapper.py \\
      --input assemblies.tsv \\
      --database-dir databases/ \\
      --output-dir results/ \\
      --gather-script utils/gather_filter_asms.sh \\
      --threads 32
  
  # Initial setup with database creation
  python orthophyl_pipeline_wrapper.py \\
      --input assemblies.tsv \\
      --database-dir databases/ \\
      --output-dir results/ \\
      --orthophyl-runs orthophyl_runs.tsv \\
      --threads 32
  
  # Resume from checkpoint
  python orthophyl_pipeline_wrapper.py \\
      --input assemblies.tsv \\
      --database-dir databases/ \\
      --output-dir results/ \\
      --resume

  # Build a database from genomes already on disk (QC runs by default; the
  # clade name/taxonomy here is user-supplied, not assigned by NCBI, unless
  # --clade-name happens to resolve against the NCBI taxdump)
  python orthophyl_pipeline_wrapper.py \\
      --genome-dir my_isolates/ \\
      --clade-name MyNovelClade \\
      --database-dir databases/ \\
      --gather-script utils/gather_filter_asms.sh \\
      --threads 32
        """
    )

    # Mode selection: batch mode (--input) vs taxon mode (--taxon)
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        '--input',
        help='Input TSV: assembly_path, taxonomy, [id]. For batch mode.'
    )
    mode_group.add_argument(
        '--taxon',
        help='Taxon name for auto-gather mode (e.g., "Methylorubrum"). For taxon mode.'
    )
    mode_group.add_argument(
        '--genome-dir',
        help='Directory of genome FASTAs already on disk (.fna/.fa/.fasta, optionally '
             '.gz). Builds a tree and database under --clade-name, QC-ing by default. '
             'For local genome-ingest mode.'
    )

    parser.add_argument(
        '--database-dir',
        required=True,
        help='Directory containing *_db databases'
    )
    parser.add_argument(
        '--output-dir',
        help='Output directory for all results (default: <database-dir>/.pipeline_runs/<taxon>_<timestamp>, '
             'where <taxon> is the --taxon name, "TaxID<num>" for a numeric TaxID, or "run" for batch mode)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help='Number of threads (default: 8)'
    )
    parser.add_argument(
        '--gather-script',
        help='Path to gather_filter_asms.sh for genome downloading'
    )
    parser.add_argument(
        '--orthophyl-runs',
        help='TSV for initial database creation (if databases don\'t exist)'
    )
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from last checkpoint'
    )
    parser.add_argument(
        '--skip-download',
        action='store_true',
        help='Skip genome downloading (use existing genomes)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be executed without running anything'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        default=0,
        help='Verbose output: -v shows stdout from subprocesses, -vv shows stdout and stderr'
    )
    parser.add_argument(
        '--low-ram',
        action='store_true',
        help='Use reduced memory mode for CheckM2 (passes --lowmem to gather_filter_asms.sh)'
    )
    parser.add_argument(
        '--use-bbmap',
        action='store_true',
        help='Use bbmap statswrapper instead of CheckM2 for genome statistics (faster, less RAM, but no completeness/contamination filtering)'
    )
    parser.add_argument(
        '--must-keep',
        help='Accessions that MUST survive QC or the run aborts with a clear per-metric '
             'report. Supply either a comma-separated list (e.g. GCF_000...,GCF_001...) '
             'or a path to a file with one accession per line.'
    )
    parser.add_argument(
        '--keep-failing-query',
        action='store_true',
        help='Let query/input genomes that fail QC through with a loud warning instead '
             'of aborting (default: a query genome failing QC aborts the run).'
    )

    # Taxon mode arguments (--taxon is in mutually_exclusive_group above)
    parser.add_argument(
        '--taxon-rank',
        choices=['species', 'genus', 'family', 'order', 'class', 'phylum'],
        help='Taxonomic rank for --taxon query (default: auto-detect)'
    )
    parser.add_argument(
        '--update-existing',
        action='store_true',
        help='Update existing database with new assemblies (taxon mode only)'
    )
    parser.add_argument(
        '--max-tree-genomes',
        type=int,
        default=2000,
        help='Regular single-tree ceiling. When a taxon downloads more raw genomes '
             'than this, the default behavior is to build ONE tree from a diverse '
             'MASH subsample of --subsample-size genomes (greedy max-min; query '
             'genomes are always kept). Default 2000.'
    )
    parser.add_argument(
        '--subsample-size',
        type=int,
        default=500,
        help='Target genome count when a taxon exceeds --max-tree-genomes: a MASH '
             'greedy max-min diverse subset of this size is used to build one tree. '
             'Sketches genomes linearly (no O(n^2) matrix), so it scales to very '
             'large taxa. Default 500.'
    )
    parser.add_argument(
        '--max-total-genomes',
        type=int,
        default=5000,
        help='Guardrail for the opt-in per-subclade partition/megatree path only. '
             'Partitioning builds an all-vs-all distance matrix that grows O(n^2) in '
             'memory (~20 GB at 50k genomes), so it is refused above this ceiling. '
             'The default subsample path never builds the matrix and is unaffected. '
             'Default 5000.'
    )
    parser.add_argument(
        '--megatree',
        action='store_true',
        help='Opt-in large-taxon strategy: instead of subsampling an oversized taxon '
             'to one tree, partition the raw set into size-bounded subclades '
             '(--subclade-size each), build a full tree per subclade, build a small '
             'BACKBONE tree from --backbone-reps diverse reps per subclade, and graft '
             'each subclade tree onto its reps -> one merged tree containing every '
             'genome. High-support bipartition disagreements are flagged (not '
             'resolved). Enforces --max-total-genomes. Overrides the default subsample.'
    )
    parser.add_argument(
        '--backbone-reps',
        type=int,
        default=5,
        help='Megatree only: number of diverse representatives each subclade '
             'contributes to the backbone tree (min(subclade_size, this)); '
             'guarantees every subclade several backbone anchors. Default 5.'
    )
    parser.add_argument(
        '--subclade-size',
        type=int,
        default=150,
        help='Megatree only: per-subclade genome ceiling passed to the partitioner '
             '(--max-size). Distinct from --max-tree-genomes (the single-tree '
             'ceiling). Default 150.'
    )
    parser.add_argument(
        '--conflict-min-support',
        type=int,
        default=90,
        help='Megatree only: support threshold for flagging a bipartition conflict '
             'between a subclade tree and the backbone (0-100 scale, e.g. IQ-TREE '
             'UFBoot). Default 90.'
    )

    # Local genome-ingest mode arguments (--genome-dir is in mutually_exclusive_group above)
    parser.add_argument(
        '--clade-name',
        help='Required with --genome-dir. Names the clade/database. Auto-resolved '
             'against the local NCBI taxdump when it is a real taxon; otherwise the '
             'database is built under this name and NOTE: it was not assigned by '
             'NCBI (use --clade-taxonomy for a full, routable lineage).'
    )
    parser.add_argument(
        '--clade-taxonomy',
        help='Optional escape hatch: a full GTDB taxonomy string '
             '(e.g. "d__Bacteria;p__...;g__MyClade"), used verbatim. Needed only '
             'when --clade-name does not resolve to a known NCBI taxon and you '
             'know the real lineage.'
    )
    parser.add_argument(
        '--clade-rank',
        choices=list(PipelineWrapper._GTDB_RANK_LETTERS),
        default='g',
        help='Rank letter at which an unresolvable --clade-name is attached '
             '(GTDB single-letter rank, d..s). Default: g (genus).'
    )
    parser.add_argument(
        '--skip-qc',
        action='store_true',
        help='Skip CheckM2 QC on --genome-dir genomes (default: QC runs). Use this '
             'only for genomes you have already quality-checked.'
    )

    args = parser.parse_args()

    # Validate argument combinations (mutually_exclusive_group handles --input vs --taxon vs --genome-dir)
    if args.update_existing and not args.taxon:
        parser.error("--update-existing requires --taxon mode.")
    if args.genome_dir and not args.clade_name:
        parser.error("--clade-name is required with --genome-dir.")
    if args.genome_dir and args.megatree:
        parser.error("--megatree is not supported with --genome-dir. Use "
                      "--subsample-size / --max-tree-genomes for large local sets.")

    # Configure logging based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    # Create and run wrapper
    wrapper = PipelineWrapper(
        input_file=args.input,
        database_dir=args.database_dir,
        output_dir=args.output_dir,
        threads=args.threads,
        gather_script=args.gather_script,
        orthophyl_runs_tsv=args.orthophyl_runs,
        resume=args.resume,
        skip_download=args.skip_download,
        dry_run=args.dry_run,
        verbose=args.verbose,
        low_ram=args.low_ram,
        use_bbmap=args.use_bbmap,
        must_keep=args.must_keep,
        keep_failing_query=args.keep_failing_query,
        taxon=args.taxon,
        taxon_rank=args.taxon_rank,
        update_existing=args.update_existing,
        max_tree_genomes=args.max_tree_genomes,
        max_total_genomes=args.max_total_genomes,
        subsample_size=args.subsample_size,
        megatree=args.megatree,
        backbone_reps=args.backbone_reps,
        subclade_size=args.subclade_size,
        conflict_min_support=args.conflict_min_support,
        genome_dir=args.genome_dir,
        clade_name=args.clade_name,
        clade_rank=args.clade_rank,
        clade_taxonomy=args.clade_taxonomy,
        skip_qc=args.skip_qc,
    )

    return wrapper.run()


if __name__ == "__main__":
    sys.exit(main())