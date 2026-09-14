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
        max_tree_genomes: int = 150
    ):
        # Validate mutually exclusive flags
        if input_file and taxon:
            raise ValueError("Cannot specify both --input and --taxon. Use one or the other.")
        
        self.input_file = Path(input_file) if input_file else None
        self.database_dir = Path(database_dir) if database_dir else None
        
        # Default output_dir to database_dir/.pipeline_runs/<name>_<timestamp> if not provided
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            run_name = self._default_run_name(taxon)
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

        # NEW: subclade partitioning. When a taxon's RAW downloaded genome set
        #   exceeds this ceiling, MASH-partition it into <= max_tree_genomes
        #   subclades and build a tree only for the subclade(s) actually needed.
        self.max_tree_genomes = max_tree_genomes

        # Script paths (relative to this wrapper)
        self.script_dir = Path(__file__).parent
        self.assembly_router = self.script_dir / "assembly_router" / "assembly_router.py"
        self.database_creator = self.script_dir / "assembly_router" / "create_hierarchical_database.py"
        self.releaf_versioner = self.script_dir / "assembly_router" / "add_releaf_version.py"
        self.subclade_partitioner = self.script_dir / "python_scripts" / "subclade_partition.py"
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

    def run(self):
        """Main execution pipeline."""
        try:
            logger.info("=" * 70)
            logger.info("ORTHOPHYL PIPELINE WRAPPER")
            if self.taxon_mode:
                logger.info(f"*** TAXON MODE: {self.taxon} ***")
                if self.update_existing:
                    logger.info("*** UPDATE MODE: Checking for new assemblies ***")
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

        # Phase 3c: lazy subclade-build route (build tree on demand, then ReLeaf)
        if routing_results.get('subclade_build_batch'):
            self._phase_subclade_build(routing_results['subclade_build_batch'])

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
        
        # Initialize or validate databases (skip in taxon create mode)
        if self.taxon_mode and not self.update_existing:
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
    
    def _create_initial_databases(self):
        """Create initial databases from orthophyl_runs.tsv."""
        cmd = [
            'python', str(self.database_creator),
            '--input', str(self.orthophyl_runs_tsv),
            '--output-dir', str(self.database_dir)
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        
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
        if self.verbose:
            logger.info(f"  Command: {' '.join(cmd)}")
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
        
        subclade_build_count = sum(
            len(v) for v in routing_results['subclade_build_batch'].values())
        logger.info(f"\nRouting Summary:")
        logger.info(f"  ReLeaf route: {len(routing_results['releaf_batch'])} assemblies")
        logger.info(f"  OrthoPhyl route: {sum(len(v) for v in routing_results['orthophyl_batch'].values())} assemblies")
        logger.info(f"    ({len(routing_results['orthophyl_batch'])} unique taxa)")
        logger.info(f"  Subclade-build route: {subclade_build_count} assemblies")
        logger.info(f"    ({len(routing_results['subclade_build_batch'])} unbuilt subclades)")

        self._write_checkpoint('routing')
        self.pipeline_status['phases']['routing'] = {
            'status': 'complete',
            'releaf_count': len(routing_results['releaf_batch']),
            'orthophyl_count': sum(len(v) for v in routing_results['orthophyl_batch'].values()),
            'subclade_build_count': subclade_build_count
        }
        
        return routing_results
    
    def _parse_routing_results(self) -> Dict:
        """Parse routing decision JSON files."""
        releaf_batch = []
        orthophyl_batch = defaultdict(list)
        # Unbuilt subclades a query routed to: build the tree on demand, then ReLeaf.
        # Grouped by subclade name so one build serves all queries that landed there.
        subclade_build_batch = defaultdict(list)

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
            elif pipeline == 'OrthoPhyl_subclade_build':
                # Lazy subclade: registered (built=false) at partition time but never
                # built because no query landed there then. This query is the first
                # to route here, so the wrapper builds its tree, then ReLeafs.
                sc_name = decision['subclade_name']
                subclade_build_batch[sc_name].append({
                    'assembly_id': decision['assembly_id'],
                    'assembly_path': decision['assembly'],
                    'subclade_name': sc_name,
                    'parent_taxon': decision.get('parent_taxon'),
                    'subclade_id': decision.get('subclade_id'),
                    'database_dir': decision['database_dir'],
                    'members_file': decision.get('members_file'),
                    'sketch_file': decision.get('sketch_file'),
                    'source_genome_dir': decision.get('source_genome_dir'),
                    'taxonomy': decision.get('query_taxonomy'),
                    'tree_method': decision.get('tree_method', 'iqtree'),
                    'tree_data': decision.get('tree_data', 'CDS'),
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
            'subclade_build_batch': dict(subclade_build_batch)
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
        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")
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
        
        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")
        
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

        # ---- Stage 2: partition (only if over the ceiling) ----
        if raw_count > self.max_tree_genomes:
            logger.info(f"  Raw count {raw_count} > max_tree_genomes "
                        f"{self.max_tree_genomes}: partitioning into subclades")
            manifest = self._partition_genomes(taxon_name, raw_dir, assemblies)
        else:
            # Synthesize a single-subclade (unpartitioned) manifest covering the
            # whole raw set -- the flow below is uniform either way.
            manifest = {
                'partitioned': False, 'parent_taxon': taxon_name,
                'max_size': self.max_tree_genomes, 'n_subclades': 1,
                'subclades': [{'subclade_id': 1, 'name': taxon_name,
                               'n_genomes': raw_count, 'members_file': None,
                               'sketch_file': None}],
                'query_assignments': {Path(a['assembly_path']).name: taxon_name
                                      for a in assemblies},
            }

        assignments = manifest.get('query_assignments', {})
        subclades = manifest['subclades']

        # ---- Stage 3: per-subclade QC + build / lazy-register ----
        if not manifest.get('partitioned'):
            # Whole raw set is one tree (classic behaviour, now with explicit QC).
            entry = subclades[0]
            self._build_subclade(
                taxon_name=taxon_name, entry=entry, raw_dir=raw_dir,
                query_assemblies=assemblies,
                taxonomy=assemblies[0]['download_taxonomy'],
                is_subclade=False)
            return

        # Partitioned: build subclades holding >=1 query; lazily register the rest.
        query_by_subclade: Dict[str, List[Dict]] = {}
        for asm in assemblies:
            name = Path(asm['assembly_path']).name
            sc = assignments.get(name)
            query_by_subclade.setdefault(sc, []).append(asm)

        for entry in subclades:
            sc_name = entry['name']
            queries_here = query_by_subclade.get(sc_name, [])
            if queries_here:
                logger.info(f"\n  Building subclade {sc_name} "
                            f"({entry['n_genomes']} raw genomes, "
                            f"{len(queries_here)} query)")
                self._build_subclade(
                    taxon_name=taxon_name, entry=entry, raw_dir=raw_dir,
                    query_assemblies=queries_here,
                    taxonomy=assemblies[0]['download_taxonomy'],
                    is_subclade=True)
            else:
                logger.info(f"\n  Registering subclade {sc_name} for lazy build "
                            f"({entry['n_genomes']} raw genomes, no query)")
                self._register_lazy_subclade(
                    taxon_name=taxon_name, entry=entry, raw_dir=raw_dir,
                    taxonomy=assemblies[0]['download_taxonomy'])

    def _phase_subclade_build(self, subclade_build_batch: Dict[str, List[Dict]]):
        """Phase 3c: build lazily-registered subclades on demand, then ReLeaf.

        Each key is an unbuilt subclade (registered built=false at partition time)
        that a query has now routed to. For each we:
          1. Re-stage the subclade's raw member FASTAs (recorded source_genome_dir).
          2. QC + run OrthoPhyl on those raw members and promote the DB entry to
             built=true (force-overwriting the placeholder). The query is NOT part
             of this build -- the tree is the subclade's own genomes.
          3. ReLeaf the waiting query assemblies onto the freshly-built tree.
        A failure in one subclade is logged and skipped so others still proceed.
        """
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 3C: SUBCLADE-BUILD ROUTE (Lazy Subclades)")
        logger.info("=" * 70)
        logger.info(f"Building {len(subclade_build_batch)} unbuilt subclades on demand")

        for sc_name, queries in subclade_build_batch.items():
            logger.info(f"\n{'=' * 60}")
            logger.info(f"Subclade: {sc_name} "
                        f"({len(queries)} query assemblies waiting)")
            logger.info(f"{'=' * 60}")
            try:
                self._process_subclade_build(sc_name, queries)
            except Exception as e:
                logger.error(f"✗ Subclade build failed for {sc_name}: {e}")
                # Continue with other subclades rather than failing the whole run.
                continue

        self.pipeline_status['phases']['subclade_build'] = {
            'status': 'complete',
            'subclades_processed': len(subclade_build_batch)
        }

    def _process_subclade_build(self, sc_name: str, queries: List[Dict]):
        """Build one lazy subclade from its raw members, then ReLeaf the queries."""
        first = queries[0]
        source_genome_dir = first.get('source_genome_dir')
        if not source_genome_dir:
            raise RuntimeError(
                f"Subclade {sc_name} has no source_genome_dir recorded; cannot "
                f"locate its raw members to build. Was it registered lazily?")
        raw_dir = Path(source_genome_dir)
        if not self.dry_run and not raw_dir.exists():
            raise FileNotFoundError(
                f"Raw member directory for subclade {sc_name} not found: {raw_dir}")

        parent_taxon = first.get('parent_taxon')
        taxonomy = first.get('taxonomy')

        # Reconstruct the partition entry _build_subclade expects. members_file /
        # sketch_file point at the in-DB copies written at lazy registration.
        entry = {
            'name': sc_name,
            'subclade_id': first.get('subclade_id'),
            'members_file': first.get('members_file'),
            'sketch_file': first.get('sketch_file'),
        }

        # ---- Build the subclade tree from ITS OWN genomes (no query genomes). ----
        # force=True promotes the built=false placeholder DB to a real built entry.
        logger.info(f"\n  Building subclade {sc_name} from raw members in {raw_dir}")
        self._build_subclade(
            taxon_name=parent_taxon or sc_name,
            entry=entry,
            raw_dir=raw_dir,
            query_assemblies=[],
            taxonomy=taxonomy,
            is_subclade=True,
            force=True)

        # ---- ReLeaf the waiting queries onto the freshly-built subclade. ----
        db_dir = Path(first['database_dir'])
        tree_method = first.get('tree_method', 'iqtree')
        tree_data = first.get('tree_data', 'CDS')

        input_dir = self.releaf_dir / sc_name / "input_genomes"
        if not self.dry_run:
            input_dir.mkdir(parents=True, exist_ok=True)
            for asm in queries:
                src = Path(asm['assembly_path'])
                dst = input_dir / f"{asm['assembly_id']}.fna"
                if not dst.exists():
                    shutil.copy(src, dst)
                logger.info(f"  Prepared query for ReLeaf: {asm['assembly_id']}")

        logger.info(f"\n  ReLeaf {len(queries)} query assemblies onto {sc_name}")
        self._run_releaf(
            database_dir=db_dir,
            input_genomes=input_dir,
            output_dir=self.releaf_dir / sc_name,
            tree_method=tree_method,
            tree_data=tree_data,
            database_name=sc_name,
            n_assemblies=len(queries))

    def _build_subclade(self, taxon_name: str, entry: Dict, raw_dir: Path,
                        query_assemblies: List[Dict], taxonomy: str,
                        is_subclade: bool, force: bool = False):
        """QC one subclade's raw members, run OrthoPhyl, and create its DB entry.

        For the unpartitioned case (is_subclade=False) sc_name == taxon_name and
        the whole raw set is the member list. CheckM2 QC runs HERE (deferred from
        partition time), only on this subclade's members.

        force=True overwrites an existing DB dir for this subclade -- required when
        building a subclade that was previously registered lazily (built=false), so
        the placeholder entry is replaced by the real built=true one.
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

    def _register_lazy_subclade(self, taxon_name: str, entry: Dict,
                                raw_dir: Path, taxonomy: str):
        """Register a subclade as built=false (no QC, no tree) via the DB creator.

        Records its raw member list + sketch + source_genome_dir so a future query
        routing here can QC + build it on demand.
        """
        sc_name = entry['name']
        sc_id = entry.get('subclade_id')
        ckey = sc_name

        if self._check_checkpoint(f"database_{ckey}") and self.resume:
            logger.info(f"  ✓ Lazy subclade already registered for {sc_name} (resuming)")
            return

        cmd = [
            'python', str(self.database_creator),
            '--single-clade', sc_name, taxonomy, str(raw_dir),
            '--output-dir', str(self.database_dir),
            '--is-subclade',
            '--parent-taxon', taxon_name,
            '--register-only',
            '--n-genomes', str(entry.get('n_genomes', 0)),
            '--source-genome-dir', str(raw_dir),
        ]
        if sc_id is not None:
            cmd.extend(['--subclade-id', str(sc_id)])
        if entry.get('sketch_file'):
            cmd.extend(['--sketch-file', str(entry['sketch_file'])])
        if entry.get('members_file'):
            cmd.extend(['--members-file', str(entry['members_file'])])

        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")
        if self.dry_run:
            logger.info(f"  [DRY RUN] Would register lazy subclade {sc_name}")
            return

        log_file = self.logs_dir / f"database_{sc_name}.log"
        with open(log_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"Lazy registration failed for {sc_name}. Check log: {log_file}")
        logger.info(f"  ✓ Registered lazy subclade: {sc_name}_db (built=false)")
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
        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")
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
    # Subclade partitioning helpers (pre-QC MASH partition -> lazy build)
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

    def _partition_genomes(self, taxon_name: str, raw_genome_dir: Path,
                           query_assemblies: List[Dict]) -> Dict:
        """Run subclade_partition.py on the RAW genome set; return the manifest dict.

        Builds <orthophyl_dir>/partitions/<taxon>/ holding MASH_out, per-subclade
        .msh/.members.txt, and partition_manifest.json. Every query stem must land
        in exactly one subclade (it is a clustering leaf), which we assert.
        """
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
            '--max-size', str(self.max_tree_genomes),
            '--threads', str(self.threads),
        ]
        for asm in query_assemblies:
            cmd.extend(['--query', Path(asm['assembly_path']).name])

        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")
        if self.dry_run:
            logger.info(f"  [DRY RUN] Would partition {taxon_name} into subclades")
            # Synthesize a trivial single-subclade manifest for dry-run flow.
            return {
                'partitioned': False, 'parent_taxon': taxon_name,
                'max_size': self.max_tree_genomes, 'n_subclades': 1,
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
        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")
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
        force: bool = False
    ):
        """Create new database entry from OrthoPhyl run.

        When subclade_meta is given, the entry is written as a built subclade
        (is_subclade + parent_taxon + sketch/members recorded) via the DB creator's
        --single-clade path; otherwise the classic TSV --update path is used.

        force=True passes --force so an existing DB dir is overwritten. Needed when
        building a subclade that was previously registered lazily (built=false):
        without it the DB creator refuses (FileExistsError) and the placeholder
        entry survives instead of being promoted to built=true.
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
            ]
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

            if self.verbose:
                logger.info(f"    Command: {' '.join(cmd)}")
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
            '--update'
        ]

        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")

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

        return {'releaf_batch': releaf_batch, 'orthophyl_batch': dict(orthophyl_batch),
                'subclade_build_batch': {}}
    
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

        # Partition if over the ceiling; create mode builds ALL subclades.
        if raw_count > self.max_tree_genomes:
            logger.info(f"  Raw count {raw_count} > max_tree_genomes "
                        f"{self.max_tree_genomes}: partitioning into subclades")
            manifest = self._partition_genomes(self.taxon, raw_dir, query_assemblies=[])
        else:
            manifest = {
                'partitioned': False, 'parent_taxon': self.taxon,
                'max_size': self.max_tree_genomes, 'n_subclades': 1,
                'subclades': [{'subclade_id': 1, 'name': self.taxon,
                               'n_genomes': raw_count, 'members_file': None,
                               'sketch_file': None}],
                'query_assignments': {},
            }

        taxonomy = gatherer.get_taxonomy_string()

        if not manifest.get('partitioned'):
            # Single tree: QC the whole raw set then build + taxon-flavored DB.
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

        # Partitioned: build every subclade (create mode has no single query).
        logger.info(f"  Building all {manifest['n_subclades']} subclades")
        for entry in manifest['subclades']:
            logger.info(f"\n  Building subclade {entry['name']} "
                        f"({entry['n_genomes']} raw genomes)")
            self._build_subclade(
                taxon_name=self.taxon, entry=entry, raw_dir=raw_dir,
                query_assemblies=[], taxonomy=taxonomy, is_subclade=True)

        logger.info("=" * 70)
        logger.info("TAXON MODE COMPLETE!")
        logger.info("=" * 70)
        logger.info(f"  Subclade databases created: {manifest['n_subclades']}")
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
        
        # Download new assemblies
        logger.info(f"\nDownloading {len(new_assemblies)} new assemblies...")
        download_dir = self.output_dir / "new_assemblies"
        download_dir.mkdir(parents=True, exist_ok=True)
        
        # TODO: New assemblies downloaded here via TaxonAssemblyGatherer are NOT yet
        # QC-filtered (completeness/contamination/N50). Harmonize with gather_filter_asms.sh
        # filtering in a future update. See create-mode for the filtered path.
        gatherer.download_assemblies(new_assemblies, download_dir)
        
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
            input_genomes=download_dir,
            output_dir=releaf_output,
            tree_method=tree_method,
            tree_data=tree_data,
            database_name=existing_db['clade_name']
        )
        
        # Update database metadata
        logger.info(f"\nUpdating database metadata...")
        self._update_taxon_database_metadata(
            db_dir=db_dir,
            new_assemblies=new_assemblies
        )
        
        logger.info("=" * 70)
        logger.info("TAXON UPDATE COMPLETE!")
        logger.info("=" * 70)
        logger.info(f"  Database: {existing_db['clade_name']}")
        logger.info(f"  Added assemblies: {len(new_assemblies)}")
        logger.info(f"  Total assemblies: {len(existing_accessions) + len(new_assemblies)}")
        
        self._save_final_status()
        return 0
    
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
                accessions = sorted(
                    p.stem for p in
                    list(genomes_to_keep.glob("*.fna")) + list(genomes_to_keep.glob("*.fasta"))
                )

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
    
    def _update_taxon_database_metadata(
        self,
        db_dir: Path,
        new_assemblies: List[Dict]
    ):
        """Update database metadata after adding new assemblies."""
        config_file = db_dir / "database_config.json"
        if not config_file.exists():
            logger.warning(f"  ⚠ Config file not found: {config_file}")
            return
        
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Update metadata
        existing_accessions = config.get('assembly_accessions', [])
        new_accessions = [a['accession'] for a in new_assemblies]
        config['assembly_accessions'] = existing_accessions + new_accessions
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
        default=150,
        help='Maximum genomes per tree/subclade. When a taxon downloads more raw '
             'genomes than this, MASH-partition them into size-bounded subclades '
             '(<Taxon>_1, <Taxon>_2, ...); build a tree only for the subclade(s) '
             'containing a query (others are registered for lazy build). Default 150.'
    )

    args = parser.parse_args()
    
    # Validate argument combinations (mutually_exclusive_group handles --input vs --taxon)
    if args.update_existing and not args.taxon:
        parser.error("--update-existing requires --taxon mode.")
    
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
        max_tree_genomes=args.max_tree_genomes
    )

    return wrapper.run()


if __name__ == "__main__":
    sys.exit(main())