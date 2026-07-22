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
        output_dir: Path = None,
        threads: int = 8,
        gather_script: Optional[Path] = None,
        orthophyl_runs_tsv: Optional[Path] = None,
        resume: bool = False,
        skip_download: bool = False,
        dry_run: bool = False,
        verbose: int = 0,
        low_ram: bool = False,
        use_bbmap: bool = False,
        # NEW: Taxon mode parameters
        taxon: Optional[str] = None,
        taxon_rank: Optional[str] = None,
        update_existing: bool = False
    ):
        self.input_file = Path(input_file) if input_file else None
        self.database_dir = Path(database_dir) if database_dir else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.gather_script = Path(gather_script) if gather_script else None
        self.orthophyl_runs_tsv = Path(orthophyl_runs_tsv) if orthophyl_runs_tsv else None
        self.resume = resume
        self.skip_download = skip_download
        self.dry_run = dry_run
        self.verbose = verbose
        self.low_ram = low_ram
        self.use_bbmap = use_bbmap
        
        # NEW: Taxon mode
        self.taxon = taxon
        self.taxon_rank = taxon_rank
        self.update_existing = update_existing
        self.taxon_mode = taxon is not None
        
        # Script paths (relative to this wrapper)
        self.script_dir = Path(__file__).parent
        self.assembly_router = self.script_dir / "assembly_router" / "assembly_router_multi.cmd_out3.py"
        self.database_creator = self.script_dir / "assembly_router" / "create_hierarchical_database_v2.py"
        self.releaf_versioner = self.script_dir / "assembly_router" / "add_releaf_version.py"
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
        
        # Status tracking
        self.pipeline_status = {
            'start_time': datetime.now().isoformat(),
            'phases': {},
            'summary': {}
        }
    
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
        
        # Phase 4: Results aggregation
        self._phase_aggregation()
        
        logger.info("=" * 70)
        logger.info("PIPELINE COMPLETE!")
        logger.info("=" * 70)
        
        self._save_final_status()
        return 0
    
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
        
        # Initialize or validate databases
        if not self.database_dir.exists() or not (self.database_dir / "database_index.json").exists():
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
        
        logger.info(f"\nRouting Summary:")
        logger.info(f"  ReLeaf route: {len(routing_results['releaf_batch'])} assemblies")
        logger.info(f"  OrthoPhyl route: {sum(len(v) for v in routing_results['orthophyl_batch'].values())} assemblies")
        logger.info(f"    ({len(routing_results['orthophyl_batch'])} unique taxa)")
        
        self._write_checkpoint('routing')
        self.pipeline_status['phases']['routing'] = {
            'status': 'complete',
            'releaf_count': len(routing_results['releaf_batch']),
            'orthophyl_count': sum(len(v) for v in routing_results['orthophyl_batch'].values())
        }
        
        return routing_results
    
    def _parse_routing_results(self) -> Dict:
        """Parse routing decision JSON files."""
        releaf_batch = []
        orthophyl_batch = defaultdict(list)
        
        for json_file in self.routing_dir.glob("routing_decision_*.json"):
            with open(json_file, 'r') as f:
                decision = json.load(f)
            
            if decision['pipeline'] == 'ReLeaf':
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
            'orthophyl_batch': dict(orthophyl_batch)
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
            self._run_releaf(
                database_dir=db_dir,
                input_genomes=input_dir,
                output_dir=output_dir,
                tree_method=tree_method,
                tree_data=tree_data,
                database_name=database_name
            )
            
            self._write_checkpoint(checkpoint_name)
        
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
        database_name: str
    ):
        """Execute ReLeaf for a single database."""
        cmd = [
            str(self.releaf_script),
            '--store', str(database_dir / 'orthophyl_run'),
            '--input_genomes', str(input_genomes),
            '-t', str(self.threads),
            '--tree_method', tree_method,
            '--TREE_DATA', tree_data
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
            raise RuntimeError(f"ReLeaf failed for {database_name}. Check log: {log_file}")
        
        logger.info(f"  ✓ ReLeaf complete for {database_name}")
            
        # Create new database version from ReLeaf output
        if self.releaf_versioner.exists():
            self._create_releaf_version(database_name, output_dir)
        else:
            logger.warning(f"  ⚠ ReLeaf versioner not found, skipping version creation")
    
    def _create_releaf_version(self, database_name: str, releaf_output_dir: Path):
        """Create a new database version from ReLeaf output."""
        logger.info(f"\n  Creating new database version from ReLeaf output...")
        
        # Find the database directory
        db_dir = None
        for db_path in self.database_dir.glob("*_db"):
            config_file = db_path / "database_config.json"
            if config_file.exists():
                with open(config_file, 'r') as f:
                    config = json.load(f)
                if config.get('clade_name') == database_name:
                    db_dir = db_path
                    break
        
        if not db_dir:
            logger.warning(f"  ⚠ Could not find database for {database_name}")
            return

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
            
            # Stage 1: Download genomes
            download_dir = self.orthophyl_dir / "downloads" / taxon_name
            
            # Check if download already completed successfully
            download_complete = self._verify_download_complete(download_dir, taxon_name)
            
            if download_complete and self.resume:
                logger.info(f"  ✓ Download already complete (resuming)")
                genomes_to_keep = download_dir / "genomes_to_keep"
                genome_files = list(genomes_to_keep.glob("*.fna")) + list(genomes_to_keep.glob("*.fasta"))
                logger.info(f"    Found {len(genome_files)} genomes")
            elif not self.skip_download:
                self._download_genomes(taxon_name, download_dir)
                self._write_checkpoint(f"download_{taxon_name}")
            else:
                logger.info(f"  Skipping download (--skip-download)")
                # Verify genomes exist if skipping download
                genomes_to_keep = download_dir / "genomes_to_keep"
                if not genomes_to_keep.exists():
                    raise FileNotFoundError(
                        f"--skip-download specified but no genomes found at: {genomes_to_keep}\n"
                        f"Please either:\n"
                        f"  1. Remove --skip-download to allow downloading\n"
                        f"  2. Manually place genomes in {genomes_to_keep}/"
                    )
            
            # Stage 2: Add query genomes
            genomes_to_keep = download_dir / "genomes_to_keep"
            if not genomes_to_keep.exists():
                genomes_to_keep.mkdir(parents=True, exist_ok=True)
            
            logger.info(f"\n  Adding {len(assemblies)} query genomes to input set...")
            for asm in assemblies:
                src = Path(asm['assembly_path'])
                dst = genomes_to_keep / src.name
                if not dst.exists():
                    shutil.copy(src, dst)
                    logger.info(f"    Added: {asm['assembly_id']}")
            
            # Count total genomes
            total_genomes = len(list(genomes_to_keep.glob("*.fna")))
            logger.info(f"  Total genomes for OrthoPhyl: {total_genomes}")
            
            # Stage 3: Run OrthoPhyl
            orthophyl_output = self.orthophyl_dir / "orthophyl_runs" / taxon_name
            if not self._check_checkpoint(f"orthophyl_{taxon_name}") or not self.resume:
                self._run_orthophyl(
                    input_dir=genomes_to_keep,
                    output_dir=orthophyl_output,
                    taxon_name=taxon_name,
                    assemblies=assemblies
                )
                self._write_checkpoint(f"orthophyl_{taxon_name}")
            else:
                logger.info(f"  ✓ OrthoPhyl already complete (resuming)")
            
            # Stage 4: Create database
            if not self._check_checkpoint(f"database_{taxon_name}") or not self.resume:
                self._create_database_entry(
                    taxon_name=taxon_name,
                    orthophyl_output=orthophyl_output,
                    taxonomy=assemblies[0]['download_taxonomy']
                )
                self._write_checkpoint(f"database_{taxon_name}")
            else:
                logger.info(f"  ✓ Database already created (resuming)")
        
        self.pipeline_status['phases']['orthophyl'] = {
            'status': 'complete',
            'taxa_processed': len(orthophyl_batch)
        }
    
    def _download_genomes(self, taxon_name: str, output_dir: Path):
        """Download genomes using gather_filter_asms.sh."""
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
        
        # Add optional flags
        if self.use_bbmap:
            cmd.append('--use-bbmap')
            logger.info(f"  Using bbmap statswrapper instead of CheckM")
        elif self.low_ram:
            cmd.append('--reduced_tree')
            logger.info(f"  Using CheckM --reduced_tree option (low RAM mode)")
        
        logger.info(f"  Downloading genomes for {taxon_name}...")
        if self.verbose:
            logger.info(f"    Command: {' '.join(cmd)}")
        logger.info(f"    Output: {output_dir}")
        
        if self.dry_run:
            logger.info(f"  [DRY RUN] Would download genomes for {taxon_name}")
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
            raise RuntimeError(f"OrthoPhyl failed for {taxon_name}. Check log: {log_file}")
        
        logger.info(f"  ✓ OrthoPhyl complete for {taxon_name}")
        
        # Verify queries are in tree
        tree_file = output_dir / "FINAL_SPECIES_TREES" / "SCO_strict.CDS.iqtree.treefile"
        if tree_file.exists():
            self._verify_queries_in_tree(tree_file, assemblies)
        else:
            logger.warning(f"  ⚠ Tree file not found: {tree_file}")
    
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
        taxonomy: str
    ):
        """Create new database entry from OrthoPhyl run."""
        logger.info(f"\n  Creating database entry for {taxon_name}...")
        
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
        
        # Collect ReLeaf trees
        releaf_trees = []
        for db_dir in self.releaf_dir.glob("*/"):
            tree_file = db_dir / "ReLeaf_results" / "phylogeny_with_new_genomes.nwk"
            if tree_file.exists():
                dst = trees_dir / "releaf" / f"{db_dir.name}_phylogeny.nwk"
                shutil.copy(tree_file, dst)
                releaf_trees.append(db_dir.name)
                logger.info(f"  ✓ Collected ReLeaf tree: {db_dir.name}")
        
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
        """Generate human-readable summary report."""
        report_file = self.results_dir / "pipeline_summary.txt"
        
        # Check for database versions
        database_versions = {}
        if self.releaf_versioner.exists():
            for db_name in releaf_trees:
                for db_path in self.database_dir.glob("*_db"):
                    if db_path.name.replace('_db', '') in db_name:
                        current_link = db_path / "current"
                        if current_link.exists() and current_link.is_symlink():
                            database_versions[db_name] = current_link.readlink().name
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("ORTHOPHYL PIPELINE - SUMMARY REPORT\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Input File: {self.input_file}\n")
            f.write(f"Database Directory: {self.database_dir}\n")
            f.write(f"Output Directory: {self.output_dir}\n\n")
            
            f.write("RESULTS:\n")
            f.write("-" * 70 + "\n\n")
            
            f.write(f"ReLeaf Route (matched databases): {len(releaf_trees)} databases\n")
            for db in sorted(releaf_trees):
                f.write(f"  - {db}\n")
            f.write("\n")
            
            f.write(f"OrthoPhyl Route (novel taxa): {len(orthophyl_trees)} taxa\n")
            for taxon in sorted(orthophyl_trees):
                f.write(f"  - {taxon} (new database created)\n")
            f.write("\n")
            
            f.write("OUTPUT LOCATIONS:\n")
            f.write(f"  Trees: {self.results_dir}/trees/\n")
            f.write(f"  Logs: {self.logs_dir}/\n")
            f.write(f"  Routing decisions: {self.routing_dir}/\n")
            f.write("\n")
            
            f.write("=" * 70 + "\n")
            f.write("All assemblies have been placed in phylogenetic trees!\n")
            f.write("=" * 70 + "\n")
        
        logger.info(f"  Summary report: {report_file}")
    
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
        
        # Query NCBI for assemblies
        logger.info(f"\nQuerying NCBI for {self.taxon} assemblies...")
        gatherer = TaxonAssemblyGatherer(
            taxon=self.taxon,
            rank=self.taxon_rank,
            output_dir=self.output_dir / "taxon_query"
        )

        if self.dry_run:
            logger.info("  [DRY RUN] Would query NCBI and download assemblies")
            logger.info("=" * 70)
            logger.info("PIPELINE COMPLETE (DRY RUN)!")
            logger.info("=" * 70)
            self._save_final_status()
            return 0
        
        # Get assemblies
        assemblies = gatherer.query_ncbi()
        
        if not assemblies:
            logger.error(f"✗ No assemblies found for taxon '{self.taxon}'")
            return 1
        
        logger.info(f"  ✓ Found {len(assemblies)} assemblies")
        
        # Download assemblies
        logger.info(f"\nDownloading {len(assemblies)} assemblies...")
        download_dir = self.output_dir / "downloaded_assemblies"
        download_dir.mkdir(parents=True, exist_ok=True)
        
        gatherer.download_assemblies(assemblies, download_dir)
        
        # Run OrthoPhyl on downloaded assemblies
        logger.info(f"\nRunning OrthoPhyl on {self.taxon} assemblies...")
        orthophyl_output = self.output_dir / "orthophyl_run"
        
        self._run_orthophyl(
            input_dir=download_dir,
            output_dir=orthophyl_output,
            taxon_name=self.taxon,
            assemblies=[]  # No query assemblies in create mode
        )
        
        # Create database with metadata
        logger.info(f"\nCreating database for {self.taxon}...")
        self._create_taxon_database(
            taxon_name=self.taxon,
            orthophyl_output=orthophyl_output,
            gatherer=gatherer,
            assemblies=assemblies
        )
        
        logger.info("=" * 70)
        logger.info("TAXON MODE COMPLETE!")
        logger.info("=" * 70)
        logger.info(f"  Database created: {self.taxon}_db")
        logger.info(f"  Assemblies: {len(assemblies)}")
        
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
        assemblies: List[Dict]
    ):
        """Create database with taxon metadata."""
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
                
                # Add taxon metadata
                config['source_taxon_name'] = taxon_name
                config['source_taxid'] = gatherer.taxid
                config['source_rank'] = gatherer.taxon_rank
                config['assembly_accessions'] = [a['accession'] for a in assemblies]
                config['n_assemblies_at_creation'] = len(assemblies)
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
    
    parser.add_argument(
        '--input',
        required=True,
        help='Input TSV: assembly_path, taxonomy, [id]'
    )
    parser.add_argument(
        '--database-dir',
        required=True,
        help='Directory containing *_db databases'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for all results'
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
        help='Use reduced memory mode for CheckM (passes --reduced_tree to gather_filter_asms.sh)'
    )
    parser.add_argument(
        '--use-bbmap',
        action='store_true',
        help='Use bbmap statswrapper instead of CheckM for genome statistics (faster, less RAM, but no completeness/contamination filtering)'
    )
    
    # NEW: Taxon mode arguments
    parser.add_argument(
        '--taxon',
        help='Taxon name for auto-gather mode (e.g., "Methylorubrum"). Mutually exclusive with --input.'
    )
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
    
    args = parser.parse_args()
    
    # Validate argument combinations
    if args.taxon and args.input:
        parser.error("--taxon and --input are mutually exclusive. Use --taxon for auto-gather mode or --input for batch mode.")
    
    if not args.taxon and not args.input:
        parser.error("Either --taxon or --input must be specified.")
    
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
        taxon=args.taxon,
        taxon_rank=args.taxon_rank,
        update_existing=args.update_existing
    )
    
    return wrapper.run()


if __name__ == "__main__":
    sys.exit(main())