#!/usr/bin/env python3
"""
Multi-Database Assembly Router - Query multiple databases and find best match

This script queries all available databases in a directory and routes assemblies
to the best matching database (ReLeaf) or suggests creating a new OrthoPhyl run
for novel taxa with commands using gather_filter_asms.sh.

Usage:
    python assembly_router.py \
        --assembly genome.fna \
        --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \
        --database-dir /path/to/databases/ \
        --gather-filter-script utils/gather_filter_asms.sh \
        --output-dir results/

Author: Generated for HGTool
Date: 2024
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from datetime import datetime
import re

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# MASH sketch parameters for subclade tie-break comparisons -- MUST match
#   script_lib/functions.sh and python_scripts/subclade_partition.py, or
#   distances between independently-produced sketches are meaningless.
MASH_K = "17"
MASH_S = "5000"


def pick_nearest_subclade(candidates: List[Dict], dist_fn) -> Tuple[Optional[Dict], Optional[float]]:
    """
    Pick the candidate dict with the smallest dist_fn(candidate) value.

    Pure function (no mash calls) so it is directly unit-testable by injecting
    a distance oracle -- mirrors greedy_maxmin's dist_fn pattern in
    python_scripts/subsample_genomes.py.

    dist_fn(candidate) -> float distance, or None if that candidate is unusable
    (e.g. missing/unreadable sketch) and should be skipped.

    Ties are broken by ascending `clade_name` for determinism across
    filesystems/runs (candidates are visited in that order and only a STRICT
    improvement replaces the current best).

    Returns (winner, distance), or (None, None) if no candidate is usable.
    """
    best = None
    best_dist = None
    for c in sorted(candidates, key=lambda c: c.get('clade_name') or ''):
        d = dist_fn(c)
        if d is None:
            continue
        if best is None or d < best_dist:
            best = c
            best_dist = d
    return best, best_dist


class GTDBTaxonomy:
    """Parse and manipulate GTDB taxonomy strings."""
    
    RANKS = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    RANK_NAMES = {
        'd': 'domain',
        'p': 'phylum',
        'c': 'class',
        'o': 'order',
        'f': 'family',
        'g': 'genus',
        's': 'species'
    }
    
    def __init__(self, taxonomy_string: str):
        self.raw = taxonomy_string.strip()
        self.levels = self._parse_taxonomy(self.raw)
    
    def _parse_taxonomy(self, taxonomy_string: str) -> Dict[str, str]:
        """Parse taxonomy into dictionary."""
        levels = {}
        parts = taxonomy_string.split(';')
        
        for part in parts:
            part = part.strip()
            if not part:
                continue
            
            match = re.match(r'^([dpcofgs])__(.*)$', part)
            if match:
                rank = match.group(1)
                name = match.group(2).strip()
                levels[rank] = name if name else None
        
        return levels
    
    def get_rank(self, rank: str) -> Optional[str]:
        return self.levels.get(rank)
    
    def get_most_specific_rank(self) -> Optional[str]:
        for rank in reversed(self.RANKS):
            if self.levels.get(rank):
                return rank
        return None
    
    def get_rank_name(self, rank: str) -> str:
        return self.RANK_NAMES.get(rank, rank)
    
    def is_within_clade(self, clade: 'GTDBTaxonomy', clade_rank: str) -> bool:
        """Check if this taxonomy belongs to specified clade."""
        for rank in self.RANKS:
            if rank == clade_rank:
                return self.get_rank(rank) == clade.get_rank(rank)
            
            rank_idx = self.RANKS.index(rank)
            clade_idx = self.RANKS.index(clade_rank)
            
            if rank_idx < clade_idx:
                if self.get_rank(rank) != clade.get_rank(rank):
                    return False
        return True
    
    def get_taxonomy_string_at_rank(self, rank: str) -> str:
        """Get taxonomy string up to specified rank."""
        parts = []
        for r in self.RANKS:
            value = self.get_rank(r)
            parts.append(f"{r}__{value if value else ''}")
            if r == rank:
                break
        return ';'.join(parts)
    
    def __str__(self):
        return self.raw


class MultiDatabaseRouter:
    """Routes assemblies by querying multiple databases."""
    
    def __init__(
        self,
        database_dir: Path,
        output_dir: Path,
        gather_filter_script: Optional[Path] = None,
        threads: int = 8,
        placement: str = 'subclade',
    ):
        self.database_dir = Path(database_dir)
        self.output_dir = Path(output_dir)
        self.gather_filter_script = Path(gather_filter_script) if gather_filter_script else None
        self.threads = threads
        # 'subclade' (default): prefer the dense per-subclade tree for best local
        #   resolution. 'backbone': prefer the sparse megatree overview tree.
        self.placement = placement

        # Load all databases
        self.databases = self._load_databases()
        
        if not self.databases:
            raise ValueError(f"No databases found in {database_dir}")
        
        logger.info(f"Loaded {len(self.databases)} databases")
        for db in self.databases:
            note = " [user-supplied taxonomy: not NCBI-assigned]" if db['taxonomy_source'] == 'user_supplied' else ""
            logger.info(f"  - {db['clade_name']} ({db['rank_name']}, {db['n_genomes']} genomes){note}")
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup file logging
        log_file = self.output_dir / f"routing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logger.addHandler(file_handler)
    
    def _load_databases(self) -> List[Dict]:
        """Load all database configs from directory."""
        databases = []
        
        # Check for database_index.json
        index_file = self.database_dir / "database_index.json"
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)
            for db_info in index['databases']:
                db_dir = Path(db_info['database_dir'])
                if db_dir.exists():
                    databases.append(self._load_single_database(db_dir))
        else:
            # Scan for *_db directories
            for db_dir in self.database_dir.glob("*_db"):
                config_file = db_dir / "database_config.json"
                if config_file.exists():
                    databases.append(self._load_single_database(db_dir))
        
        return databases
    
    def _load_single_database(self, db_dir: Path) -> Dict:
        """Load a single database config."""
        config_file = db_dir / "database_config.json"
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        tax = GTDBTaxonomy(config['clade_taxonomy'])

        return {
            'db_dir': db_dir,
            'clade_name': config['clade_name'],
            'clade_taxonomy': config['clade_taxonomy'],
            'clade_rank': config['clade_rank'],
            'rank_name': config['clade_rank_name'],
            'n_genomes': config['n_genomes'],
            'taxonomy_obj': tax,
            # Subclade fields (backward-compatible defaults for pre-subclade DBs).
            'is_subclade': config.get('is_subclade', False),
            'parent_taxon': config.get('parent_taxon'),
            'subclade_id': config.get('subclade_id'),
            'sketch_file': config.get('sketch_file'),
            'members_file': config.get('members_file'),
            'source_genome_dir': config.get('source_genome_dir'),
            'built': config.get('built', True),
            'is_backbone': config.get('is_backbone', False),
            'taxonomy_source': config.get('taxonomy_source', 'ncbi'),
            'qc_applied': config.get('qc_applied', True),
        }
    
    def find_matching_databases(self, query_taxonomy: str) -> List[Dict]:
        """Find all databases that match the query taxonomy."""
        query_tax = GTDBTaxonomy(query_taxonomy)
        matches = []
        
        for db in self.databases:
            clade_tax = db['taxonomy_obj']
            clade_rank = db['clade_rank']
            
            if query_tax.is_within_clade(clade_tax, clade_rank):
                matches.append({
                    'database': db,
                    'specificity': GTDBTaxonomy.RANKS.index(clade_rank)  # Higher = more specific
                })
        
        # Sort by specificity (most specific first)
        matches.sort(key=lambda x: x['specificity'], reverse=True)
        
        return matches
    
    def route_assembly(
        self,
        assembly_path: Path,
        taxonomy: str,
        assembly_id: Optional[str] = None
    ) -> Dict:
        """Route assembly to best matching database or suggest new OrthoPhyl run."""
        if assembly_id is None:
            assembly_id = Path(assembly_path).stem
        
        logger.info(f"\nRouting assembly: {assembly_id}")
        logger.info(f"Query taxonomy: {taxonomy}")
        
        # Find matching databases
        matches = self.find_matching_databases(taxonomy)
        
        if matches:
            # Use most specific match
            best_match = matches[0]['database']

            # Megatree tie-break: several DBs can share the SAME taxonomy string
            # -- the backbone plus its subclades <Taxon>_1, _2, ... -- so they all
            # match at the same top specificity and taxonomy alone cannot tell
            # them apart. Discriminate by --placement, then by MASH distance
            # among same-parent subclades.
            top_spec = matches[0]['specificity']
            top_dbs = [m['database'] for m in matches if m['specificity'] == top_spec]
            subclade_dbs = [d for d in top_dbs if d.get('is_subclade')]
            backbone_dbs = [d for d in top_dbs if d.get('is_backbone')]

            if subclade_dbs or backbone_dbs:
                if self.placement == 'backbone' and backbone_dbs:
                    best_match = backbone_dbs[0]
                    logger.info(f"  --placement backbone: routing to backbone {best_match['clade_name']}")
                elif subclade_dbs:
                    if len(subclade_dbs) > 1:
                        chosen, dist = self._route_subclade_by_mash(assembly_path, subclade_dbs)
                        if chosen is not None:
                            best_match = chosen
                        elif backbone_dbs:
                            logger.warning("  No subclade sketch usable for MASH tie-break; "
                                            "falling back to backbone")
                            best_match = backbone_dbs[0]
                        else:
                            best_match = subclade_dbs[0]
                    else:
                        best_match = subclade_dbs[0]
                elif backbone_dbs:
                    # placement=='subclade' but only a backbone is available.
                    best_match = backbone_dbs[0]

            logger.info(f"✓ MATCH FOUND: {best_match['clade_name']}")
            logger.info(f"  Rank: {best_match['rank_name']}")
            logger.info(f"  Database: {best_match['db_dir'].name}")
            logger.info(f"  Genomes: {best_match['n_genomes']}")

            if len(matches) > 1:
                logger.info(f"  (Also matches {len(matches)-1} other database(s) at broader levels)")

            # An unbuilt subclade cannot serve ReLeaf yet -- its tree does not
            # exist. Emit a build-then-releaf decision the wrapper will execute.
            if best_match.get('is_subclade') and not best_match.get('built', True):
                return self._route_to_subclade_build(
                    assembly_path, assembly_id, taxonomy, best_match)

            return self._route_to_releaf(assembly_path, assembly_id, taxonomy, best_match)
        else:
            logger.info("✗ NO MATCH FOUND in any database")
            logger.info("  → New OrthoPhyl run required")
            
            return self._route_to_orthophyl(assembly_path, assembly_id, taxonomy)
    
    def _route_to_releaf(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str,
        database: Dict
    ) -> Dict:
        """Generate ReLeaf routing decision."""
        # Prepare ReLeaf input
        releaf_input_dir = self.output_dir / "releaf_input"
        releaf_input_dir.mkdir(parents=True, exist_ok=True)
        
        import shutil
        dest = releaf_input_dir / f"{assembly_id}.fna"
        if not dest.exists():
            shutil.copy(assembly_path, dest)
        
        # Read database config to check available options
        config_file = database['db_dir'] / 'database_config.json'
        available_methods = ['iqtree']  # default
        available_data_types = ['CDS']  # default
        
        if config_file.exists():
            with open(config_file, 'r') as f:
                config = json.load(f)
            available_methods = config.get('available_tree_methods', ['iqtree'])
            available_data_types = config.get('available_data_types', ['CDS'])
        
        # Choose best available options (prefer iqtree and CDS)
        tree_method = 'iqtree' if 'iqtree' in available_methods else available_methods[0] if available_methods else 'iqtree'
        tree_data = 'CDS' if 'CDS' in available_data_types else available_data_types[0] if available_data_types else 'CDS'
        
        # Build ReLeaf command with available options
        releaf_cmd = (
            f"./ReLeaf.sh \\\n"
            f"    --store {database['db_dir'] / 'orthophyl_run'} \\\n"
            f"    --input_genomes {releaf_input_dir} \\\n"
            f"    -t {self.threads} \\\n"
            f"    --tree_method {tree_method} \\\n"
            f"    --TREE_DATA {tree_data}"
        )
        
        # Add note if using non-default options
        if tree_method != 'iqtree' or tree_data != 'CDS':
            releaf_cmd += f"\n\n# Note: Using tree_method={tree_method}, TREE_DATA={tree_data}\n"
            releaf_cmd += f"# Available in database: methods={available_methods}, data_types={available_data_types}"
        
        decision = {
            'pipeline': 'ReLeaf',
            'reason': f"Taxonomy matches {database['clade_name']} at {database['rank_name']} level",
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': taxonomy,
            'matched_database': database['clade_name'],
            'matched_rank': database['rank_name'],
            'database_dir': str(database['db_dir']),
            'database_genomes': database['n_genomes'],
            'command': releaf_cmd,
            'tree_method': tree_method,
            'tree_data': tree_data,
            'available_methods': available_methods,
            'available_data_types': available_data_types
        }
        
        self._save_decision(decision)
        return decision

    # ------------------------------------------------------------------ #
    # Subclade routing (MASH sequence distance among same-parent subclades)
    # ------------------------------------------------------------------ #

    def _run_mash(self, cmd) -> str:
        """Run a mash command (shell=False) and return its stdout. Wrapped for mocking."""
        import subprocess
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE,
                                universal_newlines=True)
        return result.stdout

    def _route_subclade_by_mash(
        self, assembly_path: Path, subclade_dbs: List[Dict]
    ) -> Tuple[Optional[Dict], Optional[float]]:
        """
        Pick the subclade DB whose member set contains the query's NEAREST genome.

        Sketches the query with the SAME params as the per-subclade sketches
        (mash -k 17 -s 5000), then `mash dist query.msh <subclade sketch>` for each
        candidate. Each sketch is multi-genome, so `mash dist` emits one line per
        member; we take the MIN distance over members (nearest member) and choose
        the subclade with the smallest such distance via pick_nearest_subclade.
        Returns (None, None) if no distance could be computed (caller keeps its
        default).
        """
        import tempfile
        import shutil as _shutil

        tmp_dir = Path(tempfile.mkdtemp(prefix="router_mash_"))
        try:
            query_prefix = tmp_dir / "query"
            self._run_mash([
                "mash", "sketch", "-k", MASH_K, "-s", MASH_S,
                "-p", str(self.threads), "-o", str(query_prefix), str(assembly_path),
            ])
            query_msh = str(query_prefix) + ".msh"

            def dist_fn(db):
                sketch = db.get('sketch_file')
                if not sketch or not Path(sketch).exists():
                    logger.warning(f"  Subclade {db['clade_name']} has no sketch; skipping in MASH tie-break")
                    return None
                out = self._run_mash(["mash", "dist", query_msh, str(sketch)])
                # Each line: <ref> <query> <dist> <p-value> <shared-hashes>
                dmin = None
                for line in out.splitlines():
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    try:
                        d = float(parts[2])
                    except ValueError:
                        continue
                    if dmin is None or d < dmin:
                        dmin = d
                if dmin is not None:
                    logger.info(f"  MASH nearest-member dist to {db['clade_name']}: {dmin:.4f}")
                return dmin

            best_db, best_dist = pick_nearest_subclade(subclade_dbs, dist_fn)
            if best_db is not None:
                logger.info(f"  → Subclade selected by MASH: {best_db['clade_name']} (dist {best_dist:.4f})")
            return best_db, best_dist
        finally:
            _shutil.rmtree(tmp_dir, ignore_errors=True)

    def _route_to_subclade_build(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str,
        database: Dict
    ) -> Dict:
        """
        Emit a decision to lazily build an unbuilt subclade before ReLeaf.

        The chosen subclade was registered at partition time (built=false) with a
        member list + sketch but no tree. The wrapper executes this by QC-ing +
        running OrthoPhyl on the subclade's raw members, then routing the query
        through ReLeaf against the freshly-built database.
        """
        decision = {
            'pipeline': 'OrthoPhyl_subclade_build',
            'reason': (f"Query nearest to subclade {database['clade_name']} "
                       f"(parent {database.get('parent_taxon')}), which is registered "
                       f"but not yet built; build its tree then ReLeaf"),
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': taxonomy,
            'matched_database': database['clade_name'],
            'matched_rank': database['rank_name'],
            'database_dir': str(database['db_dir']),
            'parent_taxon': database.get('parent_taxon'),
            'subclade_id': database.get('subclade_id'),
            'subclade_name': database['clade_name'],
            # The subclade's OWN recorded taxonomy -- NOT the query's -- so a
            # rebuild does not silently overwrite the subclade's taxonomy with
            # whatever a single query happened to carry.
            'subclade_taxonomy': database['clade_taxonomy'],
            'members_file': database.get('members_file'),
            'sketch_file': database.get('sketch_file'),
            'source_genome_dir': database.get('source_genome_dir'),
            'command': (
                f"# This subclade is registered but not yet built.\n"
                f"# The pipeline wrapper will QC + build its tree, then run ReLeaf.\n"
                f"# Members: {database.get('members_file')}\n"
                f"# Source genomes: {database.get('source_genome_dir')}"
            ),
        }
        self._save_decision(decision)
        return decision

    def _route_to_orthophyl(
        self,
        assembly_path: Path,
        assembly_id: str,
        taxonomy: str
    ) -> Dict:
        """Generate OrthoPhyl routing decision with gather_filter_asms.sh."""
        query_tax = GTDBTaxonomy(taxonomy)
        specific_rank = query_tax.get_most_specific_rank()
        
        # Determine what taxonomic level to download
        download_rank = specific_rank
        download_value = query_tax.get_rank(specific_rank) if specific_rank else None
        
        # If at species level and species is undefined, go up to genus
        if specific_rank == 's' and not download_value:
            download_rank = 'g'
            download_value = query_tax.get_rank('g')
        
        # If still no value, go up to family
        if not download_value and download_rank in ['s', 'g']:
            download_rank = 'f'
            download_value = query_tax.get_rank('f')
        
        rank_name = query_tax.get_rank_name(download_rank) if download_rank else 'unknown'
        
        # Prepare OrthoPhyl input directory
        orthophyl_input_dir = self.output_dir / "orthophyl_new_genomes"
        orthophyl_input_dir.mkdir(parents=True, exist_ok=True)
        
        import shutil
        dest = orthophyl_input_dir / Path(assembly_path).name
        if not dest.exists():
            shutil.copy(assembly_path, dest)
        
        # Build taxonomy string for gather_filter_asms.sh
        # The script can use taxon name directly
        
        # Build commands
        commands = []
        
        # Step 1: Download and filter related genomes using gather_filter_asms.sh
        if self.gather_filter_script and self.gather_filter_script.exists():
            # gather_filter_asms.sh usage: script taxon output_dir threads
            # The taxon can be a taxonomy string or taxon name
            gather_cmd = (
                f"# Step 1: Download and filter related genomes from NCBI\n"
                f"# This will download genomes, run CheckM QC, and filter by quality\n"
                f"{self.gather_filter_script} \\\n"
                f"    {download_value} \\\n"
                f"    {orthophyl_input_dir} \\\n"
                f"    {self.threads}\n"
                f"\n"
                f"# Note: The script will create {orthophyl_input_dir}/genomes_to_keep/\n"
                f"# with high-quality, non-redundant genomes.\n"
                f"# Quality filters (adjust in script if needed):\n"
                f"#   - Completeness >= 95%\n"
                f"#   - Contamination <= 1.0%\n"
                f"#   - Duplication <= 2%\n"
                f"#   - Removes RefSeq/GenBank redundancy\n"
                f"#   - Filters by assembly stats (N50, GC content, length)"
            )
            commands.append(gather_cmd)
            
            # Update orthophyl input directory to the filtered genomes
            orthophyl_input_actual = f"{orthophyl_input_dir}/genomes_to_keep"
        else:
            gather_cmd = (
                f"# Step 1: Download related genomes manually\n"
                f"# Search NCBI for: {download_value} ({rank_name})\n"
                f"# Recommended: Use utils/gather_filter_asms.sh for automatic download & QC\n"
                f"# Place genome files in: {orthophyl_input_dir}/"
            )
            commands.append(gather_cmd)
            orthophyl_input_actual = str(orthophyl_input_dir)
        
        # Step 2: Run OrthoPhyl
        orthophyl_cmd = (
            f"# Step 2: Run OrthoPhyl on expanded genome set\n"
            f"# Note: Your query genome ({Path(assembly_path).name}) should be in {orthophyl_input_actual}\n"
            f"./OrthoPhyl.sh \\\n"
            f"    -g {orthophyl_input_actual} \\\n"
            f"    -o {self.output_dir / 'orthophyl_output'} \\\n"
            f"    -t {self.threads} \\\n"
            f"    --tree_method iqtree \\\n"
            f"    --TREE_DATA CDS \\\n"
            f"    --use_partitions true"
        )
        commands.append(orthophyl_cmd)
        
        # Step 3: Create database
        tax_string_for_download = query_tax.get_taxonomy_string_at_rank(download_rank)
        db_cmd = (
            f"# Step 3: Create new database entry\n"
            f"# Add this line to your orthophyl_runs.tsv:\n"
            f"{download_value}\t{self.output_dir / 'orthophyl_output'}\t{tax_string_for_download}\n"
            f"\n"
            f"# Then rebuild the database index:\n"
            f"python OP_database_tool.py \\\n"
            f"    --input orthophyl_runs.tsv \\\n"
            f"    --output-dir {self.database_dir} \\\n"
            f"    --update"
        )
        commands.append(db_cmd)
        
        full_command = "\n\n".join(commands)
        
        decision = {
            'pipeline': 'OrthoPhyl',
            'reason': 'Novel taxonomy not represented in any database',
            'assembly': str(assembly_path),
            'assembly_id': assembly_id,
            'query_taxonomy': taxonomy,
            'download_taxonomy': tax_string_for_download,
            'download_rank': rank_name,
            'download_value': download_value,
            'suggestion': f"Create new database for {rank_name} '{download_value}'",
            'command': full_command
        }
        
        self._save_decision(decision)
        return decision
    
    def _save_decision(self, decision: Dict):
        """Save routing decision to files."""
        decision_file = self.output_dir / f"routing_decision_{decision['assembly_id']}.json"
        
        with open(decision_file, 'w') as f:
            json.dump(decision, f, indent=2)
        
        logger.info(f"\n→ Decision saved: {decision_file}")
        
        # Human-readable summary
        summary_file = self.output_dir / f"routing_summary_{decision['assembly_id']}.txt"
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("ASSEMBLY ROUTING DECISION (Multi-Database Query)\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Assembly ID: {decision['assembly_id']}\n")
            f.write(f"Assembly: {decision['assembly']}\n")
            f.write(f"Query Taxonomy: {decision['query_taxonomy']}\n")
            f.write(f"\nPipeline: {decision['pipeline']}\n")
            f.write(f"Reason: {decision['reason']}\n\n")
            
            if decision['pipeline'] == 'ReLeaf':
                f.write("Match Details:\n")
                f.write(f"  Database: {decision['matched_database']}\n")
                f.write(f"  Matched at: {decision['matched_rank']} level\n")
                f.write(f"  Database directory: {decision['database_dir']}\n")
                f.write(f"  Contains: {decision['database_genomes']} genomes\n")
                f.write(f"\nReLeaf Configuration:\n")
                f.write(f"  Tree method: {decision.get('tree_method', 'iqtree')}\n")
                f.write(f"  Data type: {decision.get('tree_data', 'CDS')}\n")
                if 'available_methods' in decision:
                    f.write(f"  Available methods: {', '.join(decision['available_methods'])}\n")
                if 'available_data_types' in decision:
                    f.write(f"  Available data types: {', '.join(decision['available_data_types'])}\n")
            elif decision['pipeline'] == 'OrthoPhyl_subclade_build':
                f.write("Match Details (unbuilt subclade):\n")
                f.write(f"  Subclade: {decision.get('subclade_name')}\n")
                f.write(f"  Parent taxon: {decision.get('parent_taxon')}\n")
                f.write(f"  Database directory: {decision.get('database_dir')}\n")
                f.write(f"  Members list: {decision.get('members_file')}\n")
                f.write(f"  Source genomes: {decision.get('source_genome_dir')}\n")
            else:
                f.write("Action Required:\n")
                f.write(f"  {decision['suggestion']}\n")
                f.write(f"  Download taxonomy: {decision['download_taxonomy']}\n")
                f.write(f"  Target rank: {decision['download_rank']}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("COMMAND(S) TO RUN\n")
            f.write("=" * 70 + "\n\n")
            f.write(decision['command'] + "\n")
        
        logger.info(f"→ Summary saved: {summary_file}")
    
    def batch_route(self, input_table: Path) -> List[Dict]:
        """Route multiple assemblies from table."""
        logger.info(f"\nBatch routing from {input_table}")
        
        decisions = []
        
        with open(input_table, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('#') and not first_line.lower().startswith('assembly'):
                f.seek(0)
            
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    logger.warning(f"Line {line_num}: insufficient fields")
                    continue
                
                assembly_path = Path(fields[0]).expanduser()  # Expand ~ to home directory
                taxonomy = fields[1].strip('"')  # Remove quotes if present
                assembly_id = fields[2] if len(fields) > 2 else None
                
                try:
                    decision = self.route_assembly(assembly_path, taxonomy, assembly_id)
                    decisions.append(decision)
                except Exception as e:
                    logger.error(f"Failed to route {assembly_path}: {e}")
        
        self._generate_batch_summary(decisions)
        return decisions
    
    def _generate_batch_summary(self, decisions: List[Dict]):
        """Generate batch summary."""
        summary_file = self.output_dir / "batch_routing_summary.txt"
        
        n_releaf = sum(1 for d in decisions if d['pipeline'] == 'ReLeaf')
        n_orthophyl = sum(1 for d in decisions if d['pipeline'] == 'OrthoPhyl')
        n_subclade_build = sum(1 for d in decisions if d['pipeline'] == 'OrthoPhyl_subclade_build')

        # Group by database
        by_database = {}
        for d in decisions:
            if d['pipeline'] == 'ReLeaf':
                db = d['matched_database']
                if db not in by_database:
                    by_database[db] = []
                by_database[db].append(d['assembly_id'])
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("BATCH ROUTING SUMMARY (Multi-Database)\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Available databases: {len(self.databases)}\n")
            for db in self.databases:
                f.write(f"  - {db['clade_name']} ({db['rank_name']}, {db['n_genomes']} genomes)\n")
            f.write("\n")
            
            f.write(f"Total assemblies: {len(decisions)}\n")
            f.write(f"  → ReLeaf (matched): {n_releaf}\n")
            f.write(f"  → OrthoPhyl (novel): {n_orthophyl}\n")
            f.write(f"  → OrthoPhyl subclade build (unbuilt subclade): {n_subclade_build}\n\n")

            if n_releaf > 0:
                f.write("ReLeaf Routing by Database:\n")
                f.write("-" * 70 + "\n")
                for db_name, assemblies in sorted(by_database.items()):
                    f.write(f"\n{db_name} ({len(assemblies)} assemblies):\n")
                    for asm in assemblies:
                        f.write(f"  - {asm}\n")
                f.write("\n")

            if n_orthophyl > 0:
                f.write("OrthoPhyl Assemblies (novel taxa):\n")
                f.write("-" * 70 + "\n")
                for d in decisions:
                    if d['pipeline'] == 'OrthoPhyl':
                        f.write(f"  {d['assembly_id']:30s} - {d['download_rank']}: {d['download_value']}\n")
                f.write("\n")

            if n_subclade_build > 0:
                f.write("OrthoPhyl Subclade Builds (unbuilt subclades to build then ReLeaf):\n")
                f.write("-" * 70 + "\n")
                for d in decisions:
                    if d['pipeline'] == 'OrthoPhyl_subclade_build':
                        f.write(f"  {d['assembly_id']:30s} - {d.get('subclade_name')} "
                                f"(parent {d.get('parent_taxon')})\n")
                f.write("\n")

        logger.info(f"\n→ Batch summary: {summary_file}")
        print(f"\nRouting complete: {n_releaf} → ReLeaf, {n_orthophyl} → OrthoPhyl, "
              f"{n_subclade_build} → OrthoPhyl subclade build")


def main():
    parser = argparse.ArgumentParser(
        description="Multi-database assembly router with automatic best-match selection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Route single assembly (queries all databases)
  python assembly_router.py \\
      --assembly genome.fna \\
      --taxonomy "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__VAXT01;s__" \\
      --database-dir /path/to/databases/ \\
      --output-dir results/
  
  # With gather_filter_asms.sh for downloading genomes
  python assembly_router.py \\
      --assembly genome.fna \\
      --taxonomy "d__Bacteria;..." \\
      --database-dir /path/to/databases/ \\
      --gather-filter-script utils/gather_filter_asms.sh \\
      --output-dir results/
  
  # Batch mode
  python assembly_router.py \\
      --batch assemblies.tsv \\
      --database-dir /path/to/databases/ \\
      --gather-filter-script utils/gather_filter_asms.sh \\
      --output-dir results/
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--assembly', help='Single assembly FASTA')
    input_group.add_argument('--batch', help='TSV: assembly_path, taxonomy, [id]')
    
    parser.add_argument('--taxonomy', help='GTDB taxonomy (for --assembly)')
    parser.add_argument('--assembly-id', help='Assembly identifier')
    parser.add_argument(
        '--database-dir',
        required=True,
        help='Directory containing multiple *_db databases'
    )
    parser.add_argument(
        '--output-dir',
        default='routing_output',
        help='Output directory'
    )
    parser.add_argument(
        '--gather-filter-script',
        help='Path to gather_filter_asms.sh for downloading genomes'
    )
    parser.add_argument('-t', '--threads', type=int, default=8, help='Threads')
    parser.add_argument(
        '--placement',
        choices=['subclade', 'backbone'],
        default='subclade',
        help=(
            "For megatree databases with tied matches: 'subclade' (default) places the "
            "query into the most MASH-similar dense subclade tree; 'backbone' places it "
            "into the sparse backbone tree for a broad overview."
        )
    )

    args = parser.parse_args()

    if args.assembly and not args.taxonomy:
        parser.error("--taxonomy required with --assembly")

    try:
        router = MultiDatabaseRouter(
            database_dir=Path(args.database_dir),
            output_dir=Path(args.output_dir),
            gather_filter_script=Path(args.gather_filter_script) if args.gather_filter_script else None,
            threads=args.threads,
            placement=args.placement
        )
    except Exception as e:
        logger.error(f"Failed to initialize: {e}")
        return 1
    
    try:
        if args.batch:
            decisions = router.batch_route(Path(args.batch))
        else:
            decision = router.route_assembly(
                Path(args.assembly).expanduser(),  # Expand ~ to home directory
                args.taxonomy.strip('"'),  # Remove quotes if present
                args.assembly_id
            )
            decisions = [decision]
        
        print("\n" + "=" * 70)
        print("ROUTING COMPLETE")
        print("=" * 70)
        
        for decision in decisions:
            print(f"\n{decision['assembly_id']}:")
            print(f"  Pipeline: {decision['pipeline']}")
            if decision['pipeline'] == 'ReLeaf':
                print(f"  Database: {decision['matched_database']} ({decision['matched_rank']})")
            elif decision['pipeline'] == 'OrthoPhyl_subclade_build':
                print(f"  Subclade: {decision.get('subclade_name')} "
                      f"(parent {decision.get('parent_taxon')}, build then ReLeaf)")
            else:
                print(f"  Action: {decision['suggestion']}")
            print(f"\n  Commands saved to: routing_summary_{decision['assembly_id']}.txt")
        
        return 0
        
    except Exception as e:
        logger.error(f"Routing failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())