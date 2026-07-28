"""Tree comparison utilities using ete3 for integration tests."""

from pathlib import Path
from typing import List, Set, Tuple

from ete3 import Tree


def clean_taxon_name(name: str) -> str:
    """
    Clean taxon name to match OrthoPhyl's CLEAN_N_COPY_GENOMES behavior.
    
    OrthoPhyl strips certain characters from genome filenames:
    - Parentheses: ( and )
    - Colons: :
    - Equals patterns: _=_
    
    This function applies the same transformations so test expectations
    match the actual processed names in OrthoPhyl outputs.
    
    Args:
        name: Raw taxon name (e.g., from filename stem).
        
    Returns:
        Cleaned taxon name matching OrthoPhyl's processing.
        
    Example:
        >>> clean_taxon_name("KU551270.1_Neottia_(PE)")
        'KU551270.1_Neottia_PE'
    """
    # Remove characters that OrthoPhyl strips from filenames
    cleaned = name.replace("(", "").replace(")", "")
    cleaned = cleaned.replace(":", "")
    cleaned = cleaned.replace("_=_", "")
    return cleaned


def load_tree(path: str) -> Tree:
    """
    Load a Newick tree file.
    
    Args:
        path: Path to the Newick tree file.
        
    Returns:
        ete3.Tree object.
        
    Raises:
        FileNotFoundError: If the tree file doesn't exist.
    """
    tree_path = Path(path)
    if not tree_path.exists():
        raise FileNotFoundError(f"Tree file not found: {path}")
    
    # format=1 is flexible Newick format (internal node names allowed)
    return Tree(str(tree_path), format=1)


def normalize_labels(tree: Tree) -> Tree:
    """
    Normalize leaf labels by stripping whitespace and quotes.
    
    This handles common variations in taxon naming (e.g., spaces, quotes)
    that don't affect biological identity but can break string comparisons.
    
    Args:
        tree: ete3.Tree object (modified in place).
        
    Returns:
        The same tree object (for chaining).
    """
    for leaf in tree.iter_leaves():
        # Strip whitespace and common quote characters
        leaf.name = leaf.name.strip().strip('"').strip("'")
    return tree


def rf_distance(tree1: Tree, tree2: Tree, unrooted: bool = True) -> Tuple[int, int, int, int, int]:
    """
    Compute Robinson-Foulds distance between two trees.
    
    Args:
        tree1: First tree.
        tree2: Second tree.
        unrooted: If True, compare as unrooted trees (default).
        
    Returns:
        Tuple of (rf, max_rf, common_leaves, tree1_leaves, tree2_leaves):
            - rf: Robinson-Foulds distance (number of differing bipartitions).
            - max_rf: Maximum possible RF distance for these trees.
            - common_leaves: Number of leaves in common.
            - tree1_leaves: Total leaves in tree1.
            - tree2_leaves: Total leaves in tree2.
    """
    result = tree1.robinson_foulds(tree2, unrooted_trees=unrooted)
    
    # ete3 versions differ: some return 5 values, newer ones return 6
    # Handle both cases for compatibility
    if len(result) == 6:
        # Newer ete3: (rf, max_rf, common, t1_leaves, t2_leaves, effective_tree_size)
        # We ignore the 6th value (effective_tree_size)
        return result[:5]
    else:
        # Older ete3: (rf, max_rf, common, t1_leaves, t2_leaves)
        return result


def same_topology(tree1: Tree, tree2: Tree, tolerance: int = 0) -> bool:
    """
    Check if two trees have the same topology within a tolerance.
    
    Args:
        tree1: First tree.
        tree2: Second tree.
        tolerance: Maximum allowed RF distance (default 0 = exact match).
        
    Returns:
        True if RF distance ≤ tolerance, False otherwise.
    """
    rf, max_rf, common, t1_leaves, t2_leaves = rf_distance(tree1, tree2)
    
    # If trees have different leaf sets, they can't have the same topology
    if t1_leaves != t2_leaves or common != t1_leaves:
        return False
    
    return rf <= tolerance


def get_leaf_names(tree: Tree) -> Set[str]:
    """
    Extract all leaf names from a tree.
    
    Args:
        tree: ete3.Tree object.
        
    Returns:
        Set of leaf names.
    """
    return {leaf.name for leaf in tree.iter_leaves()}


def assert_taxa_present(tree: Tree, expected_taxa: List[str]) -> None:
    """
    Assert that all expected taxa are present in the tree.
    
    Args:
        tree: ete3.Tree object.
        expected_taxa: List of taxon names that must be present.
        
    Raises:
        AssertionError: If any expected taxa are missing.
    """
    leaf_names = get_leaf_names(tree)
    missing = set(expected_taxa) - leaf_names
    
    if missing:
        raise AssertionError(
            f"Missing taxa in tree: {sorted(missing)}\n"
            f"Tree contains: {sorted(leaf_names)}"
        )


def assert_taxa_count(tree: Tree, expected_count: int) -> None:
    """
    Assert that the tree has exactly the expected number of taxa.
    
    Args:
        tree: ete3.Tree object.
        expected_count: Expected number of leaves.
        
    Raises:
        AssertionError: If the leaf count doesn't match.
    """
    leaf_count = len(list(tree.iter_leaves()))
    
    if leaf_count != expected_count:
        leaf_names = get_leaf_names(tree)
        raise AssertionError(
            f"Expected {expected_count} taxa, found {leaf_count}\n"
            f"Tree contains: {sorted(leaf_names)}"
        )
