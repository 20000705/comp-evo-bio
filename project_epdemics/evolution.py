# evolution.py
from __future__ import annotations
from typing import List, Optional, Tuple
import numpy as np
from numpy.random import Generator, default_rng
from trees import Tree


def branch_lengths_time(tree: Tree) -> List[float]:
    """
    Compute time-based branch lengths: dt = time(child) - time(parent).
    Returns a list L where L[i] is the branch length (parent[i] -> i).
    Root length is set to 0.0.
    """
    n = tree.n_nodes()
    L = [0.0] * n
    for i, p in enumerate(tree.parent):
        if p == -1:
            L[i] = 0.0
        else:
            dt = tree.time[i] - tree.time[p]
            if dt < 0:
                raise ValueError("Tree times must be nondecreasing along edges (parent_time <= child_time).")
            L[i] = float(dt)
    return L


def branch_lengths_mutations(
    tree: Tree,
    mu: float,
    rng: Optional[Generator] = None,
) -> List[float]:
    """
    Simulate mutation-count branch lengths using a molecular clock:
      k ~ Poisson(mu * dt) on each edge.

    Returns a list M where M[i] is the mutation count on (parent[i] -> i).
    Root length is set to 0.0.
    """
    rng = default_rng() if rng is None else rng

    dt_lengths = branch_lengths_time(tree)
    n = tree.n_nodes()
    M = [0.0] * n
    for i in range(n):
        if tree.parent[i] == -1:
            M[i] = 0.0
        else:
            lam = mu * dt_lengths[i]
            M[i] = float(rng.poisson(lam))
    return M


def export_newicks_time_and_mut(
    tree: Tree,
    mu: float,
    rng: Optional[Generator] = None,
    precision: int = 6,
) -> Tuple[str, str]:
    """
    Convenience: return (newick_time, newick_mutations).
    """
    L_time = branch_lengths_time(tree)
    L_mut  = branch_lengths_mutations(tree, mu=mu, rng=rng)

    newick_time = tree.to_newick(branch_lengths=L_time, precision=precision)
    newick_mut  = tree.to_newick(branch_lengths=L_mut,  precision=precision)

    return newick_time, newick_mut