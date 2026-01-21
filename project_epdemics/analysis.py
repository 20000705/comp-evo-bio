from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional
import numpy as np
import matplotlib.pyplot as plt

from sim import ODEResult, GillespieResult
from trees import Tree


def plot_ode(result: ODEResult, labels: List[str], title: str = "") -> None:
    t = result.t
    y = result.y  # (dim, len(t))

    if len(labels) != y.shape[0]:
        raise ValueError(f"labels length {len(labels)} must match y.shape[0]={y.shape[0]}")

    for i, lab in enumerate(labels):
        plt.plot(t, y[i], label=lab)

    plt.xlabel("time")
    plt.ylabel("count (or fraction)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_gillespie(result: GillespieResult, labels: List[str], title: str = "") -> None:
    t = result.times
    y = result.states  # (n_steps, dim)

    if len(labels) != y.shape[1]:
        raise ValueError(f"labels length {len(labels)} must match y.shape[1]={y.shape[1]}")

    for i, lab in enumerate(labels):
        plt.step(t, y[:, i], where="post", label=lab)

    plt.xlabel("time")
    plt.ylabel("count")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()


@dataclass
class TreeStats:
    n_nodes: int
    n_leaves: int
    sackin: int
    colless: Optional[int]  # None if not strictly binary


def summarize_tree(tree: Tree) -> TreeStats:
    n_nodes = tree.n_nodes()
    n_leaves = len(tree.leaves())
    sackin = tree.sackin_index()

    try:
        colless = tree.colless_index()
    except Exception:
        colless = None

    return TreeStats(n_nodes=n_nodes, n_leaves=n_leaves, sackin=sackin, colless=colless)


def compare_distributions(values_a: np.ndarray, values_b: np.ndarray, name: str = "metric") -> None:
    """
    Lightweight comparison: show mean/std and visualize distribution overlap.

    Notes
    -----
    - For publication, you may later add statistical tests (KS, MWU) and effect sizes.
    - Uses shared histogram bins so plots are comparable.
    """
    values_a = np.asarray(values_a, dtype=float)
    values_b = np.asarray(values_b, dtype=float)

    mean_a = float(values_a.mean()) if len(values_a) else float("nan")
    mean_b = float(values_b.mean()) if len(values_b) else float("nan")
    std_a = float(values_a.std(ddof=1)) if len(values_a) > 1 else 0.0
    std_b = float(values_b.std(ddof=1)) if len(values_b) > 1 else 0.0

    print(f"{name} A: mean={mean_a:.4f}, std={std_a:.4f}, n={len(values_a)}")
    print(f"{name} B: mean={mean_b:.4f}, std={std_b:.4f}, n={len(values_b)}")

    # Shared bins for fair comparison
    pooled = np.concatenate([values_a, values_b]) if (len(values_a) and len(values_b)) else (values_a if len(values_a) else values_b)
    bins = np.histogram_bin_edges(pooled, bins=20) if len(pooled) else 20

    plt.hist(values_a, bins=bins, alpha=0.5, label="A")
    plt.hist(values_b, bins=bins, alpha=0.5, label="B")
    plt.xlabel(name)
    plt.ylabel("count")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Quick robust summary
    plt.boxplot([values_a, values_b], labels=["A", "B"])
    plt.ylabel(name)
    plt.title(f"{name}: A vs B")
    plt.tight_layout()
    plt.show()