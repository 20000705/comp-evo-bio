from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class Tree:
    """
    Minimal rooted tree representation (rooted, directed away from root).

    Representation
    --------------
    Nodes are indexed 0..n-1.

    parent[i] : int
        Parent node index of node i, or -1 if i is the root.

    time[i] : float
        Node time (e.g., infection time / coalescent time / sampling time).
        Assumption: time is nondecreasing along root -> leaf
        (i.e., parent time <= child time).
        NOTE: The current metrics use topology only. Time is useful later
        for branch lengths (child_time - parent_time).

    label[i] : Optional[str]
        Optional node labels (often important for leaves in phylogenetics).

    Notes
    -----
    - This structure is useful for both transmission trees and evolutionary trees.
      The difference is interpretation of edges (transmission vs ancestry) and
      what "time" represents (infection time vs divergence time).
    """
    parent: List[int]
    time: List[float]
    label: Optional[List[str]] = None
    validate: bool = True  # turn off for speed if needed

    def __post_init__(self):
        n = len(self.parent)

        if len(self.time) != n:
            raise ValueError("parent and time must have the same length.")
        if self.label is not None and len(self.label) != n:
            raise ValueError("label must have the same length as parent/time.")

        if self.validate:
            self._validate_tree()

    # ----------------- Basic structure helpers -----------------

    def n_nodes(self) -> int:
        return len(self.parent)

    def root(self) -> int:
        """Return the unique root node index."""
        roots = [i for i, p in enumerate(self.parent) if p == -1]
        if len(roots) != 1:
            raise ValueError(f"Expected exactly 1 root, found {len(roots)}.")
        return roots[0]

    def children_lists(self) -> List[List[int]]:
        """Adjacency list of children for each node."""
        n = self.n_nodes()
        ch = [[] for _ in range(n)]
        for i, p in enumerate(self.parent):
            if p != -1:
                ch[p].append(i)
        return ch

    def leaves(self) -> List[int]:
        """Nodes with no children."""
        ch = self.children_lists()
        return [i for i in range(self.n_nodes()) if len(ch[i]) == 0]
    
    def to_newick(
        self,
        branch_lengths: Optional[List[float]] = None,
        precision: int = 6,
    ) -> str:
        """
        Export to Newick format.

        Parameters
        ----------
        branch_lengths : Optional[List[float]]
            If provided, branch_lengths[i] is the length of the branch (parent[i] -> i).
            Convention: branch_lengths[root] is ignored (can be 0.0).
            If None, exports topology-only Newick (no ':length').

        precision : int
            Number of decimal places for branch lengths.

        Notes
        -----
        This keeps Tree generic: you can export time-based trees or mutation-based trees
        depending on what branch_lengths you pass in.
        """
        ch = self.children_lists()
        r = self.root()
        n = self.n_nodes()

        if branch_lengths is not None and len(branch_lengths) != n:
            raise ValueError("branch_lengths must have the same length as number of nodes")

        def node_name(u: int) -> str:
            return self.label[u] if self.label is not None else f"n{u}"

        def bl(u: int) -> str:
            if branch_lengths is None:
                return ""
            if self.parent[u] == -1:
                return ""  # root has no incoming branch
            return f":{branch_lengths[u]:.{precision}f}"

        def rec(u: int) -> str:
            kids = ch[u]
            name = node_name(u)
            if not kids:
                return f"{name}{bl(u)}"
            inside = ",".join(rec(v) for v in kids)
            return f"({inside}){name}{bl(u)}"

        return rec(r) + ";"


    # ----------------- Validation -----------------

    def _validate_tree(self) -> None:
        """
        Validate that:
        - exactly one root exists
        - no cycles exist
        - all nodes are reachable from the root (connected)
        - optional time monotonicity holds: time[parent] <= time[child]
        """
        r = self.root()  # ensures exactly one root

        ch = self.children_lists()
        n = self.n_nodes()

        # DFS to check reachability and cycles
        state = [0] * n  # 0=unvisited, 1=visiting, 2=done

        def dfs(u: int):
            state[u] = 1
            for v in ch[u]:
                if state[v] == 1:
                    raise ValueError("Cycle detected in parent pointers.")
                if state[v] == 0:
                    dfs(v)
            state[u] = 2

        dfs(r)

        if any(s == 0 for s in state):
            raise ValueError("Tree is disconnected: some nodes are not reachable from the root.")

        # Check time monotonicity along edges
        for child, par in enumerate(self.parent):
            if par != -1 and self.time[par] > self.time[child]:
                raise ValueError("Time must be nondecreasing along root->leaf (parent_time <= child_time).")

    # ----------------- Newick export -----------------

    def to_newick(self, include_branch_lengths: bool = False, precision: int = 6) -> str:
        """
        Export to Newick format.

        include_branch_lengths:
            If True, add branch lengths using (time[child] - time[parent]).
            This is useful for evolutionary trees (genetic time / divergence time).
        """
        ch = self.children_lists()
        r = self.root()

        def node_name(u: int) -> str:
            return self.label[u] if self.label is not None else f"n{u}"

        def branch_len(u: int) -> str:
            p = self.parent[u]
            if p == -1:
                return ""  # root has no parent branch
            bl = self.time[u] - self.time[p]
            return f":{bl:.{precision}f}"

        def rec(u: int) -> str:
            kids = ch[u]
            name = node_name(u)

            if not kids:
                # Leaf
                if include_branch_lengths:
                    return f"{name}{branch_len(u)}"
                return name

            inside = ",".join(rec(v) for v in kids)
            if include_branch_lengths:
                return f"({inside}){name}{branch_len(u)}"
            return f"({inside}){name}"

        return rec(r) + ";"

    # ----------------- Topology metrics -----------------

    def depths(self) -> List[int]:
        """Compute topological depth (#edges from root) for each node."""
        n = self.n_nodes()
        d = [-1] * n
        r = self.root()
        d[r] = 0

        for i in range(n):
            if d[i] != -1:
                continue

            # Walk upward until hitting a node whose depth is known
            path = []
            u = i
            while u != -1 and d[u] == -1:
                path.append(u)
                u = self.parent[u]

            base = d[u] if u != -1 else 0  # should not happen in connected tree
            for k, node in enumerate(reversed(path), start=1):
                d[node] = base + k

        return d

    def sackin_index(self) -> int:
        """Sackin index: sum of leaf depths (topology only)."""
        d = self.depths()
        return sum(d[i] for i in self.leaves())

    def colless_index(self) -> int:
        """
        Colless index for strictly binary trees.

        For each internal node, add |nL - nR| where nL,nR are the number of leaves
        in the left and right subtrees.
        """
        ch = self.children_lists()

        # Ensure strict binary for internal nodes
        for u in range(self.n_nodes()):
            if len(ch[u]) not in (0, 2):
                raise ValueError("Colless index requires a strictly binary tree.")

        def leaf_count(u: int) -> int:
            kids = ch[u]
            if not kids:
                return 1
            return leaf_count(kids[0]) + leaf_count(kids[1])

        def rec(u: int) -> int:
            kids = ch[u]
            if not kids:
                return 0
            L, R = kids
            nL = leaf_count(L)
            nR = leaf_count(R)
            return abs(nL - nR) + rec(L) + rec(R)

        return rec(self.root())