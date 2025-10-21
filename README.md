# JC69 Tree Likelihood Calculator

## Overview

Given:

* a phylogenetic tree in **Newick** format with branch lengths, and
* a DNA sequence alignment in **FASTA** or **PHYLIP** format,

the program computes the **total log-likelihood** of the tree under **JC69**.

Core steps:

1. Parse the Newick tree into a binary tree structure (multifurcations are handled by left-heavy binarization).
2. Read and validate the alignment (FASTA or PHYLIP).
3. Match sequence labels to tree tip labels (1:1 check).
4. Traverse the tree in post-order and compute conditional likelihoods at each node (Felsenstein pruning).
5. Combine root partials with base frequencies to compute the total log-likelihood.

---

## How to Run

### Basic Commands

**PHYLIP:**

```bash
python main.py --tree tree.nwk --data data.phy --format phylip
```

**FASTA:**

```bash
python main.py --tree tree.nwk --data data.fa --format fasta
```

### Useful Flags

* `--log`  : Output the log-likelihood (default behavior prints a single number).
* `--check` : Print diagnostics (number of tips, sites, and tip order) before computation.
* `--debug` : Verbose mode; prints intermediate values during pruning.

### Example

```bash
python main.py \
  --tree tests/case1.nwk \
  --data tests/case1.phy \
  --format phylip \
  --log --check
```

Example output:

```
TIPS=6 SITES=12 labels_ok=True
Tip order: ['A', 'B', 'C', 'D', 'E', 'F']
-30.022675
```

---

## Numerical Stability Choice

This implementation uses **per-node scaling** during pruning to ensure numerical stability, even for large alignments (e.g., ≥ 10,000 sites).

* Each node’s partials are divided by their sum, and the log of the scaling factor is accumulated per site.
* At the root, the total log-likelihood includes both the root site log-likelihood and the accumulated scaling factors.
* This removes floating-point underflow without log-space recursion.

---

## Tests

Automated tests live under `tests/`:

* **Case 1 — Larger tree (PHYLIP):**

  * Validates pruning on 6 tips × 12 sites.
  * Checks diagnostics, finite negative logL, and clean exit.
* **Case 2 — Missing data (FASTA):**

  * Includes `N` and gap-only sequences; verifies correct handling and valid likelihood.
* **Case 3 — Label mismatch (FASTA):**

  * Extra/missing labels; expects non-zero exit with informative error.

Run individually:

```bash
python tests/case1-test.py
python tests/case2-test-missing.py
python tests/case3-test-badlabel.py
```

---

## Known Limitations

* Assumes uniform base frequencies; no empirical frequencies.
* Only JC69 model is implemented; HKY/GTR not included unless extended.
* Tree is treated as bifurcating during traversal (multifurcations are binarized left-heavy).

---

## Possible Extensions

* **HKY/GTR** models with `--model` flag and empirical base frequencies from the alignment.
* **Discrete-Γ rate variation** (k=4 categories, shape `α`).
* **log-space** pruning for arbitrarily large alignments.

---

## 📚 References

* Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution*, 17(6), 368–376.
* Yang, Z. (2014). *Molecular Evolution: A Statistical Approach*. Oxford University Press.
* Course lecture materials on likelihood calculation.
