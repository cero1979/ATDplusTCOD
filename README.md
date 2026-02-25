# Canonical Poset–TDA Pipeline for Binary Linear Codes  
## Syndrome Cayley Graphs, Canonical Ranks, and Rank–Truncated Order Complexes

This repository accompanies an experimental implementation of a **canonical Topological Data Analysis (TDA) pipeline** for **binary linear codes**, organized around the combinatorics of **syndrome Cayley graphs** and the induced **graded poset** structure on \(\mathbb{F}_2^m\). The core artifact is the notebook:

- `ExpATDTC.ipynb` — *Experimental Notebook: Syndrome Cayley Graph and Rank-Truncated Order Complexes*

The notebook follows the pipeline described in **Ramirez Ovalle (2025)** and uses examples motivated by the study of binary linear codes via **coset-leader Hasse diagrams** (Delgado Castillo, 202x).

---

## 1. Mathematical setting

Let \(C \subseteq \mathbb{F}_2^n\) be a binary linear code with parameters \([n,k]\), and let  
\[
H \in \mathbb{F}_2^{m\times n}, \qquad m=n-k,
\]
be a parity-check matrix for \(C\). The **syndrome space** is \(\mathbb{F}_2^m\), with \(|\mathbb{F}_2^m|=2^m\).

Write the columns of \(H\) as \(\{h_i\}_{i=1}^n \subseteq \mathbb{F}_2^m\) (not necessarily distinct; the implementation uses the set of distinct columns as generators).

---

## 2. Pipeline overview (canonical construction)

### Step 1 — Syndrome Cayley graph
Define the (undirected) **syndrome Cayley graph** \(G_H\) by:
- **Vertices:** \(V(G_H)=\mathbb{F}_2^m\)
- **Edges:** \(s \sim s+h_i\) for each generator \(h_i\) (addition mod 2)

This graph encodes the action of adding columns of \(H\) on syndromes.

### Step 2 — Canonical rank function via BFS
Define the **canonical rank** (a.k.a. canonical distance-to-zero) by:
\[
\mathrm{rk}(s) := d_{G_H}(0,s),
\]
the shortest-path distance in \(G_H\) from the zero syndrome to \(s\).

In the coding-theoretic interpretation (as stated in the notebook), \(\mathrm{rk}(s)\) coincides with the minimum Hamming weight among error vectors \(e \in \mathbb{F}_2^n\) satisfying \(He^T=s\), i.e. the minimum weight in the corresponding coset.

**Algorithmic realization:** BFS on \(G_H\) from \(0\).  
**Complexity:** \(O(n\cdot 2^m)\) time, \(O(2^m)\) space.

### Step 3 — Canonical cover relation (graded poset)
Using the rank function, define a **graded poset** structure on \(\mathbb{F}_2^m\) by the cover relation:
\[
u \prec v \iff \exists\, h_i:\ v=u+h_i \ \text{and}\ \mathrm{rk}(v)=\mathrm{rk}(u)+1.
\]
This selects precisely those Cayley edges that increase rank by one, yielding a directed acyclic graph which serves as a **Hasse diagram** for the canonical poset.

### Step 4 — Rank truncation and order-complex filtration
For \(t \in \mathbb{Z}_{\ge 0}\), define the rank-truncated subposet:
\[
P_{\le t} := \{ s \in \mathbb{F}_2^m \;:\; \mathrm{rk}(s)\le t \}.
\]
Let \(\Delta(P_{\le t})\) denote its **order complex**, whose simplices correspond to finite chains
\[
s_0 < s_1 < \cdots < s_\ell \quad \text{in } P_{\le t}.
\]
This yields a filtration:
\[
K_t := \Delta(P_{\le t}), \qquad t=0,1,2,\dots
\]
A chain is inserted at filtration value \(\max_j \mathrm{rk}(s_j)\).

### Step 5 — Unpunctured vs punctured variants
The notebook computes persistence for two related filtrations:

- **Unpunctured filtration:** includes the zero syndrome \(0\). In many examples, \(0\) acts as an apex producing a cone-like structure, and the complex is expected to be contractible (hence “trivial” persistence in positive degrees).

- **Punctured filtration:** removes the apex \(0\) (denoted “punctured=True” in the implementation). This can destroy the cone structure and allow nontrivial homology to appear.

### Step 6 — Persistent homology
Persistent homology is computed from the filtration using **GUDHI** (`gudhi.SimplexTree`), producing persistence pairs, diagrams, and barcodes.

---

## 3. What the notebook implements

### 3.1 Core functions
The notebook defines seven reusable functions forming the computational backbone:

1. `bfs_canonical_ranks(H)`  
   BFS on \(G_H\) \(\Rightarrow\) canonical ranks \(\mathrm{rk}(\cdot)\).

2. `build_cover_relation(H, ranks)`  
   Canonical cover pairs \(u\prec v\) with rank increase \(+1\).

3. `build_hasse_graph(ranks, covers)`  
   Directed graph representation of the Hasse diagram (for visualization).

4. `build_cayley_graph(H)`  
   Explicit Cayley graph on \(\mathbb{F}_2^m\) (primarily for visualization).

5. `build_order_complex(ranks, covers, max_rank=None, punctured=False)`  
   Builds the filtered order complex as a `gudhi.SimplexTree`. Supports:
   - truncation via `max_rank`
   - puncturing the zero syndrome via `punctured=True`

6. `find_coset_leaders(H)`  
   Brute-force enumeration of all \(e\in\mathbb{F}_2^n\) to compute coset leaders and minimum weights (verification only).  
   **Complexity:** \(O(2^n)\) (intended for small \(n\)).

7. `rank_distribution_table(ranks)`  
   Summary statistics: number of syndromes at each rank; covering radius estimate \(\rho=\max \mathrm{rk}\).

### 3.2 Visualization and reporting
Helper routines include:
- Cayley-graph plots with nodes colored by rank
- Layered Hasse-diagram plots (by rank)
- Persistence diagrams and barcodes
- Tables of persistence pairs and Betti numbers along the filtration

---

## 4. Example study: Hamming \([7,4,3]\)

The notebook includes a worked example for the Hamming code \([7,4,3]\), used as a verification benchmark.

Key facts (as recorded in the notebook):
- \(m=3\), hence \(2^m=8\) syndromes.
- The columns of \(H\) exhaust all nonzero vectors of \(\mathbb{F}_2^3\).
- Expected covering radius: \(\rho=1\).
- Rank distribution: one element at rank \(0\), seven at rank \(1\).
- The Hasse diagram is a “flat star” of depth \(1\).
- The punctured complex at \(t=1\) yields 7 isolated vertices, so \(\beta_0=7\) at that level.

The notebook also checks agreement between:
- BFS-derived ranks \(\mathrm{rk}(s)\), and
- brute-force minimum coset weights (for small \(n\))

---

## 5. Dependencies

The notebook uses:

- `numpy`
- `networkx`
- `matplotlib`
- `pandas`
- `gudhi` (required for persistent homology)

Install (example):

```bash
pip install numpy networkx matplotlib pandas gudhi
```

If `gudhi` is not installed, the notebook sets a flag (`HAS_GUDHI=False`) and skips persistence computations while printing a warning.

---

## 6. Remarks on complexity / scalability

- The BFS rank computation scales as \(O(n\cdot 2^m)\) and is feasible for moderate \(m\).
- The order-complex construction can become combinatorially large because it implicitly enumerates chains in the poset.
- The brute-force coset-leader verification scales as \(O(2^n)\) and is intended only for small blocklength \(n\).

---

## 7. References (as cited in the notebook)

- Ramirez Ovalle (2025): canonical Poset–TDA pipeline.
- Delgado Castillo (202x): coset-leader Hasse diagrams and related constructions in coding theory.

---

## License

Add your license here (e.g., MIT/BSD/GPL) if applicable.
