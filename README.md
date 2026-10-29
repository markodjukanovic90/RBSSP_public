# RBSSP_public

This repository contains the source code, benchmark instances, and raw computational results accompanying the manuscript:

> **Optimization Techniques for Solving the Minimum Red--Blue Separating Set Problem with an Application to Human--AI Text Separation**  
> **Marko Đukanović**, **Staša Vujičić Stanković**, **Žiga Zebec**, and **Christian Blum**

The repository is intended to support reproducibility of the computational study and to provide implementations of the optimization methods considered for the **Minimum Red--Blue Separating Set Problem (MRBSSP)**.

## Problem overview

Let \(G=(V,E)\) be an undirected graph with two disjoint sets of distinguished vertices:

- \(R \subseteq V\): red vertices,
- \(B \subseteq V\): blue vertices.

For a vertex \(v\in V\), let \(N[v]\) denote its closed neighborhood. A vertex \(v\) separates a red--blue pair \((r,b)\in R\times B\) if \(N[v]\) contains exactly one of \(r\) and \(b\).

For each red--blue pair, define

\[
D(r,b)=\{v\in V : |N[v]\cap\{r,b\}|=1\}.
\]

A set \(S\subseteq V\) is a feasible red--blue separating set if

\[
S\cap D(r,b)\neq\emptyset
\qquad
\forall (r,b)\in R\times B.
\]

The objective of MRBSSP is to find a feasible separating set of minimum cardinality.

Equivalently, MRBSSP can be viewed as a structured Set Cover / Hitting Set problem in which red--blue pairs are the elements to be covered and graph vertices are the candidate covering objects.

## Methods implemented

The repository contains implementations used in the computational study of several exact and heuristic approaches.

### Classical ILP

The direct integer linear programming formulation introduces one binary variable \(x_v\) for each vertex and one covering constraint for each red--blue pair:

\[
\min \sum_{v\in V} x_v
\]

subject to

\[
\sum_{v\in D(r,b)} x_v \ge 1,
\qquad
\forall (r,b)\in R\times B.
\]

The implementation also uses preprocessing/kernelization and a heuristic warm start.

### Logic-Based Benders Decomposition (LBBD)

The LBBD approach starts from a relaxed master problem containing the decision variables, objective, essential-vertex fixings, and a valid upper-bound cutoff obtained from GRASP.

The red--blue covering constraints are generated dynamically. After each master solve, a logical feasibility oracle identifies uncovered red--blue pairs and adds the corresponding violated constraints.

A configurable parameter \(q_{\mathrm{LBBD}}\) limits the number of violated cuts added per LBBD iteration.

### Branch-and-Benders-Cut (BBC)

BBC embeds the same logical separation mechanism directly into a single branch-and-bound search tree through lazy constraints.

Whenever an integer candidate is found, the callback checks whether all red--blue pairs are separated. Violated separation constraints are then added lazily without restarting the search tree.

The parameter \(q_{\mathrm{BBC}}\) controls the maximum number of lazy cuts generated in one callback invocation.

### GRASP

A Greedy Randomized Adaptive Search Procedure (GRASP) is used both:

- as a standalone heuristic, and
- as a source of feasible warm starts and upper bounds for the exact approaches.

The construction phase uses a randomized restricted candidate list, followed by a local-search phase based on vertex removal while preserving feasibility.

### Variable Neighborhood Search (VNS)

The computational study also includes a VNS implementation/baseline for comparison with GRASP and the exact methods.

## Preprocessing / kernelization

The exact approaches use common preprocessing rules before optimization.

The main reductions include:

1. **Essential vertices**  
   If \(|D(r,b)|=1\), the unique vertex in \(D(r,b)\) must belong to every feasible solution.

2. **Row dominance**  
   If \(D_i \subset D_j\), then the row corresponding to \(D_j\) is dominated and can be removed.

3. **Column dominance**  
   If the set of red--blue pairs covered by a vertex \(u\) is a strict subset of the set covered by another vertex \(v\), then \(u\) can be removed.

These reductions are applied iteratively while keeping the surviving separation sets synchronized with the surviving vertex set.

## Repository structure

```text
RBSSP_public/
├── instances/
├── raw_results/
└── src/
```

### `instances/`

Contains the benchmark instances used in the computational experiments.

The manuscript considers two main benchmark families:

- **Random graph instances**, used for systematic algorithmic evaluation and parameter tuning;
- **NLP-derived instances**, constructed for the application to human--AI text separation.

### `raw_results/`

Contains raw outputs from the computational experiments reported in the manuscript.

These results are intended to support verification and reproduction of the reported comparisons, including solution quality, lower bounds, optimality gaps, optimality certificates, and method-specific statistics such as LBBD iterations and BBC callback activity.

### `src/`

Contains the source code implementing the optimization and heuristic approaches used in the study, including the common graph-processing and preprocessing routines.

## Human--AI text separation application

The second benchmark family is derived from an NLP application in which human-written and AI-generated texts are represented through graph-based vocabulary information.

In this setting, red and blue vertices correspond to the two text classes, while candidate separator vertices represent vocabulary-based features. Solving MRBSSP therefore identifies a minimum set of features sufficient to distinguish every human--AI pair represented in the graph.

This application provides an interpretable feature-selection perspective in addition to the purely combinatorial benchmark study.

## Experimental focus

The computational experiments compare the methods using several complementary criteria, including:

- feasible solution quality / upper bounds;
- lower-bound strength;
- optimality gap;
- number of proven optimal instances;
- running time;
- LBBD master re-solves;
- BBC callback/separation effort;
- sensitivity to \(q_{\mathrm{LBBD}}\) and \(q_{\mathrm{BBC}}\);
- performance on both random and NLP-derived instances.

The repository includes the raw computational material corresponding to these experiments.

## Software requirements

The code is written in Python and the exact approaches rely on IBM ILOG CPLEX through DOcplex / the CPLEX Python API.

Typical dependencies used by the implementations include:

- Python 3;
- `networkx`;
- `docplex`;
- IBM ILOG CPLEX;
- `pulp` for selected auxiliary functionality.

A working CPLEX installation and a valid CPLEX Python environment are required for the exact ILP, LBBD, and BBC implementations.

## Reproducibility

To reproduce the computational study:

1. clone this repository;
2. install the required Python packages and configure IBM ILOG CPLEX;
3. select an instance from `instances/`;
4. run the desired implementation from `src/`;
5. compare the produced output with the corresponding files in `raw_results/`.

Because solver paths, cluster configurations, and time limits may depend on the local computing environment, users may need to adjust these settings before running the exact methods.

## Manuscript

This repository accompanies the manuscript:

> **Optimization Techniques for Solving the Minimum Red--Blue Separating Set Problem with an Application to Human--AI Text Separation**

Authors:

- **Marko Đukanović** — University of Nova Gorica / IZUM  
  Corresponding author: `marko.dukanovic@ung.si`
- **Staša Vujičić Stanković** — `stasa.vujicic.stankovic@matf.bg.ac.rs`
- **Žiga Zebec** — IZUM  
  `ziga.zebec@izum.si`
- **Christian Blum** — Artificial Intelligence Research Institute (IIIA-CSIC)  
  `christian.blum@iiia.csic.es`

If you use the code, instances, or experimental results from this repository, please cite the accompanying manuscript.

## Citation

A BibTeX entry can be completed once the final bibliographic information for the manuscript is available.

```bibtex
@article{DjukanovicEtAl_RBSSP,
  title   = {Optimization Techniques for Solving the Minimum Red--Blue Separating Set Problem with an Application to Human--AI Text Separation},
  author  = {Marko Dukanovic and Stasa Vujicic Stankovic and Ziga Zebec and Christian Blum},
  note    = {Manuscript},
  year    = {2026}
}
```

## Contact

For questions concerning the repository or the computational experiments, please contact:

**Marko Đukanović**  
`marko.dukanovic@ung.si`

## License

Please add the intended software/data license to the repository and update this section accordingly.
