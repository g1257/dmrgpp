---
title: 'DMRG++: DMRG calculations across time, frequency, and temperature'
software_repository_url: 'https://github.com/dmrgpp-project/dmrgpp'
archive_doi: '10.11578/dc.20260601.1'
tags:
  - C++
  - condensed matter
  - physics
authors:
  - name: Gonzalo Alvarez
    corresponding: true # (This is how to denote the corresponding author)
    orcid: 0000-0002-1496-8261
#    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Daniel Arndt
    orcid: 0000-0001-8773-4901
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Damien Lebrun-Grandie
    orcid: 0000-0003-1952-7219
    affiliation: "1"
  - name: Peter Doak
    orcid: 0000-0001-6039-9752
    affiliation: 1
  - name: Steven Hahn
    orcid: 0000-0002-2018-7904
    affiliation: 1
  - name: Thomas A. Maier
    orcid: 0000-0002-1424-9996
    affiliation: 1
  - name: Andrey Prokopenko
    orcid: 0000-0003-3616-5504
    affiliation: 1
affiliations:
 - name: Oak Ridge National Laboratory
   index: 1
   ror: 01qz5mb56
date: 13 August 2026
bibliography: dmrgpp.bib

# Optional fields heres
---

# Summary

DMRG++ is an open-source C++ application for simulations of strongly correlated
quantum systems using the density matrix renormalization group (DMRG)
[@White:1992; @White:1993]. The DMRG
algorithm is used to study many-body Hamiltonians in condensed-matter physics,
quantum chemistry, and Hamiltonian lattice gauge theory. DMRG++ focuses on
condensed-matter models and provides an input-driven framework for selecting the
physical model, geometry, interactions, symmetries, and numerical method.

DMRG++ computes ground states and correlation functions as well as
frequency-dependent, time-dependent, and finite-temperature quantities. It
supports built-in geometry kinds together with input-defined connections and
coupling values. Researchers can therefore study different strongly correlated
systems without writing a separate DMRG implementation for each Hamiltonian.

# Scientific background

Many problems in condensed-matter physics, quantum chemistry, and Hamiltonian
lattice gauge theory can be expressed as an eigenvalue problem,

$$H|\psi\rangle=E|\psi\rangle,$$

where the Hamiltonian $H$ describes interacting local degrees of freedom. These
degrees of freedom may be lattice sites, electronic orbitals, or lattice-gauge
links and sites. If each of $L$ local degrees of freedom has $d$ states, the
dimension of the full Hilbert space is

$$\dim \mathcal{H}=d^L.$$

Fixing particle number or other conserved quantities reduces this space, but
its dimension still grows combinatorially. Directly storing an
arbitrary state or diagonalizing the full Hamiltonian therefore becomes
impractical after few sites.

Numerical methods avoid this cost by sampling or restricting the many-body
space. Auxiliary-field quantum Monte Carlo, for example, can have polynomial
cost when its sign or phase problem is absent, but polynomial scaling does not hold
generally [@Motta:2018]. DMRG instead uses the restricted entanglement of many
low-energy states. Ground states of local, gapped one-dimensional Hamiltonians
satisfy an entanglement area law and can be approximated efficiently by matrix
product states [@Hastings:2007; @Schollwock:2011].

In the conventional formulation used by DMRG++, the lattice is divided into
system and environment blocks. At each step, the reduced density matrix ranks
the block states by their contribution to the target state. DMRG retains the
states with the largest density-matrix eigenvalues and discards the rest, so an
exponentially large Hilbert space need not be represented explicitly
[@White:1992]. Repeated sweeps across the lattice optimize the retained basis.
The sum of the discarded eigenvalues and convergence with the number of retained
states provide practical diagnostics of truncation error. In the equivalent
modern description, DMRG variationally optimizes a matrix product state (MPS),
and the number of retained states is its bond dimension [@Schollwock:2011].

DMRG is most effective for chains, where a moderate bond dimension often gives
high accuracy and the cost at fixed bond dimension is polynomial in the chain
length. Ladders can be mapped to a one-dimensional ordering, but the entanglement
and required bond dimension generally grow exponentially with the ladder width,
that is, its short or transverse dimension [@Stoudenmire:2012]. DMRG++ supplies
the model, geometry, symmetry, measurement, and solver infrastructure needed to
apply this procedure to a range of interacting Hamiltonians.

# Statement of need

Research applications of DMRG often require changes to the Hamiltonian, lattice
connectivity, conserved quantities, target states, and observables. Implementing
these changes in a problem-specific code can duplicate work on basis
construction, sweeps, measurements, and validation. Researchers studying
strongly correlated lattice Hamiltonians therefore need reusable solver
infrastructure that can accommodate different physical models and calculation
protocols without requiring a new DMRG implementation for each study.

DMRG++ [@Alvarez:2009] addresses this need as an input-driven C++ application implementing
conventional DMRG [@White:1992]. Users select from model
implementations compiled into the application and specify interaction-dependent
connections and coupling matrices in the input. Representative models include
one-, two-, and three-band Hubbard variants, Kondo models, Heisenberg and Kitaev
spin models, and Holstein and Hubbard-Holstein electron-phonon models. Built-in
geometry kinds and general connection matrices or coupling values supplied
through the input share the same solver machinery, and calculations can use the
conserved quantum-number sectors supported by each model.

The shared infrastructure supports ground-state and static-correlation
calculations, time evolution, frequency-dependent response, and
finite-temperature observables [@Alvarez:2011; @Alvarez:2013;
@Nocera:2016a; @Nocera:2016b]. Researchers can therefore change a Hamiltonian
or calculation protocol while retaining the same simulation and analysis
workflow.

# State of the field

The broad applicability of DMRG has produced a diverse software ecosystem. A
recent survey identified more than 50 implementations and compared 37 packages
in detail [@Sehlstedt:2026]. These span conventional DMRG and MPS or
tensor-network formulations, application-oriented programs and general-purpose
libraries, and several scientific domains. Their differing approaches to
Hamiltonian construction, symmetries, targeting, and numerical backends reflect
distinct research workflows rather than a single hierarchy of implementations.

Independent packages nevertheless reproduce substantial core functionality.
Sehlstedt et al. found overlap in symmetry handling and high-performance
strategies, while most packages provide their own Lanczos- or Davidson-based
eigensolver [@Sehlstedt:2026]. Such specialization serves individual domains but
increases the work needed to maintain numerical infrastructure and adopt new
backends. Shared numerical components, standard interfaces, and collaboration
offer complementary ways to reduce that duplication.

Within this landscape, DMRG++ is a free and open-source C++ implementation of
conventional DMRG. Development began at Oak Ridge National Laboratory in 2008
to solve condensed-matter problems with systematically monitored truncation
error; the first software paper appeared in 2009 [@Alvarez:2009]. DMRG++ is an
application-oriented code configured through input files, not a general-purpose
tensor-network library. It uses generic programming and a comparatively small
software dependency set [@Sehlstedt:2026]. Accelerator-oriented batched kernels
were originally developed in the separate `DmrgppPluginSc` repository and were
integrated into the main DMRG++ repository in 2026.^[The integration is recorded
by repository commit
[`1b818160`](https://github.com/dmrgpp-project/dmrgpp/commit/1b8181600faa7829accf836de28180bee866d1a7),
“Bring in files from DmrgppPluginSc,” dated April 7, 2026.]

# Software design

DMRG++ is an input-driven C++17 scientific application that combines generic
programming with runtime configuration. The principal `dmrg` executable reads
and validates a single input file, after which `DmrgRunner` assembles the types
and objects needed for the calculation. This includes choosing real or complex
arithmetic, a Hamiltonian matrix-vector implementation, the geometry, the
physical model, and the targeting method. The calculation has four main layers:
geometry, model, DMRG engine, and targeting and measurement.

The model and geometry layers separate the local physics from the connections
between sites. A model defines its local Hilbert space, operators, Hamiltonian
terms, parameters, and conserved quantum numbers. The geometry defines the site
ordering, connections, and coupling values for each interaction term. Thus, the
same model can use different built-in geometry kinds or general connection
matrices and coupling values supplied through the input, while the geometry
machinery is shared by different models. Geometry implementations themselves
remain compiled components. Model implementations derive from a common
interface and are selected by `ModelSelector`. New models are added as C++
classes, registered with the selector, and compiled into DMRG++.

`DmrgSolver` implements the model-independent growth and finite-sweep stages. It
constructs system and environment blocks, obtains target states through Lanczos
diagonalization, truncates each block, and transforms states and operators into
the retained basis. Bases carry both operators and quantum numbers, allowing the
calculation to be divided into symmetry sectors and represented with
block-sparse matrices. Lanczos accesses the Hamiltonian through a common
matrix-vector interface. DMRG++ can store the sparse Hamiltonian, apply it
on-the-fly through the model, or evaluate its Kronecker-product decomposition
without constructing the full superblock matrix. The conventional numerical
path uses threaded execution, BLAS, and LAPACK. Through PsimagLite, Kokkos and
Kokkos Kernels will provide portable kernels in the future. For now, an optional integrated batched path can
use MAGMA as an accelerator backend.

The targeting layer determines which states contribute to the reduced density
matrix during a sweep. Implementations of the `TargetingBase` interface reuse
the same model, geometry, basis, and sweep code for ground-state calculations,
frequency-dependent response, correction-vector methods, real-time evolution,
Chebyshev expansions, ancilla-based finite-temperature evolution, minimally
entangled typical thermal states (METTS) [@Stoudenmire:2010], and
expression-based targeting. The
ancilla path reuses the time-step targeting infrastructure for purification,
while METTS has a dedicated targeting implementation. Measurements follow a
similar separation. Selected operator expressions can be evaluated in situ
during the main run, while the separate `observe` executable reads saved
simulation data to calculate one-point and multipoint correlations afterward.

DMRG++ writes model and solver parameters, the geometry, an encoded copy of the
input, energies, target states, renormalized bases, and checkpoint data to HDF5.
These records support restart and recovery as well as post-processing without
repeating the full simulation. Program and source-revision information is also
reported with each run. CMake organizes common implementation code into internal
library targets linked by the command-line executables; these `make` targets avoid code
duplication but are not presented as a public programming interface. The `dmrg`
executable runs simulations, `observe` post-processes saved results, `manyOmegas`
coordinates calculations at multiple frequencies, and `introspect` examines
model bases, operators, and Hamiltonian terms. Regression tests exercise
representative input files and compare selected energies and observables with
expected results. The model, targeting, and matrix-vector interfaces provide
compiled extension points while retaining the same solver, input,
and test infrastructure.

# Research impact statement

The project received research support during two distinct periods between 2011
and 2020, and a new period of research funding began in 2025. [FIXME: Add funding
agencies, programs, award numbers, other funding details, and corrected exact
dates.] Continued maintenance has enabled its models, algorithms, numerical
kernels, and analysis tools to evolve.

Over time, the scientific scope of DMRG++ expanded beyond ground states to
nonequilibrium dynamics, frequency-dependent response, and finite-temperature
observables.

As of 2026, the seven authors of this paper maintain DMRG++ as part-time
developers. The authors know of approximately 15 active users and collaborators;
this figure is not intended as a complete census of everyone who has run the
software. More than 80 publications cite the original DMRG++ implementation
paper [@Alvarez:2009].^[Citation records for DOI
`10.1016/j.cpc.2009.02.016` were checked in Crossref and OpenAlex on July 16,
2026. Crossref reported 85 citations, while an OpenAlex citing-works query
returned 94 records, or 92 after deduplication by DOI. Because database coverage
and record handling differ, we report the conservative lower bound “more than
80 publications.”] This aggregate count provides a conservative measure of the
software's visibility in the research literature.

DMRG++ is free and open-source software under the UT-Battelle open-source
license, and its development takes place publicly at
<https://github.com/dmrgpp-project/dmrgpp>. The repository exposes the commit
history, issues, pull requests, contributor activity, and continuous-integration
workflows. Users can therefore inspect changes, report defects, propose new
features, and contribute code through the same public development process.

# AI usage disclosure

Hermes Agent (Nous Research), using the GPT-5.6-sol model, assisted with
planning, drafting, revising, and checking this manuscript. The authors reviewed
and edited all AI-assisted material and take responsibility for the accuracy and
final content of the paper. AI assistance was also used in selected DMRG++
commits, which identify that assistance in their commit messages.

# Acknowledgements

# References
