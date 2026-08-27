# GBEK reference and figure reproduction

This directory contains a Python reference implementation of the parts of the
GBEK method needed to reproduce and validate results from:

> M. Gramsch, M. Balzer, M. Eckstein, and M. Kollar, *Phys. Rev. B* **88**,
> 235106 (2013).

The production C++ implementation is in `cincuenta/src`, principally
`NeqBathDecomposition.h`, `ImpuritySolverNeqGBEK.h`, and
`NeqDmftSolver.h`.

The Python and C++ implementations serve different roles:

- The Python code is an independent, small-system reference for the GBEK
  equations, bath decomposition, real-time propagation, and observables.
- `cincuenta` is the production C++ implementation.
- Comparison tools check that the production output agrees with the reference
  where the two runs represent the same physical setup.

The intended product is a reproducible set of paper-figure regenerations and
an accompanying LaTeX report built from those figures.

## Contents

### Reference implementation

- `gbek_ed.py` constructs the Fock-space basis and fermion operators.
- `gbek_dynamics.py` constructs the Hamiltonian, propagates states, and
  evaluates two-time Green functions.
- `gbek_cholesky.py` implements the causal rank-`L` Cholesky decomposition of
  the lesser hybridization.
- `gbek_selfconsistency.py` combines these pieces into the Bethe-lattice DMFT
  self-consistency calculation.  `run_self_consistency()` is its main entry
  point.
- `gbek_energy.py` evaluates kinetic, interaction, and total energies.
- `atomic_limit_reference.py` provides analytic atomic-limit checks.

The reference focuses on the atomic-limit start used for the relevant paper
calculations: no first bath, a cosine hopping ramp, Bethe-lattice
self-consistency, and spin-seed averaging.  It is not a general-purpose DMFT
solver.

### Figure and comparison tools

- `run_fig7_scan.py`, `run_fig8_scan.py`, `run_fig9_scan.py`, and
  `run_fig10_scan.py` generate the data for Figs. 7--10; `plot_docc_scan.py`
  and `plot_energy_scan.py` plot those data.
- `plot_fig3_errstep.py`, `plot_fig4_hybridization.py`, and `plot_errstep.py`
  generate the Fig. 3/4 hybridization diagnostics.
- `compare_reference.py`, `compare_V_rows.py`, `plot_cpp_vs_python_fig9.py`,
  and `plot_cpp_vs_python_fig10.py` compare Python-reference and C++ output.
- `cross_check_gbek_hamiltonian.py`, `cross_check_gbek_propagation.py`, and
  `cross_check_seed_scheme.py` derive the fixed values used by the C++ Catch2
  tests in `cincuenta/src/tests/test_ImpuritySolverNeqGBEK.cpp`.

The remaining `check_*`, `scan_*`, `investigate_*`, and `plot_*` scripts are
focused diagnostics.  Their module docstrings describe their inputs and
outputs.

## Tests

Run the Python reference smoke tests after changing its core modules:

```sh
./run_self_tests.sh
```

These checks cover the ED operators, propagation, Green functions, and the
causal Cholesky decomposition against analytic or independently constructed
results.  C++ production behavior is tested with the Catch2 suite under
`cincuenta/src/tests`.

## Regenerating figures and the report

Generated plots, data files, fetched paper figures, and PDFs are ignored by
Git.  Recreate them locally from the scripts and parameters in this directory:

```sh
./regenerate_plots.sh            # all figure groups
./regenerate_plots.sh --group-a  # Python reproductions of Figs. 7--10
./regenerate_plots.sh --group-b  # atomic-limit and C++/Python comparisons
./regenerate_plots.sh --report   # figures plus report.pdf and report-summary.pdf
```

Group B configures/builds `cincuenta` and runs its atomic-limit GBEK input
when required C++ output is missing or stale.  Set `DMRGPP_BUILD_DIR` to use
an existing build directory instead of this checkout's `build/` directory.

`report.tex` is the full LaTeX report and `report-summary.tex` is its concise
counterpart.  The report build requires `latexmk`, Ghostscript for EPS
conversion, and network access when `fetch_arxiv_figures.sh` needs to obtain
the paper figures.  It is not part of the CMake build.

## Direct reference run

For a small reference calculation:

```sh
uv run --with numpy --with scipy python3 gbek_selfconsistency.py --probe
uv run --with numpy --with scipy python3 gbek_selfconsistency.py \
  --L 3 --N 100 --dt 0.04 --U 2.0 --tq 0.25 \
  --out gbek-atomic-limit-exact-lesser
```

Use `compare_reference.py` or the plotting tools above to compare the result
with a matching `cincuenta` run.
