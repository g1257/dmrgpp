# GBEK atomic-limit reference generator

Standalone reference implementation (numpy/scipy only, no dependency on
cincuenta, Ainur, or LanczosPlusPlus) that computes the "exact"
-i*Λ^+_<(t,t') hybridization for the atomic-limit hopping quench setup of

    Gramsch, Balzer, Eckstein, Kollar, PRB 88, 235106 (2013), Sec. VI.

This is the quantity plotted in the top-left panel of the paper's Fig. 3: a
self-consistently converged Weiss field, obtained by exact diagonalization of
a small SIAM, used as ground truth against which cincuenta's production
rank-L Cholesky second-bath approximation (`NeqBathDecomposition.h`,
`ImpuritySolverNeqGBEK.h`, elsewhere in the `cincuenta` C++ tree) is
validated.

## Why this exists

Getting an "exact" self-consistent Weiss field equivalent to the one used in
the paper proved impossible without significant modification to cincuenta:
producing it there means going through the full Ainur/LanczosPlusPlus
machinery, which has several incidental limitations at the atomic limit
(nBath=0, impurity seeded with a single spin) --

1. the Ainur vector grammar cannot parse an empty `[]` Connectors literal,
2. a hardcoded half-filling constraint (`nup+ndown==nsites`), and
3. Pauli-forbidden zero-dimensional Fock sectors when seeding the impurity
   with a single spin (no valid particle-added/removed sector exists).

Constraints 1 and 3 were resolved on the C++ side once this Python
reference existed to guide and validate a closed-form atomic-limit bypass
(`ImpuritySolverNeqExactDiag::solveAtomicLimit`, nBath=0 only). Constraint 2
remains a general limitation of cincuenta's equilibrium/near-equilibrium
exact-diagonalization solver (removing it is a tracked follow-up, unrelated
to the atomic-limit setup specifically).

So this directory reimplements the atomic-limit Weiss-field calculation
directly and independently of the production C++ code, to produce a
trustworthy comparison target and, over time, ended up growing Python
implementations of the rest of the GBEK pipeline (Cholesky decomposition,
propagation, self-consistency) as well, so cincuenta's C++ and this
independent Python implementation could be cross-checked against each
other. In practice bugs were found and fixed in both directions, not just
Python informing C++.

## What this is

**This is not a general-purpose GBEK solver**, although many of the paper's
results can be reproduced with just this Python code. For lack of a better
name, this directory is referred to as the "py-GBEK reference
implementation." It only handles the one setup used throughout: atomic-limit
start (Λ^- = 0 identically, no first bath), a cosine hopping ramp,
Bethe-lattice self-consistency, and the alpha/beta spin-seed averaging the
paper uses to restore particle-hole symmetry.

## Physics implemented

- **Atomic limit**: no first bath. `Λ = Λ^+` entirely.
- **Hopping (bandwidth) quench**: `v(t)` cosine ramp 0->1 over `[0, tq]`,
  constant afterward. `U` fixed throughout (this is a hopping quench, not an
  interaction quench).
- **Bethe-lattice self-consistency**: `Λ(t,t') = -i * hop(t) G(t,t') hop(t')`,
  `hop(t) = t_star_f * v(t)`.
- **Second bath**: `L` pairs of bath sites (one initially doubly-occupied,
  one initially empty), coupled via a rank-`L` causal Cholesky decomposition
  of `Λ^<(t,t')` -- see `gbek_cholesky.py`.
- **G-sigma averaging**: `G_sigma = 1/2 (G_alpha,0sigma + G_beta,0sigma)`.
  By the up<->down / alpha<->beta symmetry of this setup,
  `G_beta,0sigma = G_alpha,0,-sigma`, so only the alpha (impurity seeded
  spin-up) trajectory needs to be propagated; averaging its up- and
  down-spin correlators gives the same result as averaging the alpha and
  beta seeds, at half the cost.

## Files

### Core library

- `gbek_ed.py` -- Fock-space basis construction and fermion operators
  (bit-manipulation based, independent ED implementation).
- `gbek_dynamics.py` -- Hamiltonian construction, real-time propagation
  (hand-rolled truncated-Taylor-series propagator, with an automatic
  fallback to `scipy.sparse.linalg.expm_multiply` if a step doesn't
  converge), and two-time Green's function extraction via forward/backward
  propagation.
- `gbek_cholesky.py` -- rank-`L` causal Cholesky decomposition
  (`cholesky_causal`), implemented directly from the paper's Eq. 56-63.
  Also keeps two deliberately-incorrect variants,
  `cholesky_causal_buggy_fixed_window()` and `cholesky_causal_linear_only()`,
  purely as regression fixtures the Catch2 tests compare against.
- `gbek_energy.py` -- kinetic/interaction/total energy observables
  (paper Fig. 9/10).
- `gbek_selfconsistency.py` -- wires the above into the DMFT
  self-consistency loop; `run_self_consistency()` is the main entry point,
  `dump_lesser()` writes the two-time `t t' Re Im` file format shared with
  cincuenta's own `KadanoffBaym::dump`.
- `atomic_limit_reference.py` -- closed-form analytic atomic-limit Green's
  functions, cross-checked against a direct 4-state ED.
- `gbek_colormap.py` -- diverging colormap sampled from the paper's own
  Fig. 3 colorbar, for visual consistency with the published figures.
- `provenance.py` -- stamps generated plots/data with the git commit and a
  hash of any dirty tracked files, so a result can be traced back to the
  exact code that produced it.

### Paper-figure reproduction (pure Python, no cincuenta needed)

`run_fig7_scan.py` / `run_fig8_scan.py` / `run_fig9_scan.py` /
`run_fig10_scan.py` reproduce the self-consistency runs behind paper Figs.
7-10 (double occupation vs. DMFT iteration and `L_bath`; energy
conservation; double occupation vs. `U`), writing `.npz` data consumed by
`plot_docc_scan.py` (Figs. 7/8) and `plot_energy_scan.py` (Figs. 9/10).
`plot_fig4_hybridization.py`, `plot_errstep.py`, and `plot_fig3_errstep.py`
reproduce Fig. 4 and Fig. 3's bottom-left `err^step(t)` panel directly from
a `gbek_selfconsistency.py::dump_lesser` output.

### Comparison against cincuenta's C++ output

- `compare_reference.py` -- overlays this Python reference against
  cincuenta's own rank-L Cholesky output (`ImpuritySolverNeqGBEK::dumpPlusBath`).
- `compare_neq_delta_lesser.py` (`cincuenta/TestSuite/`, one directory up)
  -- plots a single two-time `t t' Re Im` file as a 2D colormap in the
  paper's Fig. 3 style; works equally on this directory's own
  `dump_lesser()` output or a cincuenta C++ dump, since both share the
  format.
- `plot_cpp_vs_python_fig9.py` / `plot_cpp_vs_python_fig10.py` -- overlay
  cincuenta's C++ `docc`/energy observable dumps
  (`ImpuritySolverNeqGBEK::dumpDoccAndEnergy`) against this directory's own
  Fig. 9/10 reproduction.
- `compare_V_rows.py` -- row-by-row comparison of the Cholesky factor
  `V(n,p)` between cincuenta's online C++ computation and an offline batch
  replay of `gbek_cholesky.py::cholesky_causal()` on the same target, to
  localize exactly where two implementations diverge, if they do.

### C++ test-fixture derivation

`cross_check_gbek_hamiltonian.py`, `cross_check_gbek_propagation.py`, and
`cross_check_seed_scheme.py` are independent, from-scratch derivations used
to generate the hardcoded expected values in
`cincuenta/src/tests/test_ImpuritySolverNeqGBEK.cpp`'s Catch2 tests. If the
Hamiltonian construction, propagation, or seeding scheme is ever
intentionally changed, rerun the relevant script and update that test's
expected values from its printed output -- these scripts are the permanent
source of truth for those fixtures, not just one-off checks.
`check_alpha_beta_docc_symmetry.py` verifies the `Gbeta_0sigma =
Galpha_0,-sigma` simplifying assumption the double-occupation figures
(Figs. 7/8) rely on; `cross_check_beta_trajectory.py` performs the
equivalent check for the underlying Green's function itself.

### One-off diagnostics

`check_fixedpoint_residual.py`, `scan_t3_activation.py`, `scan_t3_t4_t5.py`,
`plot_errstep_t3scan.py`, `investigate_L4_tmax4.py`, `check_cholesky_step.py`,
and `quantify_lambda_minus_leak.py` were written to answer specific
questions during past investigations (documented in each script's own
docstring) and are kept for provenance and as reusable templates for
similar future investigations -- they are not a regression suite (see
"Self-tests" below) and most print a comparison for a human to read rather
than asserting pass/fail.

### Report build

`report.tex` / `report-summary.tex`, `fetch_arxiv_figures.sh`, and
`regenerate_plots.sh` -- see "Building the LaTeX report" below.

## Validation

Each core-library piece is checked against an independent exact/analytic
result (see "Self-tests" below for how to run these):

- ED operators: atomic-limit ground-state energy = `-U/4` per site (matches
  the paper's stated atomic-limit energy).
- Propagation + two-time `G^<`: matches an analytic single-particle
  (`U=0`, static hopping) oscillation and Green's function to machine
  precision.
- Cholesky decomposition: reproduces the C++ Catch2 tests' expected
  behavior on exact-rank targets (exact reconstruction at `L = rank`,
  strictly worse at `L < rank`, row 0 identically zero), and matches a
  from-the-paper's-equations derivation on a genuinely full-rank target.
- Full self-consistency loop (`L=3`, `N=100`, `dt=0.04`, `tmax=4`, the
  paper-matching grid): converges to `max|dLambda| < 1e-8` in ~15
  iterations. The diagonal `Λ(t,t)` is flat at its particle-hole-symmetric
  exact value `0.5*hop(t)^2` across the entire time range -- a genuine
  invariant of the physics, independent of Cholesky rank `L`. The full
  matrix is positive-semidefinite to numerical precision.

## Self-tests

`gbek_ed.py`, `gbek_dynamics.py`, `gbek_cholesky.py`, and
`atomic_limit_reference.py` each carry a smoke test in their own
`if __name__ == "__main__":` block, validating against an independent
analytic/exact result (see "Validation" above for what each one checks).
Nothing runs these automatically -- `./run_self_tests.sh` is the "run
everything, fail loudly on the first problem" entry point; run it after
changing any of these modules, or before trusting a fresh checkout. This is
a smoke-test harness, not a full regression suite: the one-off diagnostic
scripts documented above are not included in it.

## Usage

    cd cincuenta/TestSuite/gbek_reference
    uv run --with numpy --with scipy python3 gbek_selfconsistency.py --probe
    uv run --with numpy --with scipy python3 gbek_selfconsistency.py \
        --L 3 --N 100 --dt 0.04 --U 2.0 --tq 0.25 --out gbek-atomic-limit-exact-lesser

    # compare against cincuenta's own rank-L Cholesky approximation:
    uv run --with numpy --with scipy --with matplotlib python3 compare_reference.py \
        gbek-atomic-limit-exact-lesser /path/to/gebk-fig3-L3-plus-bath-lesser --tmax 4.0

    # check the Cholesky update step itself against a run's own dumped data:
    uv run --with numpy --with scipy python3 check_cholesky_step.py \
        /path/to/gebk-fig3-L3-weiss-delta-lesser /path/to/gebk-fig3-L3-plus-bath-lesser \
        --V <fitted V_alpha, comma-separated> --eps <fitted eps_alpha, comma-separated> \
        --beta <FicticiousBeta> --L 3

    # or, plot the exact reference alone with the existing tooling:
    python3 ../compare_neq_delta_lesser.py gbek-atomic-limit-exact-lesser \
        --tmax 4.0 --title "GEBK Fig. 3: exact -i Λ^<_+ (atomic limit reference)"

## Regenerating plots

**No `.png`, `.npz`, or `.provenance.txt` file in this directory is
committed to the repo** (see `.gitignore` here) -- binary blobs don't
diff/review usefully and bit-rot in git history; the source (scripts +
parameters) is what's checked in. Regenerate anything you need locally:

    ./regenerate_plots.sh            # everything (groups A+B)
    ./regenerate_plots.sh --group-a  # just Fig. 7/8 (pure Python, fast)
    ./regenerate_plots.sh --group-b  # the older C++-dependent validation
                                      # plots (builds+runs cincuenta if the
                                      # prerequisite dumps aren't already
                                      # in build/)
    ./regenerate_plots.sh --report   # build report.pdf (see below)

Group A (`fig7_docc.png`, `fig8_docc.png`) is pure Python:
`run_fig7_scan.py --L 2`, `run_fig7_scan.py --L 4`, `run_fig8_scan.py`,
`investigate_L4_tmax4.py`, then `plot_docc_scan.py --figure 7 --L 2,4` /
`--figure 8`.

!!!! Nothing below here will work until the follow cincuenta GBEK PR
lands !!!!

Group B depends on actual `cincuenta` C++ dumps existing in `build/`
first: `inputNeqAtomicLimitGBEKL3.ain` produces the
`build/atomic-limit-gbek-L3-*` files used by `plot_atomic_limit_2d.py`,
`plot_collapse_evidence_summary.py`, `plot_errstep_t3scan.py`,
`scan_t3_activation.py`, and (as the `approx` argument) `compare_reference.py`;
`inputNeqGBEKFig3L3.ain` produces the `build/gebk-fig3-L3-*` files used by
`plot_fig3l3_post_fix.py`, `plot_collapse_evidence_summary.py`, and (as
the `prefix` argument) `quantify_lambda_minus_leak.py`. `compare_reference.py`
also needs the pure-Python exact reference (`gbek_selfconsistency.py --L 3
--N 100 --dt 0.04 --U 2.0 --tq 0.25 --out gbek-atomic-limit-exact-lesser`,
also in `regenerate_plots.sh`). See each script's own header for its exact
expected inputs.

## Building the LaTeX report

`report.tex` is the canonical, reproducible write-up of this effort (the
content previously maintained by hand as a Claude Artifact). It embeds
groups A+B's plots plus the paper's own Figs. 3, 4, 7, 8, 9, 10, fetched
directly from the arXiv e-print source (`fetch_arxiv_figures.sh`, arXiv
1306.6315) rather than checked-in crops -- neither the fetched paper
figures nor the built PDF are committed (see `.gitignore`), same
reasoning as the plots above.

    ./regenerate_plots.sh --group-a --group-b   # or already-present plots
    ./regenerate_plots.sh --report

Requires network access (arXiv fetch), a TeX install with `latexmk`, and
`ghostscript` (`brew install ghostscript`) for the EPS-to-PDF conversion
step. Not part of the cmake build.

`report-summary.tex` builds `report-summary.pdf` from the same source
with the development-history callouts omitted -- a status/capabilities-only
copy for sharing with people who aren't interested in development-process
detail (see `report.tex`'s `\PIVERSION` header comment). Both PDFs build
together via `--report` above, or individually with `latexmk -pdf
report.tex` / `latexmk -pdf report-summary.tex`.
