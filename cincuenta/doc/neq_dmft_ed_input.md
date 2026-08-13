# Non-equilibrium DMFT input, exact-diagonalization (ED) solver

This guide is for someone running `cincuenta` to do a non-equilibrium DMFT
(neq-DMFT) interaction-quench calculation, using the exact-diagonalization
(ED) impurity solver, **without** needing to read the C++. It assumes you can
already run an equilibrium `cincuenta` job (see
[`cincuenta_design.md`](cincuenta_design.md) and `cincuenta/README.md`).

If you only need to get a run going quickly, skip to
[Worked example](#worked-example) and come back to the parameter tables as
needed.

## The physical picture

`cincuenta` first solves the **equilibrium** DMFT problem at U = U₀ to get a
converged bath (the `{V_α, ε_α}` in `cincuenta_design.md`'s Anderson model).
Then, in a second stage, it quenches the interaction from `U_i` to `U_f` at
t = 0 and follows the resulting real-time dynamics on the Keldysh contour,
**with the bath held fixed at its equilibrium value** — this is the "fixed
bath" approximation described in the paper this input format is modeled on,
Aoki, Tsuji et al., PRB **88**, 235106 (2013). (A local copy of this paper
is available to developers via the repo's `CLAUDE.md` overlay; it is not
part of the git-tracked repo itself.)

Non-equilibrium mode is switched on simply by the presence of `TmaxNeq=` in
the input file with a positive value. Everything above `TmaxNeq=` in your
`.ain` file is the ordinary equilibrium input (same keywords you already
know: `FicticiousBeta=`, `NumberOfBathPoints=`, `ImpuritySolver=`, etc.).
Everything below is neq-specific.

There are three neq **solvers**, chosen by `NeqSolver=`:

| `NeqSolver=` | Solver class | What it does |
|---|---|---|
| *(omitted)* | `ImpuritySolverNeqExactDiag` | **This guide's subject.** Exact diagonalization on the fixed equilibrium bath. |
| `"tdmrg"` | `ImpuritySolverNeqTdmrg` | Time-dependent DMRG instead of ED. Not covered here — ask a developer, or see the class itself. |
| `"gbek"` | `ImpuritySolverNeqGBEK` | Adds a second, low-rank ("Cholesky") bath on top of ED, per Gramsch/Balzer/Eckstein/Kollar, PRB **88**, 235106 (2013). |

This guide covers the default (no `NeqSolver=` line) **ED** path. The GBEK
scheme reuses several of the same keywords (`NeqBathRank=`, `NeqAtomicLimit=`)
but has its own subtleties documented in
`cincuenta/TestSuite/gbek_reference/README.md` — read that if you set
`NeqSolver="gbek"`.

## Required equilibrium input (recap)

The neq stage is bolted onto an ordinary equilibrium run, so all the usual
required equilibrium keywords still apply and are still required:
`FicticiousBeta=`, `ChemicalPotential=`, `Matsubaras=`, `LatticeGf=`,
`NumberOfBathPoints=`, `DmftNumberOfIterations=`, `DmftTolerance=`,
`ImpuritySolver=` (must be `"exactdiag"` for the neq-ED path to make sense —
the equilibrium bath fit is what the neq stage's bath comes from).

Two keywords used by the ED solver itself, read directly (not through
`ParamsDmftSolver`/`ParamsNeqDmftSolver`), so they're easy to miss if you're
scanning the `Params*.h` files instead of the actual `.ain` examples:

| Keyword | Type | Meaning |
|---|---|---|
| `TargetElectronsUp=` | int | Number of up-spin electrons in the target sector |
| `TargetElectronsDown=` | int | Number of down-spin electrons in the target sector |

**Half-filling constraint:** the ED neq solver hardcodes the impurity
chemical potential to `mu_imp = -U/2` (particle-hole symmetry). It will
`err()` at run time unless `TargetElectronsUp + TargetElectronsDown ==
NumberOfBathPoints + 1` (i.e., exactly half-filling of the full star-geometry
system: 1 impurity + `NumberOfBathPoints` bath sites). There is no way to run
away from half filling with this solver today.

The equilibrium run also needs its usual DMRG-omega-scan block
(`RootOutputname=`, `FiniteLoopsGs=`, `OmegaBegin=`, etc.) even when
`ImpuritySolver="exactdiag"` — `ParamsDmftSolver`'s constructor and the
downstream `DmftSolver` machinery require those fields to parse regardless of
which solver ultimately consumes them. Copy them from one of the working
examples below rather than trying to omit them.

## Neq-specific input

All of the following go below the equilibrium block, and are parsed by
`ParamsNeqDmftSolver` (`cincuenta/src/ParamsNeqDmftSolver.h`) unless noted.

### Required

| Keyword | Type | Meaning |
|---|---|---|
| `HubbardU=` | real | U_i, the **initial** (equilibrium) interaction. Same keyword the equilibrium block already needs — not duplicated, just reused. |
| `HubbardUFinal=` | real | U_f, the interaction after the t = 0 quench. Set equal to `HubbardU=` for a no-quench sanity check. |
| `TmaxNeq=` | real | Total real-time window to simulate, t ∈ [0, TmaxNeq]. **Presence of this line (with a positive value) is what turns on neq mode at all.** |
| `NtNeq=` | int | Number of real-time steps. `dt = TmaxNeq / NtNeq` is fixed automatically — you don't set `dt` directly. |

### Optional (all default sensibly if omitted)

| Keyword | Type | Default | Meaning |
|---|---|---|---|
| `NeqDmftIter=` | int | same as equilibrium `DmftNumberOfIterations=` | Inner DMFT self-consistency iterations *at each time step*. The loop itself always runs (`NeqDmftSolver::timeStep()` is shared by every neq solver), but each iteration calls `prepareTimeStep()`, which is a no-op for plain ED (the bath is fixed) — so on the plain-ED path, extra iterations are harmless but do nothing; only `NeqSolver="gbek"` (which updates its Cholesky bath in `prepareTimeStep()`) actually needs more than 1. |
| `NeqDmftTolerance=` | real | same as equilibrium `DmftTolerance=` | Convergence tolerance for the above. |
| `NeqSolver=` | string | *(empty → plain ED)* | `"tdmrg"` or `"gbek"`; see solver table above. This label is declared in the Ainur grammar but **read directly in `cincuenta.cpp`**, not by `ParamsNeqDmftSolver` — if you go looking for it in `ParamsNeqDmftSolver.h` you won't find it. |
| `NeqBathRank=` | int | 0 | GBEK second-bath rank L. Irrelevant unless `NeqSolver="gbek"`. |
| `BandwidthFinal=` | real | 0 (no quench) | Bethe-lattice bandwidth *for t > 0*, i.e. a **hopping** quench on top of (or instead of) the interaction quench. 0 means the lattice bandwidth stays at its equilibrium `LatticeGf=` value for all t. |
| `QuenchShape=` | string | `"step"` | Ramp shape for the hopping quench: `"step"`, `"cosine"`, or `"tanh"`. Only matters if `BandwidthFinal=` is nonzero. |
| `QuenchDuration=` | real | 0 (instantaneous) | Ramp duration t_q for the hopping quench. |
| `NeqOutputPrefix=` | string | `""` | Prefix prepended to output Green's-function filenames (see [Output files](#output-files)). |
| `NeqAtomicLimit=` | int (0/1) | 0 | See [NeqAtomicLimit gotcha](#neqatomiclimit-only-affects-the-gbek-path) below — **read this before using it.** |

## Gotchas

### `NeqAtomicLimit=` only affects the GBEK path

`NeqAtomicLimit=1` is meant to start the neq run from the true atomic limit
(no bath at t = 0, so the hybridization Λ⁻(t,t′) is identically zero — the
setup in GBEK Sec. VI). Looking at how it actually flows through
`cincuenta.cpp`: the flag is read into `neqParams.neqAtomicLimit`, and an
"empty bath" substitution (`neqBathParams`) is built from it — **but that
substitution is only ever used on the `NeqSolver="gbek"` branch.** Both the
default (bare ED) branch and the `NeqSolver="tdmrg"` branch pass the
equilibrium-fitted bath (`dmftSolver.bathResult()`) unconditionally,
ignoring `NeqAtomicLimit=` entirely.

In practice this means:
- With `NeqSolver="gbek"` and `NeqAtomicLimit=1`: works as intended — see
  `inputNeqAtomicLimitGBEKL3.ain` for a full worked example (also below).
- With no `NeqSolver=` (plain ED) and `NeqAtomicLimit=1`: **the flag does
  nothing.** `ImpuritySolverNeqExactDiag::solve()` dispatches to
  `solveAtomicLimit()` (the closed-form single-Hubbard-atom path) purely
  based on `bathParams.size() == 0` — and on the plain-ED branch,
  `cincuenta.cpp` always passes the equilibrium-fitted bath
  (`dmftSolver.bathResult()`), never the empty one, regardless of
  `NeqAtomicLimit=`. Setting `NumberOfBathPoints=0` will **not** work around
  this either: the equilibrium bath-fit minimizer requires at least one bath
  parameter pair and errors out before the neq stage is reached (confirmed by
  the comment in `inputNeqAtomicLimitGBEKL3.ain`, which needs
  `NumberOfBathPoints=1` for exactly this reason even though it discards the
  eq fit result). **There is currently no supported way to reach the plain-ED
  atomic limit from the input file.** If you need it, use
  `NeqSolver="gbek"` with `NeqAtomicLimit=1` (which does work, and itself
  bottoms out in the same `solveAtomicLimit()`), or ask a developer to wire
  `NeqAtomicLimit=` into the default branch in `cincuenta.cpp`.

This is existing behavior, not something this doc changes — flagged here so
you don't spend an afternoon debugging why `NeqAtomicLimit=1` had no visible
effect on a plain ED run. It's also a real input-ergonomics gap worth fixing
in the code (a declared flag that's silently inert on two of three solver
branches), separate from this documentation pass.

### `HubbardU=` is shared between the two stages

There's only one `HubbardU=` line in the file. It's U_i for the neq quench
*and* the equilibrium bath-fit's interaction — you don't write it twice.
`HubbardUFinal=` is the only new interaction keyword.

### DMRG-only fields are still required even for `ImpuritySolver="exactdiag"`

See [Required equilibrium input](#required-equilibrium-input-recap) above.
This trips people up because it looks redundant, but omitting it will fail
input parsing before you even get to the neq stage.

## Output files

After the neq run finishes, `dumpGreenFunctions()` writes (via
`KadanoffBaym::dump`, with `<prefix>` = `NeqOutputPrefix=` value or empty):

| File | Contents |
|---|---|
| `<prefix->green-retarded` | Retarded impurity Green's function G^R(t, t′) |
| `<prefix->green-lesser` | Lesser impurity Green's function G^<(t, t′) |
| `<prefix->green-matsubara-t` | Mixed Matsubara/real-time component |
| `<prefix->weiss-green-retarded/lesser/matsubara-t` | Same three components for the Weiss field G₀ (the bath-only, U = 0 reference) |
| `<prefix->weiss-delta-...` | Hybridization Λ(t, t′) components. (The filenames say "delta" — that's `KadanoffBaym::dump`'s literal, unrenamed output-file naming; the physical quantity is Λ, not Δ.) |

Two boundary-condition identities are checked automatically in CI
(`cincuenta/src/CMakeLists.txt`, `neqTsuji0` test): `G^R(0,0) = -i` and
`G^M(τ=0) = -n_imp`. If your own run's outputs don't satisfy these, something
upstream is wrong before you get to physics interpretation.

For the U_i = U_f = 0 sanity check (Σ = 0 exactly ⟹ G_imp must equal the
Weiss field G₀), use `cincuenta/TestSuite/compare_gimp_weiss.py` — see
`inputU0NeqExactDiag.ain` for the reference setup.

## Worked example

The reduced/fast Tsuji reference test,
`cincuenta/TestSuite/inputs/inputNeqTsujiExactDiag.ain` (a full-resolution
version, `inputNeqTsuji.ain`, is what's actually registered in CTest —
identical physics, `NtNeq=200` instead of `50`):

```
##Ainur1.0

# ---- Equilibrium parameters ------------------------------------------------
FicticiousBeta=16;
ChemicalPotential=0.;
Matsubaras=200;
LatticeGf="energy,semicircular,4";
NumberOfBathPoints=5;

DmftNumberOfIterations=20;
DmftTolerance=0.001;
ImpuritySolver="exactdiag";

FitOptions=particleholesymmetric;
MinParamsDelta=0.01;
MinParamsMaxIter=10000;
MinParamsDelta2=0.01;
MinParamsTolerance=1e-3;
MinParamsVerbose=1;

TargetElectronsUp=3;
TargetElectronsDown=3;

int ImpuritySite=0;
real HubbardU=0.;

# ---- DMRG omega-scan settings (required by ParamsDmftSolver) --------------
RootOutputname="NeqTsujiExactDiag";
InfiniteLoopKeptStates=200;
matrix FiniteLoopsGs=[[@auto, 200, 0],[@auto, 200, 0]];

real OmegaBegin=-6.;
integer OmegaTotal=120;
real OmegaStep=0.1;
real OmegaDelta=0.1;

integer TridiagSteps=400;
real TridiagEps=1e-9;
TruncationTolerance="1e-10,100";
CorrectionVectorEta=0.;
GsWeight=0.1;

matrix FiniteLoopsOmega=[[@auto, 200, 2],[@auto, 200, 2]];

# ---- Non-equilibrium DMFT parameters ---------------------------------------
HubbardUFinal=2.;
TmaxNeq=5.0;
NtNeq=50;

NeqDmftIter=1;
NeqDmftTolerance=0.001;

# ---- Exact-diag non-equilibrium solver (default, NeqSolver not set) -------
# (No NeqSolver line: falls through to ImpuritySolverNeqExactDiag)
```

Note `TargetElectronsUp=3`/`TargetElectronsDown=3` against
`NumberOfBathPoints=5`: 3 + 3 = 6 = 5 + 1, satisfying the half-filling
constraint.

Physics: at half filling, 5-bath fit, an interaction quench from U_i = 0 to
U_f = 2 is followed for t ∈ [0, 5] in 50 steps.

### Sanity-check variant (U = 0, no quench)

`inputU0NeqExactDiag.ain` sets `HubbardUFinal=0.` (same as `HubbardU=0.`) —
with Σ ≡ 0 exactly, G_imp must equal the Weiss field G₀ by the Dyson
equation, which is a strong, solver-independent correctness check to run
before trusting any real quench result. Worth running on any new input setup
before changing `HubbardUFinal=` to something nonzero.

### GBEK atomic-limit variant (different solver)

`inputNeqAtomicLimitGBEKL3.ain` sets `NeqSolver="gbek"`, `NeqBathRank=3`, and
`NeqAtomicLimit=1` together, plus a hopping quench (`BandwidthFinal=`,
`QuenchShape="cosine"`, `QuenchDuration=`) — this is the GBEK Fig. 3 setup,
not the plain-ED path this guide otherwise covers. Included here mainly to
show `NeqAtomicLimit=1` used in the one context where it actually does
something (see the gotcha above). For anything beyond copying this file
verbatim, read `cincuenta/TestSuite/gbek_reference/README.md` first.
