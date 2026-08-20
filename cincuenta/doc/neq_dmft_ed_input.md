# Non-equilibrium DMFT input: the exact-diagonalization (ED) solver

This guide is for someone running `cincuenta` to do a non-equilibrium DMFT
(neq-DMFT) interaction-quench calculation, without needing to read the C++.
It assumes you can already run an equilibrium `cincuenta` job — see
[`cincuenta_design.md`](cincuenta_design.md) for that, including the
completed equilibrium **Input File Reference** table.

If you only need to get a run going quickly, skip to
[Worked examples](#worked-examples) and come back to the parameter tables as
needed.

## The physical picture

`cincuenta` first solves the **equilibrium** DMFT problem to get a converged
bath (the `{V_α, ε_α}` in `cincuenta_design.md`'s Anderson model).
Then, in an optional second stage, it can do limited neq-DMFT calculations.
It quenches either the interaction from `U_i` to `U_f`
at t = 0 or the hopping from ? to ? with an optional ramp function.
The code then follows the resulting real-time dynamics on the Keldysh
contour.
This second stage is set up by adding a `TmaxNeq=` line (with a positive value) below
your usual equilibrium input — its presence is what turns on neq mode.

There are two neq **methods**, chosen by `NeqSolver=`:

| `NeqSolver=` | Meaning |
|---|---|
| `"exactdiag"` (or omit the line — this is the default) | **This guide's subject.** Exact diagonalization. |
| `"tdmrg"` | Time-dependent DMRG instead of ED. Not covered here — not fully implemented. |

### How much bath does the quench see? `NeqBathRank=`

This is the physically important choice inside the ED method, and it's worth
understanding before you pick a value:

The equilibrium bath fit gives you a bath that's correct **for the
pre-quench system**. After the quench, the true non-equilibrium bath should
respond — its effective hybridization Λ(t,t′) changes in time. `NeqBathRank=`
(default 0) controls how much of that response the calculation captures:

- **`NeqBathRank=0`** (default): the bath stays exactly as fitted at
  equilibrium for the entire real-time run — no time-dependent response at
  all. This is crude approximation, it drops the back-reaction of the
  quench on the effective medium entirely. It is only "defensible"
  when the quench is small enough that the bath wouldn't change
  change.
- **`NeqBathRank=L>0`**: adds a second, low-rank ("Cholesky")
  bath that *does* evolve with the quench, following Gramsch, Balzer, Eckstein, Kollar,
  PRB **88**, 235106 (2013) ("GBEK"). This is the physically meaningful way
  to run neq-DMFT ED — larger L captures the bath's time-dependent
  response for longer before the finite-rank truncation error becones significant (see
  `cincuenta/TestSuite/gbek_reference/README.md` for how much `L` you need
  and for how long, for a given system).

**Recommendation: unless you have a specific reason to check the
`NeqBathRank=0` limit (e.g. as a small-quench sanity check), set
`NeqBathRank=` to a positive value.** `NeqBathRank=0` is not a "simpler
version of the same method" — it's a different, more restrictive
approximation that happens to be what you get if you don't set the parameter
at all.

## Required equilibrium input (recap)

In general doing non-equilibrium DMFT depends on starting from and
equilibrium calculation. so except in a special case covered later the usual
required equilibrium keywords still apply — see
[`cincuenta_design.md`](cincuenta_design.md#input-file-reference) for the
full table. `ImpuritySolver=` should be `"exactdiag"` (the equilibrium bath
fit is what the neq stage's bath comes from).

**`LatticeGf=` constraint for any neq run:** the neq Bethe-lattice
self-consistency only supports the semicircular DOS — `LatticeGf=` must be
exactly `"energy,semicircular,W"`. The equilibrium stage itself supports
more (e.g. `"momentum,1D,Nk"`), Support for other LatticeGf is not on
the roadmap yet.

Two keywords the ED solver itself requires directly, easy to miss if you're
scanning `Params*.h` instead of the actual `.ain` examples:

| Keyword | Type | Meaning |
|---|---|---|
| `TargetElectronsUp=` | int | Number of up-spin electrons in the target sector |
| `TargetElectronsDown=` | int | Number of down-spin electrons in the target sector |

## Neq-specific input

All of the following go below the equilibrium block.

### Required

| Keyword | Type | Meaning |
|---|---|---|
| `HubbardU=` | real | U_i, the **initial** (equilibrium) interaction. Same keyword the equilibrium block already needs — not duplicated, just reused. |
| `HubbardUFinal=` | real | U_f, the interaction after the t = 0 quench. Set equal to `HubbardU=` for a no-quench sanity check. |
| `TmaxNeq=` | real | Total real-time window to simulate, t ∈ [0, TmaxNeq]. **Presence of this line (with a positive value) is what turns on neq mode at all.** |
| `NtNeq=` | int | Number of real-time steps. `dt = TmaxNeq / NtNeq` is fixed automatically — you don't set `dt` directly. |

### Optional

| Keyword | Type | Default | Meaning |
|---|---|---|---|
| `NeqSolver=` | string | `"exactdiag"` | `"exactdiag"` or `"tdmrg"`; see method table above. Anything else is a hard error. |
| `NeqBathRank=` | int | 0 | rank of the Cholesky decompsition used for the evolving bath and bath sites is this * 2 |
| `NeqDmftIter=` | int | same as equilibrium `DmftNumberOfIterations=` | Inner self-consistency iterations at each time step, used only when `NeqBathRank>0` (updating the second bath); harmless no-op at `NeqBathRank=0`. |
| `NeqDmftTolerance=` | real | same as equilibrium `DmftTolerance=` | Convergence tolerance for the above. |
| `BandwidthFinal=` | real | 0 (no quench) | Bethe-lattice bandwidth *for t > 0*, i.e. a **hopping** quench on top of (or instead of) the interaction quench. 0 means the lattice bandwidth stays at its equilibrium `LatticeGf=` value for all t. |
| `QuenchShape=` | string | `"step"` | Ramp shape for the hopping quench: `"step"`, `"cosine"`, or `"tanh"`. Only matters if `BandwidthFinal=` is nonzero. |
| `QuenchDuration=` | real | 0 (instantaneous) | Ramp duration t_q for the hopping quench. |
| `NeqOutputPrefix=` | string | `""` | Prefix prepended to output Green's-function filenames (see [Output files](#output-files)). |
| `NeqAtomicLimit=` | int (0/1) | 0 | Start the neq run from the true atomic limit no bath at all at t = 0 (hybridization Λ⁻ ≡ 0 exactly), per GBEK Sec. VI. See [The atomic limit skips the equilibrium stage entirely](#the-atomic-limit-skips-the-equilibrium-stage-entirely) below — this changes what's *allowed* in the rest of the file, not just what's optional. **Not supported with `NeqSolver="tdmrg"`** (rejected with an error). Requires `TargetElectronsUp + TargetElectronsDown == 1`. |

### When NeqAtomicLimit = 1 the equilibrium stage is skipped

As a result the following equilibrium are not allowed in input and
result in an error.

- `NumberOfBathPoints=`
- `ImpuritySolver=`
- `FitOptions=`

`LatticeGf` is still required and for the atomic limit starting point
 "energy,semicircular,0" is the required input. This is for
 consistency with all other inputs for cincuenta.

The rest of the equilibrium-only keywords (`MinParams*=`, `InitBath*=`,
`Precision=`, and the whole DMRG omega-scan block: `RootOutputname=`,
`FiniteLoopsGs=`, `OmegaBegin=`, etc.) aren't actively rejected, but they're
never read either — just omit them. ` In the future logical
consistency of inputs may be further constrained so leaving
inconsistent and/or unused keywords in inputs is not guaranteed to be
supported in future versions.

`FicticiousBeta=`, `ChemicalPotential=`,
`Matsubaras=`, `TargetElectronsUp=`/`Down=`, `ImpuritySite=`, and `HubbardU=`
are still required — those define the real-time/Matsubara grid and the
target sector, which the neq stage genuinely needs regardless of whether an
equilibrium bath fit ran.

`ChemicalPotential` must equal U/2 since only the half filling
particle hole symmetric case is implemented.

See [Required equilibrium input](#required-equilibrium-input-recap) above and
`cincuenta_design.md`'s reference table. This trips people up because it
looks redundant, but omitting it will fail input parsing before you even get
to the neq stage. **Except under `NeqAtomicLimit=1`**, where this whole block
is omitted instead — see
[The atomic limit skips the equilibrium stage entirely](#the-atomic-limit-skips-the-equilibrium-stage-entirely).

## Output files

After the neq run finishes, files are written with `<prefix>` = your
`NeqOutputPrefix=` value (or nothing, if omitted):

| File | Contents |
|---|---|
| `<prefix->green-retarded` | Retarded impurity Green's function G^R(t, t′) |
| `<prefix->green-lesser` | Lesser impurity Green's function G^<(t, t′) |
| `<prefix->green-matsubara-t` | Mixed Matsubara/real-time component |
| `<prefix->weiss-green-retarded/lesser/matsubara-t` | Same three components for the Weiss field G₀ (the bath-only, U = 0 reference) |
| `<prefix->weiss-delta-...` | Hybridization Λ(t, t′) components. (The filenames say "delta" — that's legacy output-file naming; the physical quantity is Λ, not Δ.) |

If `NeqBathRank>0`, two more files appear: `<prefix->cholesky-V` (the raw
second-bath Cholesky factor) and `<prefix->docc-energy` (double occupation
and kinetic/interaction/total energy vs. time).

Two boundary-condition identities are checked automatically in CI
(`cincuenta/src/CMakeLists.txt`, `neqTsuji0` test): `G^R(0,0) = -i` and
`G^M(τ=0) = -n_imp`. If your own run's outputs don't satisfy these, something
upstream is wrong before you get to physics interpretation.

For the U_i = U_f = 0 sanity check (Σ = 0 exactly ⟹ G_imp must equal the
Weiss field G₀), use `cincuenta/TestSuite/compare_gimp_weiss.py` — see
`inputU0NeqExactDiag.ain` for the reference setup.

## Worked examples

### Small-quench sanity check (`NeqBathRank=0`, the default)

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

# NeqSolver not set -> "exactdiag" default, NeqBathRank not set -> 0 (fixed bath).
# Reasonable here because the quench starts from U_i=0.
```

Note `TargetElectronsUp=3`/`TargetElectronsDown=3` against
`NumberOfBathPoints=5`: 3 + 3 = 6 = 5 + 1, satisfying the half-filling
constraint.

Physics: at half filling, 5-bath fit, an interaction quench from U_i = 0 to
U_f = 2 is followed for t ∈ [0, 5] in 50 steps, with the bath fixed at its
equilibrium value — reasonable here specifically because U_i = 0 means the
"equilibrium bath" is trivial to begin with.

### Sanity-check variant (U = 0, no quench)

`inputU0NeqExactDiag.ain` sets `HubbardUFinal=0.` (same as `HubbardU=0.`) —
with Σ ≡ 0 exactly, G_imp must equal the Weiss field G₀ by the Dyson
equation, which is a strong, solver-independent correctness check to run
before trusting any real quench result.

### Physically meaningful quench (`NeqBathRank>0`)

This is the setup to reach for when the quench isn't small — read
`cincuenta/TestSuite/gbek_reference/README.md` for guidance on choosing `L`
and for how long the result stays trustworthy before finite-rank error
dominates.

Right now only the special case of the atomic limit is validated:

### True atomic limit (`NeqAtomicLimit=1`)

`inputNeqAtomicLimitGBEKL3.ain` sets `NeqAtomicLimit=1` together with
`NeqBathRank=3` — this is the GBEK Fig. 3 setup starting from an exact
atomic limit (no bath at all at t = 0) rather than the near-zero-bandwidth
approximation the previous example uses. Note how much shorter this file is
than the others in this guide — equilibrium only input can be dropped,
but `LatticeGf` is still required.

```
##Ainur1.0

LatticeGf="energy,semicircular,0";
FicticiousBeta=16;
ChemicalPotential=1.0;
Matsubaras=200;

TargetElectronsUp=1;
TargetElectronsDown=0;

int ImpuritySite=0;
real HubbardU=2.;

HubbardUFinal=2.;
TmaxNeq=4.0;
NtNeq=100;

NeqDmftIter=3;
NeqDmftTolerance=1e-8;

NeqSolver="exactdiag";
NeqBathRank=3;
NeqAtomicLimit=1;

BandwidthFinal=4.;
QuenchShape="cosine";
QuenchDuration=0.25;

NeqOutputPrefix="atomic-limit-gbek-L3";
```

`TargetElectronsUp=1`/`TargetElectronsDown=0` (single-atom sector) is
required — the neq solver itself needs it (`nup + ndown == 1` for the
atomic-limit bypass), independent of the equilibrium stage this input
doesn't have. For anything beyond copying this file verbatim, read
`cincuenta/TestSuite/gbek_reference/README.md` first.
