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

The sole production neq method is the positive-rank GBEK two-bath
exact-diagonalization scheme described below.

### How much bath does the quench see? `NeqBathRank=`

This is the physically important choice inside the ED method, and it's worth
understanding before you pick a value:

The equilibrium bath fit gives you a bath that's correct **for the
pre-quench system**. After the quench, the non-equilibrium bath must respond:
its effective hybridization Λ(t,t′) changes in time. Therefore the ED/GBEK
path requires an explicit **`NeqBathRank=L>0`**. It adds a low-rank
("Cholesky") second bath that evolves with the quench, following Gramsch,
Balzer, Eckstein, Kollar, PRB **88**, 235106 (2013) ("GBEK"). Larger `L`
captures the bath's time-dependent response for longer before finite-rank
truncation error becomes significant (see
`cincuenta/TestSuite/gbek_reference/README.md` for guidance).

`NeqBathRank=0` and an omitted `NeqBathRank=` are rejected; the former
frozen-first-bath approximation is no longer a production path.

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
| `NeqBathRank=` | int | Positive rank of the Cholesky decomposition used for the evolving second bath; it has `2 * NeqBathRank` sites. Zero and omission are rejected. |

### Optional

| Keyword | Type | Default | Meaning |
|---|---|---|---|
| `NeqDmftIter=` | int | 10 | Fixed corrector iterations at each time step, updating the second bath; it must be positive when specified. |
| `BandwidthFinal=` | real | 0 (no quench) | Bethe-lattice bandwidth *for t > 0*, i.e. a **hopping** quench on top of (or instead of) the interaction quench. 0 means the lattice bandwidth stays at its equilibrium `LatticeGf=` value for all t. |
| `QuenchShape=` | string | `"step"` | Ramp shape for the hopping quench: `"step"`, `"cosine"`, or `"tanh"`. Only matters if `BandwidthFinal=` is nonzero. |
| `QuenchDuration=` | real | 0 (instantaneous) | Ramp duration t_q for the hopping quench. |
| `NeqOutputPrefix=` | string | `""` | Prefix prepended to output Green's-function filenames (see [Output files](#output-files)). |
| `NeqAtomicLimit=` | int (0/1) | 0 | Start the neq run from the true atomic limit no bath at all at t = 0 (hybridization Λ⁻ ≡ 0 exactly), per GBEK Sec. VI. See [The atomic limit skips the equilibrium stage entirely](#the-atomic-limit-skips-the-equilibrium-stage-entirely) below — this changes what's *allowed* in the rest of the file, not just what's optional. Requires `TargetElectronsUp + TargetElectronsDown == 1`. |

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
| `<prefix->lambda-...` | Total hybridization target Λ(t, t′) components. |
| `<prefix->equilibrium-gimp-matsubara` | Equilibrium Matsubara rows: frequency, real part, imaginary part. This is emitted once before neq solver dispatch for non-atomic runs. |

If `NeqBathRank>0`, two more files appear: `<prefix->cholesky-V` (the raw
second-bath Cholesky factor) and `<prefix->docc-energy` (double occupation
and kinetic/interaction/total energy vs. time).

The focused CI regression checks rank-zero rejection and verifies the
positive-rank GBEK run in `cincuenta/TestSuite/inputs/inputNeqGBEKFast.ain`,
including the boundary condition `G^R(0,0) = -i`, finite `G^{R,<}`, a nonzero
first-bath contribution, and activation of both rank-two Cholesky columns.

## Worked examples

### Fast positive-rank GBEK run (`NeqBathRank=2`)

Use `cincuenta/TestSuite/inputs/inputNeqGBEKFast.ain`. It is the focused
end-to-end regression input: a one-site fitted equilibrium bath, `HubbardU =
HubbardUFinal = 2`, `TmaxNeq=0.3`, and a rank-two evolving second bath. Its
small grid is intended for a quick execution check, not physics benchmarks.

The input has `TargetElectronsUp=1` and `TargetElectronsDown=1` for its
one-impurity-plus-one-bath-site equilibrium system and `NeqBathRank=2`.

### Physically meaningful quench (`NeqBathRank>0`)

This is the setup to reach for when the quench isn't small — read
`cincuenta/TestSuite/gbek_reference/README.md` for guidance on choosing `L`
and for how long the result stays trustworthy before finite-rank error
dominates.

The atomic-limit case provides an additional independently validated limit:

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
