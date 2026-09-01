# Garcia--Hallberg--Rozenberg 2004 star-DMFT jobs

These are cluster workloads for comparing cincuenta's finite **star bath** with
the equilibrium Bethe-DMFT results of:

> D. J. Garcia, K. Hallberg, and M. J. Rozenberg, *Dynamical Mean Field Theory
> with the Density Matrix Renormalization Group*, Phys. Rev. Lett. 93, 246403
> (2004), [arXiv:cond-mat/0403169](https://arxiv.org/abs/cond-mat/0403169v1).

They are physics-motivated scale tests, not an exact reproduction. The paper
uses two continued-fraction chains with as many as 45 total sites. Cincuenta
fits a discrete star bath on the Matsubara axis. Equal chain and star site
counts therefore are not equivalent finite-size approximations; compare
convergence of Green functions, low-energy weight, and gaps instead.

## Cases and suggested order

| Case | Purpose | Bath/total sites | Kept states | Matsubaras | Real grid | Iterations | Default ranks |
|---|---|---:|---:|---:|---:|---:|---:|
| `nb7_u2p5` | pilot | 7/8 | 128 | 32 positive | 129 | 8 | 16 |
| `nb11_u2p0` | metallic control | 11/12 | 300 | 64 positive | 257 | 15 | 32 |
| `nb11_u3p5` | insulating control | 11/12 | 300 | 64 positive | 257 | 15 | 32 |
| `nb11_u2p5` | primary near-transition run | 11/12 | 300 | 64 positive | 257 | 15 | 32 |
| `nb15_u2p5` | large stress run | 15/16 | 600 | 128 positive | 513 | 20 | 64 |

Run the pilot first. Then run the two controls before interpreting the primary
`U=2.5` result. Attempt the 15-bath case only after memory and elapsed time are
known for 11 bath sites.

All inputs use:

- Bethe half-bandwidth `D=1`, represented by
  `LatticeGf="energy,semicircular,2"`;
- `FicticiousBeta=64`;
- the explicitly centered interaction convention, `ChemicalPotential=0`,
  consistent with the paper and with `ModelParams` placing `-U/2` on the
  impurity;
- equal up/down particle sectors at half filling;
- `FitOptions=particleholesymmetric`, which enforces mirrored bath energies
  and hybridizations around a zero-energy bath site;
- a constrained initial star bath whose expanded hybridizations satisfy
  `sum(V_i^2)=1/4=t^2`.

These jobs require the centered-equilibrium particle-hole fix merged into this
branch. In that corrected convention, the impurity Hamiltonian already contains
the `-U/2` shift and the fit therefore requires `ChemicalPotential=0`. Each
`InitBathVector` contains `(Nb+1)/2` independent hybridizations followed by
`(Nb-1)/2` positive bath energies; cincuenta expands it to `2*Nb` scalar
parameters describing `Nb` symmetric bath sites.

## Build

Configure a clean cluster build with the site's MPI compiler and dependencies.
For example, from the repository root:

```bash
cmake -S . -B build -DDMRGPP_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8 --target cincuenta
```

The launcher defaults to `build/cincuenta/src/cincuenta`. Set
`CINCUENTA_EXE=/absolute/path/to/cincuenta` if the cluster build is elsewhere.
Do not invoke an MPI-linked executable directly; even a one-rank run should go
through `srun` or `mpiexec`.

## Interactive or existing allocation

`run_case.sh` detects Slurm and uses `srun`; elsewhere it uses `mpiexec`.
It relies on the site's configured `srun` MPI integration; set the site-required
`SLURM_MPI_TYPE` (for example, `pmix`) in the job environment when necessary.
Every invocation creates a timestamped directory under
`runs/garcia2004_star/` and records the commit, executable/input checksums,
launcher, environment, elapsed time, and exit status.

```bash
# In a 16-rank Slurm allocation:
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u2p5 16

# In a 32-rank allocation:
cincuenta/jobs/garcia2004_star/run_case.sh nb11_u2p5 32

# Non-Slurm MPI, with an out-of-tree executable:
CINCUENTA_EXE=/path/to/build/cincuenta/src/cincuenta \
  cincuenta/jobs/garcia2004_star/run_case.sh nb11_u2p0 32
```

Useful overrides are documented by:

```bash
cincuenta/jobs/garcia2004_star/run_case.sh --help
```

In particular, `RUN_ROOT` can point at a large scratch filesystem. Keep each
run in its own directory: DMRG restart and frequency filenames are not safe to
share between simultaneous cases.

## Slurm batch submission

`slurm_job.sh` contains deliberately minimal portable directives. Supply the
site-specific account, partition, node count, memory, and wall time on the
`sbatch` command line or in a private copy. Command-line options override the
16-task defaults in the file.

```bash
# Pilot
sbatch --ntasks=16 --time=04:00:00 \
  --export=ALL,CASE=nb7_u2p5,RANKS=16 \
  cincuenta/jobs/garcia2004_star/slurm_job.sh

# Primary run (add the site's --nodes/--partition/--account/--mem options)
sbatch --ntasks=32 --time=24:00:00 \
  --export=ALL,CASE=nb11_u2p5,RANKS=32 \
  cincuenta/jobs/garcia2004_star/slurm_job.sh

# Large stress case
sbatch --ntasks=64 --time=48:00:00 \
  --export=ALL,CASE=nb15_u2p5,RANKS=64 \
  cincuenta/jobs/garcia2004_star/slurm_job.sh
```

For multiple nodes, request enough nodes for the task count according to the
site's tasks-per-node policy. The script deliberately does not guess this.

## Outputs and quick comparison

Important files in each run directory are:

- `manifest.txt`: provenance, launch command, elapsed time, and status;
- `cincuenta.log`: rank-zero solver output, final bath parameters, and DMFT
  iteration errors;
- `gimp_dmrg.txt`: final Matsubara impurity Green function;
- `latticeG_dmrg.txt`: final Matsubara local lattice Green function;
- `gimp_dmrg_real.txt`: broadened real-axis impurity Green function;
- `launcher.log`: launcher and stderr output.

Summarize one run or every child of the run root:

```bash
python3 cincuenta/jobs/garcia2004_star/summarize_runs.py \
  runs/garcia2004_star
```

The summary marks a run complete only when its manifest reports a successful
exit and all three Green-function tables have their expected row counts. It
returns a nonzero status if any selected run is incomplete; use
`--allow-incomplete` only when monitoring work still in progress.

The reported `A_nearest_zero` is `-Im G(omega)/pi` at the real-grid point
nearest zero. Treat it as a broadened finite-bath diagnostic, not a direct
critical-order parameter. The more robust initial comparison is the low-
frequency shape of `-Im G(i omega_n)` between `U=2.0`, `2.5`, and `3.5`.

## Scaling interpretation and stopping criteria

MPI distributes independent correction-vector frequencies. Rank 0 still runs
each ground-state DMRG and assembles the result, so frequency scaling cannot
remove that serial cost. Using more ranks than useful independent frequencies
will not help.

Before advancing to the next case, require:

1. successful exit and finite three-column Green-function files;
2. decreasing DMFT errors, preferably convergence below `1e-4`;
3. approximately particle-hole-symmetric final bath parameters and spectra;
4. physically distinct metallic (`U=2.0`) and insulating (`U=3.5`) controls;
5. acceptable truncation behavior and no memory pressure in the DMRG log.

If a case reaches its iteration limit without convergence, retain the full run
directory. Do not silently increase bath size or interpret the real-axis line
shape before inspecting the fit and truncation diagnostics.
