#!/usr/bin/env python3
"""
Reproduce the inputs for GBEK PRB 88, 235106 (2013) Fig. 9: <Ekin(t)>,
<Eint(t)>, <Etot(t)> ("test of energy conservation") at U=2 and U=4, tq=0.25,
dt=0.04, t_max=4 (N=100), for Lbath = 4, 6, 8 (L=2,3,4), Cholesky
decomposition only (paper: "All data were obtained with the Cholesky
decomposition").

Also runs a U=0 case (not part of the Fig. 9 plot itself) purely as a
regression anchor: the paper states <Etot>(t) = <Ekin>(t) = 0 for all t at
U=0 (see check_energy_conservation.py for the same closed-form check at
smaller N) -- worth re-checking at the actual Fig. 9 grid size before
trusting U=2/U=4 results computed the same way.

Converges run_self_consistency first (record_history=False -- Fig. 9 only
needs the converged V, not the per-iteration history), then does one more
forward pass via compute_energy_observables to extract the energy
observables from that converged V.

Saves one .npz per (U, L): fig9_energy_U{U}_L{L}.npz, keys:
  ts, Ekin_t, Eint_t, Etot_t, d_t: (N+1,) time series
  U, L, Lbath, tstar_f, tq: scalars, for the plotting script's labels
  max_im_ekin: sanity-check scalar (should be ~0 for a physical energy)

Usage:
    uv run --with numpy --with scipy python3 run_fig9_scan.py
"""
import argparse
import time

import numpy as np

from gbek_selfconsistency import run_self_consistency, compute_energy_observables
from provenance import write_provenance

# (U, L) -- Lbath = 2*L. Run smallest/cheapest combos first.
COMBOS = [
    (0.0, 2), (0.0, 3), (0.0, 4),  # validation-only anchor, not in the paper's Fig. 9
    (2.0, 2), (2.0, 3), (2.0, 4),
    (4.0, 2), (4.0, 3), (4.0, 4),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dt", type=float, default=0.04)
    ap.add_argument("--t-max", type=float, default=4.0)
    ap.add_argument("--tstar-f", type=float, default=1.0)
    ap.add_argument("--tq", type=float, default=0.25)
    ap.add_argument("--iterations", type=int, default=30,
                    help="max DMFT iterations per (U, L) (cap; tol stops early)")
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()

    N = int(round(args.t_max / args.dt))

    for U, L in COMBOS:
        Lbath = 2 * L
        print(f"\n=== U={U} L={L} (Lbath={Lbath}) N={N} ===")
        t0 = time.time()
        Lambda, V, hop_t, ts, history, _ = run_self_consistency(
            L, N, args.dt, U, args.tstar_f, args.tq,
            n_iterations=args.iterations, tol=args.tol,
            mode="cholesky", verbose=False, record_history=False)
        wall_converge = time.time() - t0
        print(f"convergence: wall time {wall_converge:.1f} s, {len(history)} iterations, "
              f"final diff={history[-1]:.3e}")
        if history[-1] >= args.tol:
            print(f"NOTE: U={U}/L={L} hit the {args.iterations}-iteration cap WITHOUT "
                  f"reaching tol={args.tol} (final diff={history[-1]:.3e}) -- a "
                  f"numerics/physics finding worth investigating, not something to "
                  f"paper over with a looser tolerance.")

        t1 = time.time()
        ts2, Ekin_t, Eint_t, Etot_t, d_t, max_im = compute_energy_observables(
            L, N, args.dt, U, args.tstar_f, args.tq, V, verbose=True)
        wall_energy = time.time() - t1
        print(f"energy pass: wall time {wall_energy:.1f} s, max|Im(Ekin)|={max_im:.3e}")

        out = f"fig9_energy_U{U:g}_L{L}.npz"
        np.savez(out, ts=ts2, Ekin_t=Ekin_t, Eint_t=Eint_t, Etot_t=Etot_t, d_t=d_t,
                 U=U, L=L, Lbath=Lbath, tstar_f=args.tstar_f, tq=args.tq,
                 max_im_ekin=max_im)
        print(f"Wrote {out}")

        prov = write_provenance(
            out, extra_files=["run_fig9_scan.py"],
            notes=f"U={U} L={L} Lbath={Lbath} N={N} dt={args.dt} "
                  f"tstar_f={args.tstar_f} tq={args.tq} "
                  f"iterations_used={len(history)} iterations_cap={args.iterations} "
                  f"final_diff={history[-1]:.3e} max_im_ekin={max_im:.3e} "
                  f"wall_time_converge_s={wall_converge:.1f} "
                  f"wall_time_energy_s={wall_energy:.1f}")
        print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
