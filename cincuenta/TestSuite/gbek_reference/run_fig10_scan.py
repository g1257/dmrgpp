#!/usr/bin/env python3
"""
Reproduce the inputs for GBEK PRB 88, 235106 (2013) Fig. 10: converged
double occupation d(t) vs U = 0, 1, 2, 4, 6, 8, at tq=0.25, dt=0.04,
t_max=4 (N=100), for Lbath = 4, 6, 8 (L=2,3,4), Cholesky decomposition only
(paper: "All data were obtained with the Cholesky decomposition").

Modeled directly on run_fig8_scan.py: uses record_history=True and takes
docc_history[-1] (the converged iteration's curve), reusing
run_self_consistency's existing code path with no new logic -- unlike
Fig. 9, this figure needs no new physics machinery.

Run in increasing-cost order (small U, small L first) so partial results are
available early.

Watch for: at large U (6, 8), the paper's own text describes a
collapse-revival oscillation with period ~1/U -- at U=8 that's ~0.125, only
~3 grid points per oscillation at dt=0.04. This is a resolution risk, not
just a slow-convergence risk; plot_docc_scan.py --figure 10 checks the
observed period against 1/U rather than assuming it's resolved.

Saves one .npz per (U, L): fig10_docc_U{U}_L{L}.npz, keys:
  d_t, ts, v: (N+1,) converged double occupation, time grid, ramp v(t)
  U, L, Lbath: scalars, for the plotting script's labels

Usage:
    uv run --with numpy --with scipy python3 run_fig10_scan.py
"""
import argparse
import time

import numpy as np

from gbek_selfconsistency import run_self_consistency, v_ramp
from provenance import write_provenance

U_VALUES = (0.0, 1.0, 2.0, 4.0, 6.0, 8.0)
L_VALUES = (2, 3, 4)  # Lbath = 4, 6, 8


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
    combos = sorted([(U, L) for U in U_VALUES for L in L_VALUES],
                     key=lambda c: (c[0], c[1]))  # smallest U, smallest L first

    for U, L in combos:
        Lbath = 2 * L
        print(f"\n=== U={U} L={L} (Lbath={Lbath}) N={N} ===")
        t0 = time.time()
        Lambda, V, hop_t, ts, history, docc_history = run_self_consistency(
            L, N, args.dt, U, args.tstar_f, args.tq,
            n_iterations=args.iterations, tol=args.tol,
            mode="cholesky", verbose=False, record_history=True)
        wall = time.time() - t0
        print(f"wall time: {wall:.1f} s, {len(history)} iterations, "
              f"final diff={history[-1]:.3e}")
        if history[-1] >= args.tol:
            print(f"NOTE: U={U}/L={L} hit the {args.iterations}-iteration cap WITHOUT "
                  f"reaching tol={args.tol} (final diff={history[-1]:.3e}) -- a "
                  f"numerics/physics finding worth investigating, not something to "
                  f"paper over with a looser tolerance.")

        d_t = docc_history[-1]
        v = np.array([v_ramp(t, args.tq) for t in ts])
        out = f"fig10_docc_U{U:g}_L{L}.npz"
        np.savez(out, d_t=d_t, ts=ts, v=v, U=U, L=L, Lbath=Lbath)
        print(f"Wrote {out}")

        prov = write_provenance(
            out, extra_files=["run_fig10_scan.py"],
            notes=f"U={U} L={L} Lbath={Lbath} N={N} dt={args.dt} "
                  f"tstar_f={args.tstar_f} tq={args.tq} "
                  f"iterations_used={len(history)} iterations_cap={args.iterations} "
                  f"final_diff={history[-1]:.3e} wall_time_s={wall:.1f}")
        print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
