#!/usr/bin/env python3
"""
Reproduce the inputs for GBEK PRB 88, 235106 (2013) Fig. 8: CONVERGED double
occupation d(t) = <n0up(t) n0dn(t)> at U=5, t_q=0.25, dt=0.04 (fixed per the
paper's Fig. 8 caption), for four (L_bath, t_max) combinations, under BOTH the
eigenvector and Cholesky decomposition of the Weiss field:

  combo A: L=2 (Lbath=4), t_max=4  (N=100)
  combo B: L=2 (Lbath=4), t_max=2  (N=50)
  combo C: L=3 (Lbath=6), t_max=2  (N=50)
  combo D: L=4 (Lbath=8), t_max=2  (N=50)

Run in increasing-cost order (B, C, A, D) so partial results are available
early. Convergence uses tol=1e-6 (unlike Fig. 7, which wants every iteration
regardless of convergence); any combo/mode that hits the iteration cap
without reaching tol is flagged explicitly, not silently accepted.

Uses record_history=True and takes docc_history[-1] (the converged
iteration's curve) rather than a special-cased re-propagation, reusing
run_self_consistency's existing code path with no new logic.

Saves one .npz per (combo, mode): fig8_docc_{combo}_{mode}.npz, keys:
  d_t: (N+1,) converged double occupation
  ts: (N+1,) time grid
  L, t_max: scalars, for the plotting script's labels

Usage:
    uv run --with numpy --with scipy python3 run_fig8_scan.py
"""
import argparse
import time

import numpy as np

from gbek_selfconsistency import run_self_consistency
from provenance import write_provenance

MODES = ("cholesky", "eigenvector")

# (combo_label, L, t_max) -- run in increasing-cost order.
COMBOS = [
    ("B", 2, 2.0),
    ("C", 3, 2.0),
    ("A", 2, 4.0),
    ("D", 4, 2.0),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dt", type=float, default=0.04)
    ap.add_argument("--U", type=float, default=5.0)
    ap.add_argument("--tstar-f", type=float, default=1.0)
    ap.add_argument("--tq", type=float, default=0.25)
    ap.add_argument("--iterations", type=int, default=30,
                    help="max DMFT iterations per combo/mode (cap; tol stops early)")
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()

    for label, L, t_max in COMBOS:
        N = int(round(t_max / args.dt))
        for mode in MODES:
            print(f"\n=== combo {label}: L={L} t_max={t_max} N={N} mode={mode} ===")
            t0 = time.time()
            Lambda, V, hop_t, ts, history, docc_history = run_self_consistency(
                L, N, args.dt, args.U, args.tstar_f, args.tq,
                n_iterations=args.iterations, tol=args.tol,
                mode=mode, verbose=False, record_history=True)
            wall = time.time() - t0
            print(f"wall time: {wall:.1f} s, {len(history)} iterations, "
                  f"final diff={history[-1]:.3e}")

            if history[-1] >= args.tol:
                print(f"NOTE: combo {label}/{mode} hit the {args.iterations}-iteration "
                      f"cap WITHOUT reaching tol={args.tol} (final diff={history[-1]:.3e}) "
                      f"-- this is a numerics/physics finding worth investigating, "
                      f"not something to paper over with a looser tolerance.")

            d_t = docc_history[-1]
            out = f"fig8_docc_{label}_{mode}.npz"
            np.savez(out, d_t=d_t, ts=ts, L=L, t_max=t_max)
            print(f"Wrote {out}")

            prov = write_provenance(
                out, extra_files=["run_fig8_scan.py"],
                notes=f"combo={label} mode={mode} L={L} N={N} t_max={t_max} "
                      f"dt={args.dt} U={args.U} tstar_f={args.tstar_f} tq={args.tq} "
                      f"iterations_used={len(history)} iterations_cap={args.iterations} "
                      f"final_diff={history[-1]:.3e} wall_time_s={wall:.1f}")
            print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
