#!/usr/bin/env python3
"""
Reproduce the inputs for GBEK PRB 88, 235106 (2013) Fig. 7: per-DMFT-iteration
double occupation d_it(t) = <n0up(t) n0dn(t)>, at U=5, L_bath=4 (L=2), t_max=4,
t_q=0.25, for BOTH the eigenvector and Cholesky decomposition of the Weiss
field, run for exactly 7 iterations (the paper's own choice) regardless of
whether either mode has numerically converged by then -- Fig. 7's whole point
is to show the shape of the iteration-to-iteration CHANGE, not just the
converged answer (that's Fig. 8).

Saves one .npz per mode: fig7_docc_L{L}_{mode}.npz, keys:
  docc_history: (7, N+1) real array, d_it(t) for it=1..7
  ts: (N+1,) time grid

Usage:
    uv run --with numpy --with scipy python3 run_fig7_scan.py
"""
import argparse
import time

import numpy as np

from gbek_selfconsistency import run_self_consistency
from provenance import write_provenance

MODES = ("cholesky", "eigenvector")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--L", type=int, default=2, help="Cholesky/eigenvector rank (Lbath = 2L)")
    ap.add_argument("--N", type=int, default=100, help="number of real-time steps")
    ap.add_argument("--dt", type=float, default=0.04)
    ap.add_argument("--U", type=float, default=5.0)
    ap.add_argument("--tstar-f", type=float, default=1.0)
    ap.add_argument("--tq", type=float, default=0.25)
    ap.add_argument("--iterations", type=int, default=7,
                    help="fixed iteration count (paper's Fig. 7 uses exactly 7)")
    args = ap.parse_args()

    for mode in MODES:
        print(f"\n=== mode={mode} ===")
        t0 = time.time()
        Lambda, V, hop_t, ts, history, docc_history = run_self_consistency(
            args.L, args.N, args.dt, args.U, args.tstar_f, args.tq,
            n_iterations=args.iterations, tol=1e-30,  # never early-break: want all `iterations`
            mode=mode, verbose=True, record_history=True)
        print(f"wall time: {time.time()-t0:.1f} s")
        print("convergence history:", [f"{h:.2e}" for h in history])

        if len(docc_history) != args.iterations:
            print(f"WARNING: expected {args.iterations} recorded iterations, "
                  f"got {len(docc_history)} -- tol triggered an early break "
                  f"despite tol=1e-30; investigate before trusting this mode's curves.")

        docc_stack = np.stack(docc_history)  # (iterations, N+1)
        out = f"fig7_docc_L{args.L}_{mode}.npz"
        np.savez(out, docc_history=docc_stack, ts=ts)
        print(f"Wrote {out}  shape={docc_stack.shape}")

        prov = write_provenance(
            out, extra_files=["run_fig7_scan.py"],
            notes=f"mode={mode} L={args.L} N={args.N} dt={args.dt} U={args.U} "
                  f"tstar_f={args.tstar_f} tq={args.tq} iterations={args.iterations} "
                  f"n_recorded={len(docc_history)} final_diff={history[-1]:.3e}")
        print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
