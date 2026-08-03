#!/usr/bin/env python3
"""
Follow-up to Fig. 8 combo A (L=2/Lbath=4, t_max=4): our own Cholesky curve
for that combo shows the same qualitative late-time upturn (after a mid-time
dip) that the paper's own Fig. 8(b) green-dashed L_bath=4/t_max=4 curve shows
-- but ours climbs much higher (~0.35 vs. the paper's ~0.06 by t=4).

This script runs the same t_max=4 protocol at L=4 (Lbath=8) instead of L=2,
to see whether the magnitude gap shrinks at higher rank (consistent with
plain rank-insufficiency, matching the paper's own "SIAM representation
breaks down" caveat for t gtrsim 3) or persists at the same disproportionate
scale (which would point at something specific to our Cholesky optimal-
update solve rather than generic rank-insufficiency). Also directly useful
groundwork for Fig. 9, which needs Lbath=8 runs regardless.

Saves fig8_docc_E_{mode}.npz in the same format as run_fig8_scan.py's
output (combo label "E", not one of the paper's original 4), so
plot_docc_scan.py's existing loading/plotting logic can be reused by adding
"E" to its COMBOS list if desired.

Usage:
    uv run --with numpy --with scipy python3 investigate_L4_tmax4.py
"""
import time

import numpy as np

from gbek_selfconsistency import run_self_consistency
from provenance import write_provenance

L, N, dt, U, tstar_f, tq = 4, 100, 0.04, 5.0, 1.0, 0.25
MODES = ("cholesky", "eigenvector")


def main():
    for mode in MODES:
        print(f"\n=== combo E: L={L} (Lbath={2*L}) t_max={N*dt} mode={mode} ===")
        t0 = time.time()
        Lambda, V, hop_t, ts, history, docc_history = run_self_consistency(
            L, N, dt, U, tstar_f, tq, n_iterations=30, tol=1e-6,
            mode=mode, verbose=False, record_history=True)
        wall = time.time() - t0
        print(f"wall time: {wall:.1f} s, {len(history)} iterations, "
              f"final diff={history[-1]:.3e}")
        if history[-1] >= 1e-6:
            print(f"NOTE: hit the 30-iteration cap WITHOUT reaching tol=1e-6 "
                  f"(final diff={history[-1]:.3e}) -- flag, don't ignore.")

        d_t = docc_history[-1]
        out = f"fig8_docc_E_{mode}.npz"
        np.savez(out, d_t=d_t, ts=ts, L=L, t_max=float(N * dt))
        print(f"Wrote {out}")
        print(f"d(t) at t=0,1,2,2.5,3,3.5,4: "
              f"{[round(float(np.interp(t, ts, d_t)), 4) for t in (0,1,2,2.5,3,3.5,4)]}")

        prov = write_provenance(
            out, extra_files=["investigate_L4_tmax4.py"],
            notes=f"combo=E (investigative, not one of the paper's 4) mode={mode} "
                  f"L={L} N={N} t_max={N*dt} dt={dt} U={U} tstar_f={tstar_f} tq={tq} "
                  f"iterations_used={len(history)} final_diff={history[-1]:.3e} "
                  f"wall_time_s={wall:.1f}")
        print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
