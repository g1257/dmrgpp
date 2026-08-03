#!/usr/bin/env python3
"""
Reproduce GBEK Fig. 3's bottom-left panel (err^step(t), Eq. 67) directly
against our own exact atomic-limit self-consistency reference
(gbek_selfconsistency.py's own dump_lesser output), at the paper's own
ranks (L=2, L=3).

IMPORTANT convention note: this reads the target with `re + 1j*im` directly
(no extra -1j factor), because gbek_selfconsistency.py::dump_lesser already
writes Lambda^<_+ itself. This is a DIFFERENT convention from
plot_errstep.py's load_lambda(), which is designed for cincuenta's raw
*-weiss-delta-lesser C++ dumps (the raw lesser component, needing an extra
-i to become Lambda) -- using that loader on a gbek_selfconsistency.py dump
silently rotates every value by an erroneous extra factor of -i, which
looks like nonsense (e.g. two different ranks appearing to give identical,
wrong reconstructions) rather than raising an error. Do not swap loaders
between the two file families.

Usage:
    uv run --with numpy --with matplotlib python3 plot_fig3_errstep.py \\
        --target gbek-atomic-limit-exact-lesser --ranks 2,3
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, reconstruct, eigenvector_decompose
from plot_errstep import errstep_curve, global_err
from provenance import write_provenance

DEFAULT_TARGET = "gbek-atomic-limit-exact-lesser"


def load_self_consistency_lambda(path):
    """Read a gbek_selfconsistency.py::dump_lesser output directly as
    Lambda^<_+ -- NO extra -1j (see module docstring)."""
    ts, re, im = read_lesser_file(path)
    return ts, re + 1j * im


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default=DEFAULT_TARGET,
                    help="gbek_selfconsistency.py dump_lesser output (NOT a C++ "
                         "weiss-delta-lesser dump -- see module docstring)")
    ap.add_argument("--ranks", default="2,3", help="comma-separated Cholesky/eigenvector ranks")
    ap.add_argument("--out", default="fig3_errstep.png")
    ap.add_argument("--title", default="cf. GBEK Fig. 3 bottom-left panel")
    args = ap.parse_args()

    ranks = [int(x) for x in args.ranks.split(",")]
    ts, lam = load_self_consistency_lambda(args.target)

    full = np.zeros_like(lam)
    N = lam.shape[0] - 1
    for n in range(N + 1):
        for j in range(N + 1):
            full[n, j] = lam[n, j] if j <= n else np.conj(lam[j, n])

    fig, ax = plt.subplots(figsize=(7, 5.8))
    colors = plt.cm.viridis(np.linspace(0.15, 0.75, len(ranks)))
    print(f"{'curve':<16} {'err^step(t=4)':>14} {'global err[A]':>14}")

    for L, c in zip(ranks, colors):
        V_ch = cholesky_causal(lam, L)
        ch_recon = reconstruct(V_ch)
        V_ev = eigenvector_decompose(full, L)
        ev_recon = reconstruct(V_ev)

        ch_curve = errstep_curve(full, ch_recon)
        ev_curve = errstep_curve(full, ev_recon)

        ax.semilogy(ts, ch_curve, "-", color=c, lw=2.0, label=f"ch rank {L}")
        ax.semilogy(ts, ev_curve, "--", color=c, lw=1.6, alpha=0.8, label=f"ev rank {L}")

        print(f"{'ch L='+str(L):<16} {ch_curve[-1]:>14.4f} {global_err(full, ch_recon):>14.4f}")
        print(f"{'ev L='+str(L):<16} {ev_curve[-1]:>14.4f} {global_err(full, ev_recon):>14.4f}")

    # Match the paper's own Fig. 3 bottom-left axis exactly: log scale from
    # 1e-6 to 1 (1e0). Curves below 1e-6 are numerical noise from the causal
    # recursion's floating-point floor, not signal; clipping them into the
    # paper's own range (rather than auto-scaling down to ~1e-18) is what
    # makes this panel directly, visually comparable to the paper's.
    ax.set_ylim(1e-6, 1)
    ax.set_xlim(0, ts[-1])
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\mathrm{err}^{\mathrm{step}}(t)$  (Eq. 67)")
    ax.set_title(args.title, fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2, which="both")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nWrote {args.out}")

    prov = write_provenance(args.out, extra_files=[args.target],
                             notes=f"target={args.target} ranks={ranks}")
    print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
