#!/usr/bin/env python3
"""Plot Re[-i Λ^<(t,t')] = Im(delta.lesser(t,t')) as a 2D colormap.

Reproduces the style of GEBK Fig. 3 (Gramsch, Balzer, Eckstein, Kollar,
Phys. Rev. B 88, 235106, 2013): the "input Weiss field" -iΛ^<_+(t,t').

In the atomic limit (Λ_- ≈ 0): Λ^<_+ ≈ Λ^< = delta.lesser.
So Re[-iΛ^<_+] = Im(delta.lesser).

File format (written by KadanoffBaym::dump, prefix-lesser):
    t  t'  Re(Λ^<)  Im(Λ^<)
one line per (t,t') pair, full (nT+1)×(nT+1) matrix row-major.

Usage
-----
    python3 compare_neq_delta_lesser.py gebk-fig3-weiss-delta-lesser \\
        [--tmax 4.0] [--title "GEBK Fig. 3: Re[-i Λ^<_+]"]
"""

import argparse
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_lesser_file(path):
    """Return (t_vals, tp_vals, data_re, data_im) arrays from a lesser file."""
    rows = []
    with open(path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 4:
                continue
            try:
                rows.append([float(c) for c in cols[:4]])
            except ValueError:
                continue
    if not rows:
        print(f"ERROR: no data in {path}", file=sys.stderr)
        sys.exit(1)
    arr = np.array(rows)
    t_all  = arr[:, 0]
    tp_all = arr[:, 1]
    re_all = arr[:, 2]
    im_all = arr[:, 3]

    # Infer nT from unique t values
    t_unique = np.unique(t_all)
    nT1 = len(t_unique)  # nT+1

    re_mat = re_all.reshape(nT1, nT1)
    im_mat = im_all.reshape(nT1, nT1)
    return t_unique, re_mat, im_mat


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("lesser_file", metavar="LESSER_FILE",
                   help="*-weiss-delta-lesser file from dumpGreenFunctions()")
    p.add_argument("--output", default=None,
                   help="Output image filename (default: <lesser_file>.png)")
    p.add_argument("--tmax", type=float, default=None,
                   help="Clip t-axis to tmax (default: use full range)")
    p.add_argument("--title", default="Re[-i Λ$^<_+$(t,t')]",
                   help="Plot title")
    p.add_argument("--vmax", type=float, default=None,
                   help="Color scale maximum (default: auto)")
    args = p.parse_args()

    t_vals, re_mat, im_mat = read_lesser_file(args.lesser_file)
    nT1 = len(t_vals)
    dt  = t_vals[1] - t_vals[0] if nT1 > 1 else 1.0
    tmax = args.tmax if args.tmax is not None else t_vals[-1]

    # Re[-i Λ^<] = Im(Λ^<)
    quantity = im_mat

    # Clip to tmax if requested
    mask = t_vals <= tmax + 1e-10
    t_plot = t_vals[mask]
    Z = quantity[np.ix_(mask, mask)]

    vmax = args.vmax if args.vmax is not None else np.max(np.abs(Z))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: Re[-i Λ^<] = Im(Λ^<) — matches GEBK Fig. 3 top-left
    ax = axes[0]
    im = ax.pcolormesh(t_plot, t_plot, Z, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    fig.colorbar(im, ax=ax, label="Re[-i Λ$^<$]")
    ax.set_xlabel("t")
    ax.set_ylabel("t'")
    ax.set_title(args.title)
    ax.set_aspect("equal")

    # Right panel: diagonal t=t' (occupation n(t) = Im(G^<_imp(t,t')))
    ax2 = axes[1]
    nT_plot = len(t_plot)
    diag = quantity[np.ix_(mask, mask)].diagonal()
    ax2.plot(t_plot, diag, "b-", lw=1.5)
    ax2.set_xlabel("t")
    ax2.set_ylabel("Im(Λ$^<$(t,t)) = n_eff(t)")
    ax2.set_title("Diagonal (effective occupation)")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    outfile = args.output if args.output else args.lesser_file + ".png"
    fig.savefig(outfile, dpi=150)
    print(f"Wrote {outfile}")

    # Print summary statistics
    print(f"nT+1 = {nT1},  dt = {dt:.4f},  tmax = {t_vals[-1]:.4f}")
    print(f"Re[-i Λ^<] range: [{np.min(quantity):.4f}, {np.max(quantity):.4f}]")
    diag_full = quantity.diagonal()
    print(f"Diagonal (n_eff) range: [{np.min(diag_full):.4f}, {np.max(diag_full):.4f}]")

    # Sign check: diagonal must be positive (Cholesky requires -iΛ^< PSD)
    if np.any(diag_full < -1e-6):
        print("WARNING: negative diagonal — check sign convention or bath decomposition",
              file=sys.stderr)
        sys.exit(1)

    print("OK: diagonal is non-negative (consistent with positive-definite -iΛ^<_+)")


if __name__ == "__main__":
    main()
