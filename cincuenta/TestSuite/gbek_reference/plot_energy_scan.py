#!/usr/bin/env python3
"""
Plot GBEK PRB 88, 235106 (2013) Fig. 9 ("test of energy conservation") from
the .npz files produced by run_fig9_scan.py, and print the actual
quantitative claim the figure makes: Etot(t) should be FLAT for t > tq
(energy conservation after the ramp), with the residual non-flatness (a
finite-Lbath Cholesky truncation artifact) shrinking as Lbath grows 4->6->8.

Layout: single panel (unlike Fig. 7/8's eigenvector/Cholesky two-panel
layout -- Fig. 9 is Cholesky only), color = quantity (Ekin green, Eint blue,
Etot red), linestyle = U (solid U=2, dashed U=4), linewidth = Lbath (thin=4,
medium=6, thick=8) -- matching the paper's own encoding.

Note: unlike plot_docc_scan.py's Fig. 7/8 "--ylim paper" preset (whose exact
pixel-read axis ranges were independently confirmed against page_12.png),
Fig. 9's exact axis range has not been independently pixel-checked here --
use --ylim explicitly (or the default autoscale) rather than trusting an
unverified preset.

Usage:
    uv run --with numpy --with matplotlib python3 plot_energy_scan.py
    uv run --with numpy --with matplotlib python3 plot_energy_scan.py --ylim -1.6,0.6
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from provenance import write_provenance

U_VALUES = (2.0, 4.0)  # the paper's actual Fig. 9 curves
L_VALUES = (2, 3, 4)   # Lbath = 4, 6, 8
QUANTITIES = [("Ekin_t", "Ekin", "green"), ("Eint_t", "Eint", "blue"), ("Etot_t", "Etot", "red")]
U_STYLES = {2.0: "-", 4.0: "--"}
L_WIDTHS = {2: 1.2, 3: 2.0, 4: 2.8}  # Lbath=4,6,8 -> thin/medium/thick


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="fig9_energy.png")
    ap.add_argument("--title", default="cf. GBEK Fig. 9: test of energy conservation")
    ap.add_argument("--ylim", default=None,
                    help="'ymin,ymax', or omit to autoscale (see module docstring)")
    args = ap.parse_args()
    ylim = tuple(float(x) for x in args.ylim.split(",")) if args.ylim else None

    data = {}
    for U in U_VALUES:
        for L in L_VALUES:
            fname = f"fig9_energy_U{U:g}_L{L}.npz"
            data[(U, L)] = np.load(fname)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    tq = None
    for U in U_VALUES:
        for L in L_VALUES:
            d = data[(U, L)]
            ts = d["ts"]
            tq = float(d["tq"])
            for key, label, color in QUANTITIES:
                ax.plot(ts, d[key], color=color, linestyle=U_STYLES[U],
                        linewidth=L_WIDTHS[L],
                        label=f"{label} U={U:g} Lbath={int(d['Lbath'])}")
    if tq is not None:
        ax.axvline(tq, color="gray", linestyle=":", linewidth=1.0)
    ax.set_xlabel("t")
    ax.set_ylabel("energy (units of v0)")
    ax.set_title(args.title)
    ax.grid(alpha=0.2)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(fontsize=6, ncol=2, loc="lower right")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")

    extra_files = ["plot_energy_scan.py"] + \
        [f"fig9_energy_U{U:g}_L{L}.npz" for U in U_VALUES for L in L_VALUES]
    prov = write_provenance(args.out, extra_files=extra_files, notes="figure=9")
    print(f"Wrote {prov}")

    print("\n--- Fig. 9 physics criteria ---")
    print("Claim: Etot(t) should be flat for t > tq (energy conservation), with "
          "residual non-flatness (Cholesky truncation artifact) SHRINKING as Lbath "
          "grows 4->6->8. Reporting actual numbers, not smoothing over a violation.")
    for U in U_VALUES:
        print(f"\nU={U:g}:")
        flatness = {}
        for L in L_VALUES:
            d = data[(U, L)]
            ts, Etot_t, tq_ = d["ts"], d["Etot_t"], float(d["tq"])
            mask = ts > tq_
            if not mask.any():
                continue
            etot_at_tq = np.interp(tq_, ts, Etot_t)
            flat = np.max(np.abs(Etot_t[mask] - etot_at_tq))
            flatness[L] = flat
            print(f"  Lbath={2*L}: max|Etot(t)-Etot(tq)| for t>tq = {flat:.4e}")
        Ls = sorted(flatness)
        shrinking = all(flatness[Ls[i]] >= flatness[Ls[i + 1]] for i in range(len(Ls) - 1))
        print(f"  shrinks monotonically with Lbath: {shrinking}"
              + ("" if shrinking else "  <-- contradicts the paper's claim, investigate"))
        max_im = max(float(data[(U, L)]["max_im_ekin"]) for L in L_VALUES)
        print(f"  max|Im(Ekin)| across Lbath: {max_im:.3e} (expect ~0)")


if __name__ == "__main__":
    main()
