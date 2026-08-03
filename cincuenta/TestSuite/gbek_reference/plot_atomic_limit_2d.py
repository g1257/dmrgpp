#!/usr/bin/env python3
"""
2D (t,t') heatmap comparison for the true-atomic-limit GBEK run: target
Re[-i*Lambda^+_<(t,t')] vs. the independently-validated Python rank-3
Cholesky reconstruction (gbek_cholesky.py::cholesky_causal) vs. cincuenta's
own C++ online reconstruction (dumped plus-bath-lesser).

Purpose: the 2026-07-07 overnight report claimed the Python reconstruction
does NOT collapse to rank 1 on this target, based on eigenvalue spectra and
a diagonal-only slice. This script shows the full (t,t') plane so that
claim can be checked visually rather than taken on an eigenvalue printout
alone -- per this session's standing rule to show 2D structure, not just a
1D diagonal/off-diagonal slice, when claiming or refuting rank collapse.

Usage:
    uv run --with numpy --with matplotlib python3 plot_atomic_limit_2d.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, reconstruct
from provenance import write_provenance

TARGET_FILE = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-weiss-delta-lesser"
CPP_FILE = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-plus-bath-lesser"
OUT = "atomic_limit_2d_rank_comparison.png"


def main():
    ts, re, im = read_lesser_file(TARGET_FILE)
    lam_target = -1j * (re + 1j * im)  # Lambda; Lambda^-=0 here, so this IS Lambda^+

    V = cholesky_causal(lam_target, L=3)
    lam_python = reconstruct(V)

    _, re2, im2 = read_lesser_file(CPP_FILE)
    lam_cpp = re2 + 1j * im2  # plus-bath-lesser is already dumped in Lambda units

    w_target = np.sort(np.linalg.eigvalsh(lam_target))[::-1]
    w_python = np.sort(np.linalg.eigvalsh(lam_python))[::-1]
    w_cpp = np.sort(np.linalg.eigvalsh(lam_cpp))[::-1]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    panels = [
        ("Target -i*Lambda^+_< (Lambda^-=0 exactly)", lam_target, w_target),
        ("Python rank-3 cholesky_causal reconstruction", lam_python, w_python),
        ("cincuenta C++ online reconstruction (plus-bath-lesser)", lam_cpp, w_cpp),
    ]
    vmax = max(np.max(np.abs(p[1].real)) for p in panels)
    for ax, (title, mat, eigs) in zip(axes, panels):
        im_ = ax.imshow(mat.real, origin="lower", extent=[ts[0], ts[-1], ts[0], ts[-1]],
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
        ax.set_title(title + f"\ntop eigs: {eigs[0]:.1f}, {eigs[1]:.2f}, {eigs[2]:.2f}, {eigs[3]:.1e}",
                     fontsize=9)
        ax.set_xlabel("t'")
        ax.set_ylabel("t")
        fig.colorbar(im_, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("Re[-i*Lambda^+_<(t,t')]: true atomic limit, L=3", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"Wrote {OUT}")

    write_provenance(OUT, extra_files=[TARGET_FILE, CPP_FILE])

    print()
    print("Eigenvalue spectra (top 6, descending):")
    print("  target: ", np.round(w_target[:6], 4))
    print("  python: ", np.round(w_python[:6], 6))
    print("  cpp:    ", np.round(w_cpp[:6], 6))
    print()
    print("max|python_recon - target| (rank-3 residual):", np.max(np.abs(lam_python - lam_target)))
    print("max|cpp_recon - target|:                     ", np.max(np.abs(lam_cpp - lam_target)))


if __name__ == "__main__":
    main()
