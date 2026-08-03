#!/usr/bin/env python3
"""
2D (t,t') heatmap check of the near-atomic (W_i=0.1) Fig3L3 target after
the swapped-conjugate bug fix in NeqBathDecomposition.h -- confirms whether
the previously-diagnosed "input-field sensitivity, not a bug" rank-1
collapse (see project_gbek_cosine_bug.md) persists once the conjugate bug
(fixed 2026-07-08) is no longer a confound.

Usage:
    uv run --with numpy --with matplotlib python3 plot_fig3l3_post_fix.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, reconstruct
from provenance import write_provenance

WEISS_DELTA_FILE = "/Users/Shared/ornldev/code/dmrgpp/build/gebk-fig3-L3-weiss-delta-lesser"
PLUS_BATH_FILE = "/Users/Shared/ornldev/code/dmrgpp/build/gebk-fig3-L3-plus-bath-lesser"
OUT = "fig3L3_near_atomic_post_fix.png"

# Equilibrium bath fit result printed by the run (NumberOfBathPoints=5).
V_FIT = np.array([0.0053999167313065, -0.0160504072383437, -0.0052095584077009,
                  0.0053999167313065, -0.0160504072383437])
EPS_FIT = np.array([0.8959557182654330, 0.5321005311792999, 0.0,
                    -0.8959557182654330, -0.5321005311792999])
BETA = 16.0
DT = 0.04


def fermi(eps):
    x = BETA * eps
    return np.where(x > 500, 0.0, np.where(x < -500, 1.0, 1.0 / (1 + np.exp(x))))


def main():
    ts, re, im = read_lesser_file(WEISS_DELTA_FILE)
    raw_total = re + 1j * im
    N = len(ts)
    tn = np.arange(N) * DT
    tau = tn[:, None] - tn[None, :]

    raw_minus = np.zeros((N, N), dtype=complex)
    for a in range(len(V_FIT)):
        raw_minus += V_FIT[a] ** 2 * 1j * fermi(EPS_FIT[a]) * np.exp(-1j * EPS_FIT[a] * tau)

    lam_target = -1j * (raw_total - raw_minus)
    V_py = cholesky_causal(lam_target, L=3)
    recon_py = reconstruct(V_py)

    _, re2, im2 = read_lesser_file(PLUS_BATH_FILE)
    lam_cpp = re2 + 1j * im2

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    panels = [
        ("Target -i*Lambda^+_< (near-atomic, W_i=0.1)", lam_target),
        ("Python rank-3 cholesky_causal reconstruction", recon_py),
        ("cincuenta C++ recon (post conjugate-bug fix)", lam_cpp),
    ]
    vmax = max(np.max(np.abs(p[1].real)) for p in panels)
    for ax, (title, mat) in zip(axes, panels):
        w = sorted(np.linalg.eigvalsh(mat))[::-1]
        im_ = ax.imshow(mat.real, origin="lower", extent=[ts[0], ts[-1], ts[0], ts[-1]],
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
        ax.set_title(title + f"\ntop eigs: {w[0]:.1f}, {w[1]:.2f}, {w[2]:.2f}, {w[3]:.1e}",
                     fontsize=9)
        ax.set_xlabel("t'")
        ax.set_ylabel("t")
        fig.colorbar(im_, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("Re[-i*Lambda^+_<(t,t')]: near-atomic Fig3L3, L=3, post conjugate-bug fix",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"Wrote {OUT}")
    write_provenance(OUT, extra_files=[WEISS_DELTA_FILE, PLUS_BATH_FILE])

    print()
    print("target top eigs:    ", sorted(np.linalg.eigvalsh(lam_target))[-6:][::-1])
    print("python recon top eigs:", sorted(np.linalg.eigvalsh(recon_py))[-6:][::-1])
    print("cpp recon top eigs:   ", sorted(np.linalg.eigvalsh(lam_cpp))[-6:][::-1])


if __name__ == "__main__":
    main()
