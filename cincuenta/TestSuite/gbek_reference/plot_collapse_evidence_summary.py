#!/usr/bin/env python3
"""
Summary evidence figure for the near-atomic (W_i=0.1) causal-Cholesky rank
collapse: eigenvalue spectra (bar chart) and diagonal-vs-t traces for BOTH
the working atomic-limit case and the collapsing near-atomic case, side by
side, so the contrast is visible in one figure.

Usage:
    uv run --with numpy --with matplotlib python3 plot_collapse_evidence_summary.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, reconstruct
from provenance import write_provenance

OUT = "collapse_evidence_summary.png"

# --- Atomic limit (working) ---
AL_TARGET = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-weiss-delta-lesser"
AL_CPP = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-plus-bath-lesser"

# --- Near-atomic (collapsing) ---
NA_TOTAL = "/Users/Shared/ornldev/code/dmrgpp/build/gebk-fig3-L3-weiss-delta-lesser"
NA_CPP = "/Users/Shared/ornldev/code/dmrgpp/build/gebk-fig3-L3-plus-bath-lesser"
V_FIT = np.array([0.0053999167313065, -0.0160504072383437, -0.0052095584077009,
                  0.0053999167313065, -0.0160504072383437])
EPS_FIT = np.array([0.8959557182654330, 0.5321005311792999, 0.0,
                    -0.8959557182654330, -0.5321005311792999])
BETA = 16.0
DT_NA = 0.04


def fermi(eps):
    x = BETA * eps
    return np.where(x > 500, 0.0, np.where(x < -500, 1.0, 1.0 / (1 + np.exp(x))))


def load_atomic_limit():
    ts, re, im = read_lesser_file(AL_TARGET)
    lam_target = -1j * (re + 1j * im)
    _, re2, im2 = read_lesser_file(AL_CPP)
    lam_cpp = re2 + 1j * im2
    return ts, lam_target, lam_cpp


def load_near_atomic():
    ts, re, im = read_lesser_file(NA_TOTAL)
    raw_total = re + 1j * im
    N = len(ts)
    tn = np.arange(N) * DT_NA
    tau = tn[:, None] - tn[None, :]
    raw_minus = np.zeros((N, N), dtype=complex)
    for a in range(len(V_FIT)):
        raw_minus += V_FIT[a] ** 2 * 1j * fermi(EPS_FIT[a]) * np.exp(-1j * EPS_FIT[a] * tau)
    lam_target = -1j * (raw_total - raw_minus)
    _, re2, im2 = read_lesser_file(NA_CPP)
    lam_cpp = re2 + 1j * im2
    return ts, lam_target, lam_cpp


def main():
    ts_al, target_al, cpp_al = load_atomic_limit()
    py_al = reconstruct(cholesky_causal(target_al, L=3))

    ts_na, target_na, cpp_na = load_near_atomic()
    py_na = reconstruct(cholesky_causal(target_na, L=3))

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # --- Eigenvalue spectra ---
    for ax, (title, target, cpp, py) in zip(
        axes[0],
        [
            ("Atomic limit (Lambda^-===0 exactly): WORKS", target_al, cpp_al, py_al),
            ("Near-atomic (W_i=0.1 approx): COLLAPSES", target_na, cpp_na, py_na),
        ],
    ):
        w_t = sorted(np.linalg.eigvalsh(target))[::-1][:6]
        w_c = sorted(np.linalg.eigvalsh(cpp))[::-1][:6]
        w_p = sorted(np.linalg.eigvalsh(py))[::-1][:6]
        x = np.arange(6)
        width = 0.27
        ax.bar(x - width, w_t, width, label="target", color="#333333")
        ax.bar(x, w_c, width, label="cincuenta C++", color="#d62728")
        ax.bar(x + width, w_p, width, label="Python (independent)", color="#1f77b4")
        ax.set_yscale("symlog", linthresh=1e-3)
        ax.set_xticks(x)
        ax.set_xticklabels([f"$\\lambda_{{{i}}}$" for i in range(6)])
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("eigenvalue")
        ax.legend(fontsize=8)

    # --- Diagonal vs t ---
    for ax, (title, ts, target, cpp, py) in zip(
        axes[1],
        [
            ("Atomic limit: reconstruction tracks target's ridge", ts_al, target_al, cpp_al, py_al),
            ("Near-atomic: reconstruction decays away from target's plateau",
             ts_na, target_na, cpp_na, py_na),
        ],
    ):
        diag_t = np.diag(target).real
        diag_c = np.diag(cpp).real
        diag_p = np.diag(py).real
        ax.plot(ts, diag_t, color="#333333", lw=2, label="target")
        ax.plot(ts, diag_c, color="#d62728", lw=1.5, ls="--", label="cincuenta C++")
        ax.plot(ts, diag_p, color="#1f77b4", lw=1.5, ls=":", label="Python (independent)")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("t = t'")
        ax.set_ylabel("Re[-i*Lambda^+_<(t,t)]")
        ax.legend(fontsize=8)

    fig.suptitle("Rank-3 Cholesky reconstruction: atomic limit (works) vs. near-atomic (collapses)",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"Wrote {OUT}")
    write_provenance(OUT, extra_files=[AL_TARGET, AL_CPP, NA_TOTAL, NA_CPP])


if __name__ == "__main__":
    main()
