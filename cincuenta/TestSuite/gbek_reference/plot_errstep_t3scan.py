#!/usr/bin/env python3
"""
2026-07-10: err^step(t) (GBEK Eq. 67, Fig. 3 bottom-left-panel style) for the
best t3-activation found in scan_t3_activation.py (t3=1.04, n3=26, err[ch]
global=0.28), overlaid against the current baseline (plain time-ordered
seeding, t3=0.12, err[ch]=0.60) and the eigenvector low-rank reference, all
on the same L=3 atomic-limit target used for that scan.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, reconstruct
from scan_t3_activation import cholesky_forced_t3, global_err
from provenance import write_provenance

TARGET = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-weiss-delta-lesser"
L = 3
BEST_N3 = 26  # t3=1.04, the scan's minimum


def load_lambda_from_delta(path):
    ts, re, im = read_lesser_file(path)
    return ts, -1j * (re + 1j * im)


def eigenvector_lowrank(lam, L):
    w, U = np.linalg.eigh(lam)
    order = np.argsort(w)[::-1]
    w, U = w[order], U[:, order]
    return (U[:, :L] * w[:L]) @ U[:, :L].conj().T


def errstep_curve(lam, A):
    N = lam.shape[0] - 1
    out = np.full(N + 1, np.nan)
    for Nidx in range(1, N + 1):
        diff = lam[1:Nidx + 1, Nidx] - A[1:Nidx + 1, Nidx]
        w = np.full(Nidx, 2.0)
        w[-1] = 1.0
        out[Nidx] = np.sqrt(np.sum(w * np.abs(diff) ** 2))
    return out


if __name__ == "__main__":
    ts, lam = load_lambda_from_delta(TARGET)
    dt = ts[1] - ts[0]

    V_baseline = cholesky_causal(lam, L)
    ch_baseline = reconstruct(V_baseline)

    V_best = cholesky_forced_t3(lam, L, BEST_N3)
    ch_best = reconstruct(V_best)

    ev = eigenvector_lowrank(lam, L)

    curve_baseline = errstep_curve(lam, ch_baseline)
    curve_best = errstep_curve(lam, ch_best)
    curve_ev = errstep_curve(lam, ev)

    print(f"{'curve':<32} {'err^step(t=4)':>14} {'global err[A]':>14}")
    print(f"{'baseline (t3=0.12)':<32} {curve_baseline[-1]:>14.4f} {global_err(lam, ch_baseline):>14.4f}")
    print(f"{'t3=1.04 (scan minimum)':<32} {curve_best[-1]:>14.4f} {global_err(lam, ch_best):>14.4f}")
    print(f"{'eigenvector':<32} {curve_ev[-1]:>14.4f} {global_err(lam, ev):>14.4f}")

    fig, ax = plt.subplots(figsize=(7, 5.8))
    ax.semilogy(ts, curve_baseline, '-', color='crimson', lw=2.0,
                label=r"$(\Lambda_+^<)^{ch}$, plain seeding (t3=0.12)")
    ax.semilogy(ts, curve_best, '-', color='darkorange', lw=2.2,
                label=r"$(\Lambda_+^<)^{ch}$, forced t3=1.04 (scan minimum)")
    ax.semilogy(ts, curve_ev, '--', color='seagreen', lw=1.6, alpha=0.8,
                label=r"$(\Lambda_+^<)^{ev}$ (eigenvector, acausal)")
    ax.axvline(1.04, color='darkorange', ls=':', lw=1, alpha=0.6)
    ax.set_xlim(0, ts[-1])
    ax.set_ylim(1e-6, max(3, ax.get_ylim()[1]))
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\mathrm{err}^{\mathrm{step}}(t)$  (Eq. 67)")
    ax.set_title("cf. GBEK Fig. 3 bottom-left panel -- t3-activation comparison, L=3", fontsize=10)
    ax.legend(fontsize=8.5, loc='lower right')
    ax.grid(alpha=0.2, which='both')
    fig.tight_layout()
    fig.savefig("errstep_t3scan_comparison.png", dpi=150)
    prov = write_provenance("errstep_t3scan_comparison.png", extra_files=[TARGET],
                             notes=f"baseline vs forced t3={BEST_N3*dt:.2f} (n3={BEST_N3}) vs eigenvector")
    print(f"\nWrote errstep_t3scan_comparison.png, {prov}")
