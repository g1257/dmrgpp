#!/usr/bin/env python3
"""
2026-07-09: Fable's proposed decisive test for the seeding-timing hypothesis
-- scan the FIXED activation time of the third causal-Cholesky column (L=3)
over rows 1..50 (t in [0.04, 2.0]) and plot final err[ch] (global relative
Frobenius error over the full t_max=4 window), looking for an interior
minimum well past row 3 (t=0.12, where strict time-ordered seeding
currently establishes it).

Distinct from the already-tried cholesky_causal_deferred(rtol=1e-6): that
threshold, measured against the running-max diagonal, fires column 3
within a few rows of t=0.12 regardless of rtol (verified: our residual
chain 3.6e-3, 1.4e-5, 6.5e-6, 1.4e-7, 2.1e-8 clears 1e-6 already at row 3),
so it never actually varied activation time the way this scan does.

Columns 0,1 (rows 1,2) always seed via the standard exact recursion.
Column 2 is forced to activate at a specified row n3 >= 3: rows 3..n3-1
are computed via forward substitution using only the 2 established
columns (matching cholesky_causal_deferred's k<L branch -- NOT the
Eq. 63 diagonal-constrained optimal update, since with fewer than L
established columns there's no joint objective to solve, just ordinary
forward substitution against known pivots); row n3 seeds column 2 via the
standard sqrt(residual) formula; rows > n3 use the ordinary Eq. 62-63
optimal update with all 3 columns.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import _solve_optimal_update, reconstruct
from provenance import write_provenance

TARGET = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-weiss-delta-lesser"
L = 3


def load_lambda_from_delta(path):
    ts, re, im = read_lesser_file(path)
    return ts, -1j * (re + 1j * im)


def global_err(lam, A):
    return np.linalg.norm(lam - A) / np.linalg.norm(lam)


def cholesky_forced_t3(lam, L, n3):
    """L must be 3 for this scan (columns 0,1 always seed at rows 1,2;
    column 2 forced to activate at row n3)."""
    assert L == 3
    nT = lam.shape[0] - 1
    V = np.zeros((nT + 1, L), dtype=complex)

    def lam_at(n, j):
        if j <= n:
            return lam[n, j]
        return np.conj(lam[j, n])

    max_diag_seen = 0.0
    for n in range(1, nT + 1):
        max_diag_seen = max(max_diag_seen, abs(lam_at(n, n).real))
        floor = 1e-10 * max(max_diag_seen, 1e-300)

        if n < n3:
            # k=2 established columns (0,1) if n>=3, else standard seeding.
            k = min(n - 1, 2)
            for p in range(k):
                val = lam_at(n, p + 1)
                for kk in range(p):
                    val -= V[n, kk] * np.conj(V[p + 1, kk])
                pivot = V[p + 1, p]
                V[n, p] = 0 if abs(pivot) < 1e-14 else val / np.conj(pivot)
            if n <= 2:
                d = lam_at(n, n).real
                for kk in range(k):
                    d -= abs(V[n, kk]) ** 2
                V[n, k] = np.sqrt(max(d, floor))
            # else: column 2 stays 0 (deferred), no diagonal fit attempted.
        elif n == n3:
            for p in range(2):
                val = lam_at(n, p + 1)
                for kk in range(p):
                    val -= V[n, kk] * np.conj(V[p + 1, kk])
                pivot = V[p + 1, p]
                V[n, p] = 0 if abs(pivot) < 1e-14 else val / np.conj(pivot)
            d = lam_at(n, n).real
            for kk in range(2):
                d -= abs(V[n, kk]) ** 2
            V[n, 2] = np.sqrt(max(d, floor))
        else:
            a = np.array([lam_at(n, m) for m in range(1, n)])
            Q = V[1:n, :]
            d = lam_at(n, n).real
            V[n, :] = _solve_optimal_update(Q, a, d, L)
    return V


if __name__ == "__main__":
    ts, lam = load_lambda_from_delta(TARGET)
    dt = ts[1] - ts[0]

    rows = [3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50]
    errs = []
    print(f"{'n3':>4} {'t3':>6} {'err[ch]':>10}")
    for n3 in rows:
        V = cholesky_forced_t3(lam, L, n3)
        ch = reconstruct(V)
        e = global_err(lam, ch)
        errs.append(e)
        print(f"{n3:4d} {n3*dt:6.2f} {e:10.4f}")

    baseline = errs[0]  # n3=3 == plain time-ordered seeding
    print(f"\nbaseline (n3=3, plain seeding): err[ch]={baseline:.4f}")
    print(f"min over scan: err[ch]={min(errs):.4f} at n3={rows[int(np.argmin(errs))]} "
          f"(t3={rows[int(np.argmin(errs))]*dt:.2f})")

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot([r * dt for r in rows], errs, 'o-')
    ax.axhline(baseline, color='gray', ls='--', lw=1, label=f'baseline (t3=0.12): {baseline:.3f}')
    ax.axhline(0.17, color='red', ls=':', lw=1, label='paper err[ch]=0.17')
    ax.set_xlabel('forced activation time t3 of column 3')
    ax.set_ylabel('err[ch] (global, full t_max=4 window)')
    ax.set_title('Fable t3-activation scan, L=3 atomic-limit target')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('t3_activation_scan.png', dpi=150)
    prov = write_provenance('t3_activation_scan.png', extra_files=[TARGET],
                             notes=f"rows scanned={rows}")
    print(f"\nWrote t3_activation_scan.png, {prov}")
