#!/usr/bin/env python3
"""
2026-07-10: generalize scan_t3_activation.py's forced-single-column-timing
experiment to L=4, L=5, testing the hypothesis (user's intuition) that
later columns (4th, 5th, ...) have their OWN natural activation times,
later than column 3's t3~1.0-1.2, that a greedy sequential scan can locate
the same way: fix earlier columns' forced times at their already-found
optima, scan the next column's forced time, repeat.

Columns 0,1 always seed at rows 1,2 (standard exact recursion, never in
question). Columns 2..L-2 use whatever forced times are passed in (already
optimized from a smaller-L scan). The LAST column (index L-1) is the one
being scanned over `candidate_rows`. Rows between an activation and the
next stay at zero in not-yet-established columns (matching
cholesky_causal_deferred's k<L branch: plain forward substitution against
already-established columns only, no diagonal-constrained joint solve).
Rows after all L columns are established use the ordinary Eq. 62-63
optimal update.
"""
import numpy as np

from compare_reference import read_lesser_file
from gbek_cholesky import _solve_optimal_update, reconstruct

TARGET_L4 = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L4-weiss-delta-lesser"
TARGET_L5 = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L5-weiss-delta-lesser"


def load_lambda_from_delta(path):
    ts, re, im = read_lesser_file(path)
    return ts, -1j * (re + 1j * im)


def global_err(lam, A):
    return np.linalg.norm(lam - A) / np.linalg.norm(lam)


def eigenvector_lowrank(lam, L):
    w, U = np.linalg.eigh(lam)
    order = np.argsort(w)[::-1]
    w, U = w[order], U[:, order]
    return (U[:, :L] * w[:L]) @ U[:, :L].conj().T


def cholesky_forced_multi(lam, L, forced_rows):
    """forced_rows: length-L list, forced_rows[0]=1, forced_rows[1]=2 always;
    forced_rows[k] for k>=2 is the row at which column k activates (must be
    strictly increasing, forced_rows[k] >= k+1)."""
    assert forced_rows[0] == 1 and forced_rows[1] == 2
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

        # how many columns are established strictly before this row
        k = sum(1 for fr in forced_rows if fr < n)
        activating = (k < L and forced_rows[k] == n)

        if k < L and not activating:
            # forward substitution against the k established columns only
            for p in range(k):
                pivot_row_p = forced_rows[p]
                val = lam_at(n, pivot_row_p)
                for kk in range(p):
                    val -= V[n, kk] * np.conj(V[pivot_row_p, kk])
                pivot = V[pivot_row_p, p]
                V[n, p] = 0 if abs(pivot) < 1e-14 else val / np.conj(pivot)
        elif activating:
            for p in range(k):
                pivot_row_p = forced_rows[p]
                val = lam_at(n, pivot_row_p)
                for kk in range(p):
                    val -= V[n, kk] * np.conj(V[pivot_row_p, kk])
                pivot = V[pivot_row_p, p]
                V[n, p] = 0 if abs(pivot) < 1e-14 else val / np.conj(pivot)
            d = lam_at(n, n).real
            for kk in range(k):
                d -= abs(V[n, kk]) ** 2
            V[n, k] = np.sqrt(max(d, floor))
        else:
            a = np.array([lam_at(n, m) for m in range(1, n)])
            Q = V[1:n, :]
            d = lam_at(n, n).real
            V[n, :] = _solve_optimal_update(Q, a, d, L)
    return V


def greedy_scan(lam, L, fixed_rows, candidate_rows, dt):
    """fixed_rows: forced_rows[0:len(fixed_rows)] already chosen (includes
    the mandatory [1,2] prefix). Scans the NEXT column over candidate_rows."""
    results = []
    for r in candidate_rows:
        forced = fixed_rows + [r]
        V = cholesky_forced_multi(lam, len(forced), forced)
        ch = reconstruct(V)
        e = global_err(lam, ch)
        results.append((r, e))
    return results


if __name__ == "__main__":
    for label, target, L in [("L=4", TARGET_L4, 4), ("L=5", TARGET_L5, 5)]:
        print(f"\n=== {label} ===")
        ts, lam = load_lambda_from_delta(target)
        dt = ts[1] - ts[0]
        ev_err = global_err(lam, eigenvector_lowrank(lam, L))

        # baseline: plain time-ordered seeding (forced_rows = [1,2,3,4,...])
        baseline_rows = list(range(1, L + 1))
        V_base = cholesky_forced_multi(lam, L, baseline_rows)
        base_err = global_err(lam, reconstruct(V_base))
        print(f"baseline (plain seeding): err[ch]={base_err:.4f}  err[ev]={ev_err:.4f}  "
              f"ratio={base_err/ev_err:.2f}x")

        # Step 1: scan t3 (column index 2), holding columns 0,1 fixed, all
        # later columns (if L>3) still at plain time order for this step.
        fixed = [1, 2]
        candidates_t3 = [3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 26, 30, 35, 40]
        # for L>3, pad remaining columns at their plain-seeding rows so the
        # matrix stays L-dimensional during this first scan
        best_t3, best_t3_err = None, np.inf
        for r in candidates_t3:
            trial = fixed + [r] + list(range(r + 1, r + 1 + (L - 3)))
            if trial[-1] > lam.shape[0] - 1:
                continue
            V = cholesky_forced_multi(lam, L, trial)
            e = global_err(lam, reconstruct(V))
            if e < best_t3_err:
                best_t3_err, best_t3 = e, r
        print(f"best t3: row={best_t3} (t={best_t3*dt:.2f})  err[ch]={best_t3_err:.4f}")

        if L == 3:
            continue

        # Step 2: fix column 2 at best_t3, scan column 3 (t4)
        fixed2 = [1, 2, best_t3]
        candidates_t4 = [r for r in [best_t3 + 1, best_t3 + 2, best_t3 + 4,
                                       best_t3 + 6, best_t3 + 8, best_t3 + 10,
                                       best_t3 + 14, best_t3 + 18, best_t3 + 24,
                                       best_t3 + 30, best_t3 + 36]
                         if r <= lam.shape[0] - 1 - (L - 4)]
        best_t4, best_t4_err = None, np.inf
        print(f"{'t4_row':>7} {'t4':>6} {'err[ch]':>10}")
        for r in candidates_t4:
            trial = fixed2 + [r] + list(range(r + 1, r + 1 + (L - 4)))
            V = cholesky_forced_multi(lam, L, trial)
            e = global_err(lam, reconstruct(V))
            print(f"{r:7d} {r*dt:6.2f} {e:10.4f}")
            if e < best_t4_err:
                best_t4_err, best_t4 = e, r
        print(f"best t4: row={best_t4} (t={best_t4*dt:.2f})  err[ch]={best_t4_err:.4f}  "
              f"ratio={best_t4_err/ev_err:.2f}x")

        if L == 4:
            continue

        # Step 3 (L=5): fix columns 2,3 at best_t3,best_t4, scan column 4 (t5)
        fixed3 = [1, 2, best_t3, best_t4]
        candidates_t5 = [r for r in [best_t4 + 1, best_t4 + 2, best_t4 + 4,
                                       best_t4 + 6, best_t4 + 8, best_t4 + 10,
                                       best_t4 + 14, best_t4 + 18, best_t4 + 24,
                                       best_t4 + 30, best_t4 + 36]
                         if r <= lam.shape[0] - 1]
        best_t5, best_t5_err = None, np.inf
        print(f"{'t5_row':>7} {'t5':>6} {'err[ch]':>10}")
        for r in candidates_t5:
            trial = fixed3 + [r]
            V = cholesky_forced_multi(lam, L, trial)
            e = global_err(lam, reconstruct(V))
            print(f"{r:7d} {r*dt:6.2f} {e:10.4f}")
            if e < best_t5_err:
                best_t5_err, best_t5 = e, r
        print(f"best t5: row={best_t5} (t={best_t5*dt:.2f})  err[ch]={best_t5_err:.4f}  "
              f"ratio={best_t5_err/ev_err:.2f}x")
