#!/usr/bin/env python3
"""
Isolate whether cincuenta's second-bath failure (Lambda^+ collapsing to ~0
well before the paper's expected L=3 accuracy horizon) is in the Cholesky
DECOMPOSITION step or the real-time PROPAGATION step.

Method: read cincuenta's own dumped total Lambda (weiss-delta-lesser) for a
run, subtract the analytic first-bath Lambda^- (computed directly from the
run's own fitted equilibrium bath parameters -- confirmed tiny, ~3e-4, via
check_analytic_deltaminus.py), and feed the result through our own
independently-validated cholesky_causal() (a faithful port of
NeqBathDecomposition.h). If the resulting reconstruction matches cincuenta's
actual plus-bath-lesser dump, the C++ Cholesky step is behaving exactly as
the algorithm dictates and the failure is a genuine property of decomposing
this specific target at rank L (not a coding bug). If it diverges from
cincuenta's own dump, there's a real discrepancy to chase down separately
from rank-truncation.

Usage:
    python3 check_cholesky_step.py <weiss-delta-lesser> <plus-bath-lesser> \\
        --V <v0,v1,...> --eps <e0,e1,...> --beta <beta> --L <rank> --dt <dt>
"""
import argparse
import numpy as np

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal


def fermi(eps, beta):
    x = beta * eps
    return np.where(x > 500, 0.0, np.where(x < -500, 1.0, 1.0 / (1 + np.exp(x))))


def analytic_lambda_minus_lesser(ts, V, eps, beta):
    """Lambda^-_<(t_n, t_j) for all n,j, from the free-propagation formula
    in NeqBathDecomposition.h::deltaMinusLesser (there without the leading
    -i; here folded in directly since that's the convention this script
    works in throughout)."""
    N = len(ts)
    f = fermi(eps, beta)
    out = np.zeros((N, N), dtype=complex)
    for n in range(N):
        for j in range(N):
            dt_nj = ts[n] - ts[j]
            out[n, j] = np.sum(V**2 * f * np.exp(-1j * eps * dt_nj))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("delta_file", help="cincuenta's weiss-delta-lesser dump (total Lambda)")
    ap.add_argument("plus_file", help="cincuenta's plus-bath-lesser dump (its own Lambda^+)")
    ap.add_argument("--V", required=True, help="comma-separated fitted hoppings")
    ap.add_argument("--eps", required=True, help="comma-separated fitted bath energies")
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--L", type=int, required=True)
    args = ap.parse_args()

    V = np.array([float(x) for x in args.V.split(",")])
    eps = np.array([float(x) for x in args.eps.split(",")])

    ts, raw_re, raw_im = read_lesser_file(args.delta_file)
    ts2, plus_re, plus_im = read_lesser_file(args.plus_file)
    assert np.allclose(ts, ts2)
    dt = ts[1] - ts[0]

    raw_less = raw_re + 1j * raw_im                  # the dump's raw lesser component
    lam_total = -1j * raw_less                       # Lambda(t,t')
    lam_minus = analytic_lambda_minus_lesser(ts, V, eps, args.beta)  # Lambda^-
    lam_plus_target = lam_total - lam_minus          # Lambda^+ : the Cholesky's target

    print("max|lam_minus| over full grid:", np.max(np.abs(lam_minus)),
          "(should be tiny -- confirms Lambda^- is negligible)")

    # cholesky_causal expects the lower triangle (j<=n) of its input filled;
    # it internally Hermitian-completes the rest via lam_at.
    lam_input = np.where(np.tril(np.ones_like(lam_plus_target, dtype=bool)),
                          lam_plus_target, 0)
    V_ours = cholesky_causal(lam_input, args.L)
    lam_plus_ours = V_ours @ V_ours.conj().T

    diag_cincuenta = plus_re.diagonal()
    diag_ours = lam_plus_ours.real.diagonal()
    diag_target = lam_plus_target.real.diagonal()

    print(f"\n{'t':>6} {'target (should be ~0.5)':>24} {'our Cholesky':>14} {'cincuenta dump':>16}")
    for n in range(0, len(ts), max(1, len(ts) // 20)):
        print(f"{ts[n]:6.2f} {diag_target[n]:24.5f} {diag_ours[n]:14.5f} {diag_cincuenta[n]:16.5f}")

    err = np.max(np.abs(diag_ours - diag_cincuenta))
    print(f"\nmax|our_reconstruction - cincuenta_dump| on diagonal: {err:.6f}")


if __name__ == "__main__":
    main()
