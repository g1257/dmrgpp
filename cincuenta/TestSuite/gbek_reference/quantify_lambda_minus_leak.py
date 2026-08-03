#!/usr/bin/env python3
"""
Compute Lambda^- = Lambda - Lambda^+_dumped directly from cincuenta's own
dumped output, for a single run.

CORRECTED HISTORICAL NOTE (see README.md's "A wrong turn" section): this
script was originally written to test the hypothesis that cincuenta's
first bath (Lambda^-) was leaking the persistent steady-state field that
should belong to the second bath (Lambda^+), causing Lambda^+ to decay too
fast. That hypothesis is WRONG. Directly computing the analytic Lambda^-
formula from a run's own fitted equilibrium bath parameters (see
check_cholesky_step.py) shows Lambda^- is tiny (~3e-4) and constant for the
whole trajectory -- nowhere near large enough to explain a missing ~0.5.

What this script actually computes, "Lambda - Lambda^+_dumped", is NOT the
literal physical Lambda^- whenever cincuenta's Cholesky reconstruction of
Lambda^+ is itself inaccurate (which, per the real bug documented in
README.md and gbek_cholesky.py, it was) -- it's just "whatever the
reconstruction fails to capture", which is a different and potentially
misleading quantity to label "Lambda^-". Kept in the repo as a record of
this wrong turn and as a general "compute Lambda^- from Lambda and Lambda^+
dumps" utility, but don't trust its original interpretation, and prefer
check_cholesky_step.py (which uses the literal analytic Lambda^- formula)
for actually diagnosing a Lambda^+ discrepancy.

Usage:
    python3 quantify_lambda_minus_leak.py <prefix>
    # e.g. prefix = /path/to/gebk-fig3-L3   (looks for
    #   <prefix>-weiss-delta-lesser and <prefix>-plus-bath-lesser)
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("prefix", help="cincuenta NeqOutputPrefix path, e.g. .../gebk-fig3-L3")
    ap.add_argument("--out", default="lambda_minus_leak_check.png")
    args = ap.parse_args()

    ts, raw_re, raw_im = read_lesser_file(f"{args.prefix}-weiss-delta-lesser")
    ts2, plus_re, plus_im = read_lesser_file(f"{args.prefix}-plus-bath-lesser")
    assert np.allclose(ts, ts2), "weiss-delta-lesser and plus-bath-lesser grids differ"

    # weiss-delta-lesser stores the raw lesser component itself (KadanoffBaym
    # convention: Re, Im columns are its real/imag parts). plus-bath-lesser
    # stores Lambda^+_< directly (dumpPlusBath's own convention). Put both on
    # the same footing as -i*X: -i*(a+ib) = b - ia.
    lam_total = raw_im - 1j * raw_re       # Lambda(t,t')
    lam_plus = plus_re + 1j * plus_im      # already Lambda^+_<(t,t') as dumped
    lam_minus = lam_total - lam_plus       # Lambda^-_<(t,t'), by linearity

    diag_total = np.diag(lam_total).real
    diag_plus = np.diag(lam_plus).real
    diag_minus = np.diag(lam_minus).real

    print(f"{'t':>6} {'Lambda<(t,t)':>14} {'Lambda+<(t,t)':>14} {'Lambda-<(t,t)':>14}")
    for n in range(0, len(ts), max(1, len(ts) // 10)):
        print(f"{ts[n]:6.2f} {diag_total[n]:14.5f} {diag_plus[n]:14.5f} {diag_minus[n]:14.5f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(ts, diag_total, label="Total: Lambda<(t,t)", lw=2)
    axes[0].plot(ts, diag_plus, "--", label="Lambda+<(t,t)  (second bath)", lw=2)
    axes[0].plot(ts, diag_minus, ":", label="Lambda-<(t,t)  (first bath)", lw=2)
    axes[0].axhline(0.5, color="gray", lw=0.5)
    axes[0].set_xlabel("t")
    axes[0].set_title("Does Lambda^- carry the field that should be in Lambda^+?")
    axes[0].legend(fontsize=8)

    im = axes[1].imshow(np.real(lam_minus), origin="lower", extent=[ts[0], ts[-1]] * 2,
                         cmap="RdBu_r")
    axes[1].set_title("Re[Lambda^-_<(t,t')]")
    plt.colorbar(im, ax=axes[1])

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
