#!/usr/bin/env python3
"""
Reproduce GBEK Fig. 4: time evolution of the rank-L "second bath"
hybridization amplitudes V+_0p(t), for both decomposition methods, plus the
eigenvalue-decay inset.

Paper text (Sec. VI.A, around Eq. 62-66):
    "FIG. 4. Time evolution of the hybridization V+_0p(t). Both panels show
    a rank L=3 approximation that corresponds to the input Weiss field
    -iLambda^<_+ displayed in Fig. 3. To obtain V^ch_0p(t), in the top
    panel, we used the low-rank Cholesky approach. In the lower panel, we
    plot the hybridization V^ev_0p(t), calculated using the low-rank
    eigenvector approximation. The inset shows the decay of the eigenvalues
    a_p of -iLambda^<_+."
    "Because of Lambda_+(0,t) = Lambda_+(t,0) = 0, we obtain V^ch_0p(0) = 0."

Target: gbek-atomic-limit-exact-lesser only (Python's own independent
atomic-limit self-consistency output -- the one that satisfies the paper's
own exact physical constraint -iLambda^<_+(t,t) = 1/2 for t > t_q under
Lambda_- = 0, Lambda = Lambda_+, Sec. VI.A). Do NOT point this at a
cincuenta C++ dump -- see project notes on the C++ target's own diagonal
drift (open, separate issue, not yet resolved).

This is a validation of the two decompositions against each other and
against the paper's qualitative description (V^ch(0)=0, rapid eigenvalue
decay), not a test of the "does the target's own diagonal exceed 0.5"
question -- Fig. 4 plots V(t), not the target's diagonal, so it cannot
speak to that separate question.

Usage:
    uv run --with numpy --with matplotlib python3 plot_fig4_hybridization.py \\
        --target gbek-atomic-limit-exact-lesser --rank 3
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, eigenvector_decompose
from provenance import write_provenance


def load_lambda(path):
    """Read a gbek_selfconsistency.py::dump_lesser output directly as
    Lambda^<_+ (no extra -1j -- see plot_fig3_errstep.py's module
    docstring for why this differs from plot_errstep.py's loader)."""
    ts, re, im = read_lesser_file(path)
    return ts, re + 1j * im


def hermitian_complete(lam):
    full = np.zeros_like(lam)
    N = lam.shape[0] - 1
    for n in range(N + 1):
        for j in range(N + 1):
            full[n, j] = lam[n, j] if j <= n else np.conj(lam[j, n])
    return full


def eigenvalues_full(full, L):
    w = np.linalg.eigvalsh(full)
    return np.sort(w)[::-1][:max(L + 2, 8)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="gbek-atomic-limit-exact-lesser",
                     help="gbek_selfconsistency.py dump_lesser output (Lambda_+, "
                          "NOT a cincuenta C++ weiss-delta-lesser dump)")
    ap.add_argument("--rank", type=int, default=3)
    ap.add_argument("--out", default="fig4_hybridization.png")
    args = ap.parse_args()

    L = args.rank
    ts, lam = load_lambda(args.target)
    full = hermitian_complete(lam)

    V_ch = cholesky_causal(lam, L)
    V_ev = eigenvector_decompose(full, L)

    print(f"V^ch_0p(t=0) for p=1..{L}: {V_ch[0, :]}  (paper: must be exactly 0)")
    print(f"max|V^ch_0p(0)| = {np.max(np.abs(V_ch[0, :])):.2e}")

    a_p = eigenvalues_full(full, L)
    print("Leading eigenvalues a_p of -i*Lambda^<_+:")
    for p, val in enumerate(a_p):
        flag = "  <-- kept (p <= L)" if p < L else "  (discarded)"
        print(f"  a_{p+1} = {val: .6f}{flag}")

    fig, axes = plt.subplots(2, 1, figsize=(7, 7.5), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, L))

    ax = axes[0]
    for p in range(L):
        ax.plot(ts, V_ch[:, p].real, "-", color=colors[p], lw=1.8, label=f"Re V^ch_0,{p+1}(t)")
        ax.plot(ts, V_ch[:, p].imag, "--", color=colors[p], lw=1.2, alpha=0.7)
    ax.set_title(f"Causal Cholesky, rank L={L}  (cf. GBEK Fig. 4 top panel)", fontsize=10)
    ax.set_ylabel(r"$V^{ch}_{0p}(t)$")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.2)
    ax.axhline(0, color="k", lw=0.5)

    ax = axes[1]
    for p in range(L):
        ax.plot(ts, V_ev[:, p].real, "-", color=colors[p], lw=1.8, label=f"Re V^ev_0,{p+1}(t)")
        ax.plot(ts, V_ev[:, p].imag, "--", color=colors[p], lw=1.2, alpha=0.7)
    ax.set_title(f"Eigenvector approximation, rank L={L}  (cf. GBEK Fig. 4 bottom panel)", fontsize=10)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$V^{ev}_{0p}(t)$")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.2)
    ax.axhline(0, color="k", lw=0.5)

    # Eigenvalue-decay inset (paper places this in the top panel)
    inset = axes[0].inset_axes([0.62, 0.08, 0.34, 0.34])
    idx = np.arange(1, len(a_p) + 1)
    inset.semilogy(idx, np.abs(a_p), "o-", color="k", ms=3, lw=1)
    inset.axvline(L + 0.5, color="r", lw=0.8, ls=":")
    inset.set_xlabel("p", fontsize=7)
    inset.set_ylabel(r"$a_p$", fontsize=7)
    inset.tick_params(labelsize=6)
    inset.set_title("eigenvalue decay", fontsize=7)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nWrote {args.out}")

    prov = write_provenance(args.out, extra_files=[args.target, "plot_fig4_hybridization.py"],
                             notes=f"target={args.target} rank={L}")
    print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
