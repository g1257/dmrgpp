#!/usr/bin/env python3
"""
Plot err^step(t) (GBEK PRB 88, 235106 (2013), Eq. 67 -- verified by actually
compiling the paper's LaTeX source and reading the real REVTeX-assigned
number from the .aux file's \\newlabel entries; see
project_gbek_paper_source_reliability memory) for a given target Lambda,
in the same style as the paper's own Fig. 3 bottom-left panel, so a plot
made here can be visually compared directly against that panel.

    err^step(A, tau) = sqrt( sum_{n=1}^{N} (2 - delta_{nN})
                              |Lambda^<_+(n,N) - A(n,N)|^2 ),   N = tau/dt

i.e. the residual on just the newly-added row/column at time step N, not
the whole matrix up to N -- this is what makes it "stepwise": it shows how
well the approximation would have looked if you'd stopped propagating at
that exact step, which is the causal algorithm's natural unit of progress.

Always plots BOTH the independently-computed Python reference curves
(gbek_cholesky.py::cholesky_causal, and a plain eigenvector low-rank
decomposition, both computed here from the target directly) AND,
optionally, cincuenta's own dumped C++ reconstruction (--cpp-recon) so the
two can be checked against each other -- this is the intended regression
check: if cincuenta's online reconstruction's err^step curve drifts away
from the Python reference's cholesky curve at the same rank, something in
NeqBathDecomposition.h changed behavior (for better or worse).

Usage:
    uv run --with numpy --with matplotlib python3 plot_errstep.py \\
        --target /path/to/..-weiss-delta-lesser \\
        [--cpp-recon /path/to/..-plus-bath-lesser] \\
        [--ranks 2,3] [--out errstep.png] [--title "..."]

Defaults reproduce this session's atomic-limit L=3 investigation.
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, reconstruct, eigenvector_decompose
from provenance import write_provenance

DEFAULT_TARGET = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-weiss-delta-lesser"


def load_lambda(path):
    """Read a *-weiss-delta-lesser dump and convert to Lambda (= -i * the
    dump's raw lesser component).

    Lambda^- == 0 for the true-atomic-limit runs this script targets, so
    this IS Lambda_+ directly (see project_neq_atomic_limit_implementation
    memory). If you ever point this at a target where Lambda^- != 0, this
    conversion is no longer valid as "Lambda_+" and the script would need
    a Lambda^- subtraction first (cf. gbek_cholesky.py's module docstring
    and the near-atomic W_i=0.1 analysis in project_gbek_cosine_bug).
    """
    ts, re, im = read_lesser_file(path)
    return ts, -1j * (re + 1j * im)


def load_recon(path):
    """Read a *-plus-bath-lesser dump: already in Lambda units, no -i needed."""
    _, re, im = read_lesser_file(path)
    return re + 1j * im


def errstep_curve(lam, A):
    """Eq. 67, vectorized over n for each N."""
    N = lam.shape[0] - 1
    out = np.full(N + 1, np.nan)
    for Nidx in range(1, N + 1):
        diff = lam[1:Nidx + 1, Nidx] - A[1:Nidx + 1, Nidx]
        w = np.full(Nidx, 2.0)
        w[-1] = 1.0
        out[Nidx] = np.sqrt(np.sum(w * np.abs(diff) ** 2))
    return out


def global_err(lam, A):
    """err[A] = ||Lambda - A|| / ||Lambda||, the paper's single-number
    global relative Frobenius error (quoted for their own Fig. 3 data as
    err[ch]=0.17, err[ev]=0.09 at rank L=3)."""
    return np.linalg.norm(lam - A) / np.linalg.norm(lam)


def eigenvector_lowrank(lam, L):
    return reconstruct(eigenvector_decompose(lam, L))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default=DEFAULT_TARGET,
                    help="*-weiss-delta-lesser dump (the exact target Lambda). "
                         "Ignored if --run is given.")
    ap.add_argument("--cpp-recon", default=None,
                    help="optional *-plus-bath-lesser dump (cincuenta's own "
                         "online reconstruction) to overlay against the Python "
                         "reference curves at --cpp-recon-rank. Ignored if --run "
                         "is given.")
    ap.add_argument("--cpp-recon-rank", type=int, default=3,
                    help="rank cincuenta was run with, for the --cpp-recon "
                         "curve's legend label only (default 3)")
    ap.add_argument("--ranks", default="2,3",
                    help="comma-separated ranks for the Python ch/ev reference "
                         "curves, all computed against --target. Ignored if "
                         "--run is given.")
    ap.add_argument("--run", action="append", default=None,
                    metavar="TARGET,RECON,RANK",
                    help="one full cincuenta run to overlay: its own "
                         "*-weiss-delta-lesser target, its own "
                         "*-plus-bath-lesser C++ reconstruction, and the rank "
                         "it was run at, comma-separated. Repeat --run for "
                         "each rank/run you want on the same plot (e.g. to "
                         "compare L=3, L=4, L=5 side by side, each against ITS "
                         "OWN self-consistent target -- these differ slightly "
                         "run to run since the online scheme's target depends "
                         "on the reconstruction quality at every prior step). "
                         "Python ch/ev reference curves are computed at each "
                         "run's own rank from its own target. Overrides "
                         "--target/--cpp-recon/--cpp-recon-rank/--ranks.")
    ap.add_argument("--out", default="errstep_paper_style.png")
    ap.add_argument("--title", default="cf. GBEK Fig. 3 bottom-left panel")
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(7, 5.8))
    print(f"{'curve':<24} {'err^step(t=4)':>14} {'min err^step(t>0.5)':>20} {'global err[A]':>14}")
    extra_files = []

    if args.run:
        runs = []
        for spec in args.run:
            target, recon, rank = spec.split(",")
            runs.append((target, recon, int(rank)))
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(runs)))

        for (target, recon, L), color in zip(runs, colors):
            ts, lam = load_lambda(target)
            V = cholesky_causal(lam, L)
            ch = reconstruct(V)
            ev = eigenvector_lowrank(lam, L)
            cpp = load_recon(recon)

            ch_curve = errstep_curve(lam, ch)
            ev_curve = errstep_curve(lam, ev)
            cpp_curve = errstep_curve(lam, cpp)

            ax.semilogy(ts, ch_curve, '-', color=color, lw=2.0,
                        label=rf"rank {L}: $(\Lambda_+^<)^{{ch}}$ (python ref)")
            ax.semilogy(ts, ev_curve, '--', color=color, lw=1.3, alpha=0.7,
                        label=rf"rank {L}: $(\Lambda_+^<)^{{ev}}$ (python ref)")
            ax.semilogy(ts, cpp_curve, ':', color=color, lw=2.4,
                        label=rf"rank {L}: cincuenta C++ online")

            print(f"{'ch L='+str(L):<24} {ch_curve[-1]:>14.4f} "
                  f"{np.nanmin(ch_curve[ts > 0.5]):>20.4f} {global_err(lam, ch):>14.4f}")
            print(f"{'ev L='+str(L):<24} {ev_curve[-1]:>14.4f} "
                  f"{np.nanmin(ev_curve[ts > 0.5]):>20.4f} {global_err(lam, ev):>14.4f}")
            print(f"{'cpp L='+str(L):<24} {cpp_curve[-1]:>14.4f} "
                  f"{np.nanmin(cpp_curve[ts > 0.5]):>20.4f} {global_err(lam, cpp):>14.4f}")
            extra_files += [target, recon]

        notes = f"runs={args.run}"
    else:
        ranks = [int(x) for x in args.ranks.split(",")]
        ts, lam = load_lambda(args.target)
        ch_colors = plt.cm.Reds(np.linspace(0.5, 0.9, len(ranks)))
        ev_colors = plt.cm.Greens(np.linspace(0.5, 0.9, len(ranks)))

        for L, cch, cev in zip(ranks, ch_colors, ev_colors):
            V = cholesky_causal(lam, L)
            ch = reconstruct(V)
            ev = eigenvector_lowrank(lam, L)

            ch_curve = errstep_curve(lam, ch)
            ev_curve = errstep_curve(lam, ev)

            ax.semilogy(ts, ch_curve, '-', color=cch, lw=1.8,
                        label=rf"$(\Lambda_+^<)^{{ch}}$, rank {L} (python ref)")
            ax.semilogy(ts, ev_curve, '--', color=cev, lw=1.6,
                        label=rf"$(\Lambda_+^<)^{{ev}}$, rank {L} (python ref)")

            print(f"{'ch L='+str(L):<24} {ch_curve[-1]:>14.4f} "
                  f"{np.nanmin(ch_curve[ts > 0.5]):>20.4f} {global_err(lam, ch):>14.4f}")
            print(f"{'ev L='+str(L):<24} {ev_curve[-1]:>14.4f} "
                  f"{np.nanmin(ev_curve[ts > 0.5]):>20.4f} {global_err(lam, ev):>14.4f}")

        if args.cpp_recon:
            cpp = load_recon(args.cpp_recon)
            cpp_curve = errstep_curve(lam, cpp)
            ax.semilogy(ts, cpp_curve, ':', color='black', lw=2.2,
                        label=f"cincuenta C++ online, rank {args.cpp_recon_rank}")
            print(f"{'cpp L='+str(args.cpp_recon_rank):<24} {cpp_curve[-1]:>14.4f} "
                  f"{np.nanmin(cpp_curve[ts > 0.5]):>20.4f} {global_err(lam, cpp):>14.4f}")

        extra_files = [args.target] + ([args.cpp_recon] if args.cpp_recon else [])
        notes = f"target={args.target} cpp_recon={args.cpp_recon} ranks={ranks}"

    ax.set_xlim(0, ts[-1])
    ax.set_ylim(1e-6, max(3, ax.get_ylim()[1]))
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\mathrm{err}^{\mathrm{step}}(t)$  (Eq. 67)")
    ax.set_title(args.title, fontsize=10)
    ax.legend(fontsize=7.5, loc='lower right', ncol=1)
    ax.grid(alpha=0.2, which='both')
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nWrote {args.out}")

    prov = write_provenance(args.out, extra_files=extra_files, notes=notes)
    print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
