#!/usr/bin/env python3
"""
Compare the exact atomic-limit reference -i*Lambda^+_<(t,t') (produced by
gbek_selfconsistency.py) against cincuenta's own rank-L Cholesky
approximation (dumped by ImpuritySolverNeqGBEK::dumpPlusBath, e.g.
"gebk-fig3-L3-plus-bath-lesser").

Usage:
    python3 compare_reference.py gbek-atomic-limit-exact-lesser \\
        gebk-fig3-L3-plus-bath-lesser --tmax 4.0

Also reused for run-vs-run comparisons where NEITHER file is the independent
Python target (e.g. two different cincuenta runs against each other) --
ALWAYS pass --ref-label/--approx-label explicitly in that case. The
defaults describe this script's original single use case (independent
Python exact target vs. cincuenta reconstruction) and are actively
misleading if left in place for any other comparison:

    python3 compare_reference.py atomic-limit-gbek-L2-plus-bath-lesser \\
        gebk-fig3-L2-plus-bath-lesser --tmax 4.0 \\
        --ref-label "True atomic limit (Lambda_-=0), rank-2 Cholesky" \\
        --approx-label "Near-atomic small bath (W_i=0.1, 5 sites), rank-2 Cholesky"
"""
import argparse
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from provenance import write_provenance
from gbek_colormap import GBEK_CMAP, GBEK_VMIN, GBEK_VMAX


def read_lesser_file(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 4:
                continue
            rows.append([float(c) for c in cols[:4]])
    rows = np.array(rows)
    ts = np.unique(rows[:, 0])
    N = len(ts)
    re = np.zeros((N, N))
    im = np.zeros((N, N))
    idx = {t: i for i, t in enumerate(ts)}
    for t, tp, r, i in rows:
        n, j = idx[t], idx[tp]
        re[n, j] = r
        im[n, j] = i
    return ts, re, im


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference", help="exact reference file (gbek_selfconsistency.py output)")
    ap.add_argument("approx", help="cincuenta's rank-L Cholesky output (dumpPlusBath)")
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--rank", type=int, default=None,
                    help="Cholesky rank L used for --approx, shown in its plot title "
                         "(ignored if --approx-label is given)")
    ap.add_argument("--ref-label", default="Exact reference",
                    help="Label for the 'reference' file, used in its panel title, "
                         "legend entries, and the residual title. This script is "
                         "reused for comparisons where 'reference' is NOT the "
                         "independent, undecomposed Python target (e.g. two cincuenta "
                         "runs against each other) -- always set this explicitly in "
                         "that case so the plot doesn't silently claim to show "
                         "something it doesn't. Default matches this script's "
                         "original, single intended use case (independent Python "
                         "exact target vs. cincuenta reconstruction).")
    ap.add_argument("--approx-label", default=None,
                    help="Label for the 'approx' file, used the same way as "
                         "--ref-label. Defaults to 'cincuenta rank-L Cholesky approx' "
                         "(rank filled in from --rank if given) to match this "
                         "script's original default wording.")
    ap.add_argument("--note", default=None,
                    help="Optional one-line annotation printed above the figure, e.g. "
                         "to quantify a parameter mentioned in --ref-label/"
                         "--approx-label (\"small\" is meaningless without a scale to "
                         "compare it to) -- e.g. --note \"Residual first bath W_i=0.1 "
                         "is 2.5%% of the final bandwidth W_f=4\".")
    ap.add_argument("--out", default="gbek_reference_comparison.png")
    args = ap.parse_args()

    rank_label = f"rank-{args.rank}" if args.rank is not None else "rank-L"
    approx_label = args.approx_label or f"cincuenta {rank_label} Cholesky approx"
    ref_label = args.ref_label

    ts_ref, re_ref, im_ref = read_lesser_file(args.reference)
    ts_app, re_app, im_app = read_lesser_file(args.approx)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5))

    # Use the paper's own colormap and value range (Fig. 3's colorbar spans
    # -0.2 to 0.5, asymmetric -- not a plain symmetric vmin=-vmax) so these
    # two panels are directly, visually comparable to the paper's own
    # figure, not just numerically.
    # Long, explicit labels (needed to avoid mislabeling when this script is
    # reused for non-default comparisons -- see --ref-label/--approx-label
    # help text) don't fit on one line at a readable size, so wrap them
    # rather than letting them overflow into neighboring panels.
    def wrap_title(text, width=34):
        return "\n".join(textwrap.wrap(text, width=width))

    im0 = axes[0, 0].imshow(re_ref, origin="lower", extent=[ts_ref[0], ts_ref[-1]] * 2,
                             vmin=GBEK_VMIN, vmax=GBEK_VMAX, cmap=GBEK_CMAP)
    axes[0, 0].set_title(wrap_title(f"{ref_label}: Re[-i Lambda^<_+(t,t')]"), fontsize=9)
    plt.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].imshow(re_app, origin="lower", extent=[ts_app[0], ts_app[-1]] * 2,
                             vmin=GBEK_VMIN, vmax=GBEK_VMAX, cmap=GBEK_CMAP)
    axes[0, 1].set_title(wrap_title(f"{approx_label}: Re[-i Lambda^<_+(t,t')]"), fontsize=9)
    plt.colorbar(im1, ax=axes[0, 1])

    # Interpolate approx onto the reference grid for a residual map, if grids differ
    if not np.allclose(ts_ref, ts_app):
        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator((ts_app, ts_app), re_app,
                                          bounds_error=False, fill_value=np.nan)
        grid_n, grid_j = np.meshgrid(ts_ref, ts_ref, indexing="ij")
        re_app_interp = interp((grid_n, grid_j))
    else:
        re_app_interp = re_app

    resid = re_ref - re_app_interp
    # Symmetric about 0 (so white == no error, unlike GBEK_CMAP's own
    # asymmetric white-point) and scaled to the same order of magnitude as
    # the exact/approx plots to its left (GBEK_VMAX), so the residual's
    # size relative to the actual values is directly visible instead of
    # auto-scaling to whatever (typically much smaller) range the residual
    # itself spans.
    resid_scale = max(abs(GBEK_VMIN), GBEK_VMAX)
    im2 = axes[0, 2].imshow(resid, origin="lower", extent=[ts_ref[0], ts_ref[-1]] * 2,
                             vmin=-resid_scale, vmax=resid_scale, cmap="RdBu_r")
    axes[0, 2].set_title(wrap_title(f"Residual ({ref_label} - {approx_label})", width=40),
                          fontsize=8)
    plt.colorbar(im2, ax=axes[0, 2])

    # Diagonal comparison
    diag_ref = np.diag(re_ref)
    diag_app = np.diag(re_app)
    axes[1, 0].plot(ts_ref, diag_ref, label=ref_label, lw=2)
    axes[1, 0].plot(ts_app, diag_app, "--", label=approx_label, lw=2)
    axes[1, 0].axhline(0.5, color="gray", ls=":", lw=1, label="exact p-h value (0.5)")
    axes[1, 0].set_xlabel("t")
    axes[1, 0].set_ylabel("Re[-i Lambda(t,t)]")
    axes[1, 0].set_title("Diagonal")
    axes[1, 0].legend(fontsize=8)

    # A couple of off-diagonal (t, t'=fixed) slices
    for tp_fixed, ax_idx in [(0.4, 1), (1.5, 2)]:
        ax = axes[1, ax_idx]
        j_ref = np.argmin(np.abs(ts_ref - tp_fixed))
        j_app = np.argmin(np.abs(ts_app - tp_fixed))
        ax.plot(ts_ref, re_ref[:, j_ref], label=ref_label, lw=2)
        ax.plot(ts_app, re_app[:, j_app], "--", label=approx_label, lw=2)
        ax.set_xlabel("t")
        ax.set_title(f"Slice at t'={tp_fixed}")
        ax.legend(fontsize=8)

    for ax in axes.flat:
        if args.tmax:
            ax.set_xlim(0, args.tmax)
            if hasattr(ax, "set_ylim") and ax in axes[0, :]:
                ax.set_ylim(0, args.tmax)

    fig.tight_layout()
    if args.note:
        fig.subplots_adjust(top=0.90)
        fig.text(0.5, 0.965, args.note, ha="center", fontsize=10,
                  style="italic", color="#333333")
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")
    prov = write_provenance(args.out, extra_files=[args.reference, args.approx],
                             notes=f"reference={args.reference} ({ref_label}) "
                                   f"approx={args.approx} ({approx_label})")
    print(f"Wrote {prov}")

    print("\nSummary:")
    print(f"  max|diag_ref - 0.5| (post t=0.4): "
          f"{np.max(np.abs(diag_ref[ts_ref > 0.4] - 0.5)):.4f}")
    print(f"  max|diag_app - 0.5| (post t=0.4): "
          f"{np.max(np.abs(diag_app[ts_app > 0.4] - 0.5)):.4f}")
    print(f"  max|residual|: {np.nanmax(np.abs(resid)):.4f}")


if __name__ == "__main__":
    main()
