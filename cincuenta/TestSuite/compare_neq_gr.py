#!/usr/bin/env python3
"""Compare G^<(t,0) and G^R(t,0) between a tDMRG run and an exact-diag run.

Usage:
    python3 compare_neq_gr.py \\
        --tdmrg-retarded  green-retarded-tdmrg \\
        --tdmrg-lesser    green-lesser-tdmrg \\
        --ed-retarded     green-retarded-ed \\
        --ed-lesser       green-lesser-ed

The tDMRG now computes the full G^R(t,0) = G^>(t,0) - G^<(t,0) directly and
stores it in green-retarded.  G^<(t,0) is stored in green-lesser.
The exact-diag stores G^R(t,0) and G^<(t,0) separately.

File format (KadanoffBaym::dump output):
    t  t'  Re(G)  Im(G)
one row per (t,t') pair, all pairs written.
"""

import argparse
import math
import sys


def read_t0_slice(path):
    """Return dict {t: (re, im)} for all rows where t'=0."""
    data = {}
    with open(path) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 4:
                continue
            try:
                t, tp = float(cols[0]), float(cols[1])
                re, im = float(cols[2]), float(cols[3])
            except ValueError:
                continue
            if abs(tp) < 1e-9:
                data[round(t, 8)] = (re, im)
    return data


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def absdiff(a, b):
    return (abs(a[0] - b[0]), abs(a[1] - b[1]))


def print_table(header, rows, col_headers):
    print(header)
    widths = [max(len(h), 6) for h in col_headers]
    fmt_h = "  ".join(f"{{:>{w}}}" for w in widths)
    print(fmt_h.format(*col_headers))
    print("-" * (sum(widths) + 2 * (len(widths) - 1)))
    for row in rows:
        t_str = f"{row[0]:>5.1f}"
        vals = "  ".join(f"{v:>{widths[i+1]}.5f}" for i, v in enumerate(row[1:]))
        print(f"{t_str}  {vals}")


def summarise(label, errors_re, errors_im):
    print(f"\n{label}")
    print(f"  Max  |dRe|={max(errors_re):.5f}  |dIm|={max(errors_im):.5f}")
    print(f"  Mean |dRe|={sum(errors_re)/len(errors_re):.5f}"
          f"  |dIm|={sum(errors_im)/len(errors_im):.5f}")
    rms_re = math.sqrt(sum(e**2 for e in errors_re) / len(errors_re))
    rms_im = math.sqrt(sum(e**2 for e in errors_im) / len(errors_im))
    print(f"  RMS  |dRe|={rms_re:.5f}  |dIm|={rms_im:.5f}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tdmrg-retarded", default="green-retarded",
                   metavar="FILE", help="tDMRG green-retarded file (stores full G^R)")
    p.add_argument("--tdmrg-lesser",   default="green-lesser",
                   metavar="FILE", help="tDMRG green-lesser file (stores G^<)")
    p.add_argument("--ed-retarded",    default="green-retarded-ed",
                   metavar="FILE", help="Exact-diag green-retarded file (stores G^R)")
    p.add_argument("--ed-lesser",      default="green-lesser-ed",
                   metavar="FILE", help="Exact-diag green-lesser file (stores G^<)")
    p.add_argument("--tolerance-re",   default=None, type=float, metavar="TOL",
                   help="If set, exit 1 when G^R mean |dRe| exceeds TOL")
    p.add_argument("--tolerance-im",   default=None, type=float, metavar="TOL",
                   help="If set, exit 1 when G^R mean |dIm| exceeds TOL")
    args = p.parse_args()

    td_gr  = read_t0_slice(args.tdmrg_retarded)   # tDMRG: G^R(t,0)
    td_gl  = read_t0_slice(args.tdmrg_lesser)      # tDMRG: G^<(t,0)
    ed_gr  = read_t0_slice(args.ed_retarded)       # ED:    G^R(t,0)
    ed_gl  = read_t0_slice(args.ed_lesser)         # ED:    G^<(t,0)

    times = sorted(set(td_gr) & set(ed_gr))
    if not times:
        print("ERROR: no common time points found -- check file paths", file=sys.stderr)
        sys.exit(1)

    # ---- G^R(t,0) comparison ------------------------------------------------
    cols = ["t", "Re tDMRG", "Im tDMRG", "Re ED G^R", "Im ED G^R", "|dRe|", "|dIm|"]
    rows_gr, dre_gr, dim_gr = [], [], []
    for t in times:
        r_td, i_td = td_gr[t]
        r_ed, i_ed = ed_gr[t]
        dr, di = abs(r_td - r_ed), abs(i_td - i_ed)
        dre_gr.append(dr); dim_gr.append(di)
        rows_gr.append((t, r_td, i_td, r_ed, i_ed, dr, di))

    print_table("\nG^R(t,0): tDMRG (full G^R = G^> - G^<) vs exact-diag",
                rows_gr, cols)
    summarise("G^R(t,0) errors", dre_gr, dim_gr)

    # ---- G^<(t,0) comparison ------------------------------------------------
    times_gl = sorted(set(td_gl) & set(ed_gl))
    if times_gl:
        cols2 = ["t", "Re tDMRG", "Im tDMRG", "Re ED G^<", "Im ED G^<", "|dRe|", "|dIm|"]
        rows_gl, dre_gl, dim_gl = [], [], []
        for t in times_gl:
            r_td, i_td = td_gl[t]
            r_ed, i_ed = ed_gl[t]
            dr, di = abs(r_td - r_ed), abs(i_td - i_ed)
            dre_gl.append(dr); dim_gl.append(di)
            rows_gl.append((t, r_td, i_td, r_ed, i_ed, dr, di))

        print_table("\nG^<(t,0): tDMRG vs exact-diag",
                    rows_gl, cols2)
        summarise("G^<(t,0) errors", dre_gl, dim_gl)

    # ---- G^>(t,0) from tDMRG (reconstructed) --------------------------------
    # G^>(t,0) = G^R(t,0) + G^<(t,0)
    times_ggt = sorted(set(td_gr) & set(td_gl) & set(ed_gr) & set(ed_gl))
    if times_ggt:
        ed_ggt = {t: add(ed_gr[t], ed_gl[t]) for t in times_ggt}
        td_ggt = {t: add(td_gr[t], td_gl[t]) for t in times_ggt}
        cols3 = ["t", "Re tDMRG", "Im tDMRG", "Re ED G^>", "Im ED G^>", "|dRe|", "|dIm|"]
        rows_ggt, dre_ggt, dim_ggt = [], [], []
        for t in times_ggt:
            r_td, i_td = td_ggt[t]
            r_ed, i_ed = ed_ggt[t]
            dr, di = abs(r_td - r_ed), abs(i_td - i_ed)
            dre_ggt.append(dr); dim_ggt.append(di)
            rows_ggt.append((t, r_td, i_td, r_ed, i_ed, dr, di))

        print_table("\nG^>(t,0): tDMRG (reconstructed = G^R + G^<) vs exact-diag",
                    rows_ggt, cols3)
        summarise("G^>(t,0) errors", dre_ggt, dim_ggt)


    if args.tolerance_re is not None or args.tolerance_im is not None:
        mean_re = sum(dre_gr) / len(dre_gr)
        mean_im = sum(dim_gr) / len(dim_gr)
        failed = False
        if args.tolerance_re is not None and mean_re > args.tolerance_re:
            print(f"\nFAIL: G^R mean |dRe|={mean_re:.5f} exceeds tolerance {args.tolerance_re}",
                  file=sys.stderr)
            failed = True
        if args.tolerance_im is not None and mean_im > args.tolerance_im:
            print(f"\nFAIL: G^R mean |dIm|={mean_im:.5f} exceeds tolerance {args.tolerance_im}",
                  file=sys.stderr)
            failed = True
        if failed:
            sys.exit(1)


if __name__ == "__main__":
    main()
