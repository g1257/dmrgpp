#!/usr/bin/env python3
"""
Row-by-row comparison of the Cholesky factor V_(n,p) between cincuenta's
C++ online computation (dumped via NeqBathDecomposition::dumpV, one row per
time step, finalized during that step's own NeqDmftIter predictor-corrector
loop and never revisited) and an offline batch replay of the identical,
independently-validated gbek_cholesky.py::cholesky_causal() algorithm fed
the same final converged target.

Purpose: pin down exactly which row n the two implementations first
diverge at, to distinguish "small numerical drift that compounds" from "a
hard branch mismatch" -- see project_gbek_cosine_bug.md's "Update
2026-07-07 overnight" section for why this comparison is needed (the two
implementations agree on the near-atomic Fig3L3 target but disagree here).

Usage:
    uv run --with numpy python3 compare_V_rows.py
"""
import numpy as np

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal

TARGET_FILE = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-weiss-delta-lesser"
V_CPP_FILE = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-cholesky-V"
L = 3


def read_V(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 4:
                continue
            rows.append([float(c) for c in cols])
    rows = np.array(rows)
    n_max = int(rows[:, 0].max())
    p_max = int(rows[:, 1].max())
    V = np.zeros((n_max + 1, p_max + 1), dtype=complex)
    for n, p, re, im in rows:
        V[int(n), int(p)] = re + 1j * im
    return V


def main():
    ts, re, im = read_lesser_file(TARGET_FILE)
    lam_target = -1j * (re + 1j * im)

    V_py = cholesky_causal(lam_target, L=L)
    V_cpp = read_V(V_CPP_FILE)

    print(f"V_py shape={V_py.shape}  V_cpp shape={V_cpp.shape}")
    n_max = min(V_py.shape[0], V_cpp.shape[0])

    print(f"{'n':>4} {'|V_py|':>28} {'|V_cpp|':>28} {'max|diff|':>12}")
    first_divergence = None
    for n in range(n_max):
        row_py = V_py[n]
        row_cpp = V_cpp[n]
        diff = np.max(np.abs(row_py - row_cpp))
        if first_divergence is None and diff > 1e-6:
            first_divergence = n
        mag_py = " ".join(f"{abs(v):.5f}" for v in row_py)
        mag_cpp = " ".join(f"{abs(v):.5f}" for v in row_cpp)
        marker = "  <-- diverges" if (first_divergence == n) else ""
        print(f"{n:4d} {mag_py:>28} {mag_cpp:>28} {diff:12.6f}{marker}")

    print()
    if first_divergence is None:
        print("No divergence found above 1e-6 threshold.")
    else:
        n0 = first_divergence
        print(f"First divergence at row n={n0} (t={n0*0.04:.2f})")
        print(f"  V_py[{n0}]  = {V_py[n0]}")
        print(f"  V_cpp[{n0}] = {V_cpp[n0]}")
        print(f"  target diag lam[{n0},{n0}] = {lam_target[n0, n0].real:.6f}")
        L_int = L
        if n0 <= L_int:
            print(f"  n0 <= L={L_int}: divergence is in the SEEDING phase (Eq. 56-58)")
        else:
            print(f"  n0 > L={L_int}: divergence is in the OPTIMAL UPDATE phase (Eq. 62-63)")
            s = n0 - 1
            Q = V_py[1:n0, :]  # Python's Q at this row (uses ALL prior py rows)
            Q_cpp = V_cpp[1:n0, :]  # what Q would be if built from cpp's own prior rows
            print(f"  max|Q_py - Q_cpp| over prior rows 1..{s}: {np.max(np.abs(Q - Q_cpp)):.3e}")


if __name__ == "__main__":
    main()
