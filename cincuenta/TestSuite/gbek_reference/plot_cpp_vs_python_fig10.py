#!/usr/bin/env python3
"""
Compare the production C++ GBEK solver's docc(t) dumps against the Python
reference's Fig. 10 reproduction (run_fig10_scan.py), for the true atomic
limit (NeqAtomicLimit=1, nBath=0), U=0,2,4 (the subset the C++ side has
data for -- Python's own Fig. 10 also covers U=1,6,8, not reproduced in
C++ here since those combos were never run through cincuenta), across all
three Lbath=4,6,8 (L=2,3,4).

One panel per Lbath, Python (solid) and C++ (dashed) overlaid directly per
U so any mismatch is immediately visible as separated line pairs, rather
than side-by-side panels as in plot_cpp_vs_python_fig9.py.

Data sources: same as plot_cpp_vs_python_fig9.py (cpp_docc_energy/, Python
npz files) -- reuses the docc(t) column already present in the C++
docc-energy dumps produced for the Fig. 9 comparison, no new C++ runs
needed.

Usage:
    uv run --with numpy --with matplotlib python3 plot_cpp_vs_python_fig10.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from provenance import write_provenance

U_VALUES = (0.0, 2.0, 4.0)
L_VALUES = (2, 3, 4)
U_COLORS = {0.0: "tab:red", 2.0: "tab:blue", 4.0: "tab:green"}


def load_cpp(U, L):
    suffix = "" if U == 2.0 else f"-U{U:g}"
    fname = f"cpp_docc_energy/atomic-limit-gbek-L{L}{suffix}-docc-energy"
    raw = np.loadtxt(fname)
    return raw[:, 0], raw[:, 1]  # ts, docc_t


def load_py(U, L):
    d = np.load(f"fig10_docc_U{U:g}_L{L}.npz")
    return d["ts"], d["d_t"]


def main():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)

    for ax, L in zip(axes, L_VALUES):
        for U in U_VALUES:
            ts_py, d_py = load_py(U, L)
            ts_cpp, d_cpp = load_cpp(U, L)
            color = U_COLORS[U]
            ax.plot(ts_py, d_py, "-", color=color, lw=2.0, label=f"U={U:g} Python")
            ax.plot(ts_cpp, d_cpp, "--", color=color, lw=1.4, label=f"U={U:g} C++")
        ax.set_title(f"Lbath={2*L}", fontsize=10)
        ax.set_xlabel("t")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=6, ncol=2, loc="upper left")

    axes[0].set_ylabel("d(t)")
    fig.suptitle("cf. GBEK Fig. 10: Python (solid) vs. C++ (dashed), U=0,2,4")
    fig.tight_layout()
    out = "cpp_vs_python_fig10.png"
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")

    extra_files = (["plot_cpp_vs_python_fig10.py"]
                   + [f"fig10_docc_U{U:g}_L{L}.npz" for U in U_VALUES for L in L_VALUES]
                   + [f"cpp_docc_energy/atomic-limit-gbek-L{L}{'' if U == 2.0 else f'-U{U:g}'}-docc-energy"
                      for U in U_VALUES for L in L_VALUES])
    prov = write_provenance(out, extra_files=extra_files,
                             notes="figure=10 C++-vs-Python, U=0,2,4 subset")
    print(f"Wrote {prov}")

    print("\n--- C++ vs Python: d(t) agreement ---")
    for U in U_VALUES:
        for L in L_VALUES:
            ts_py, d_py = load_py(U, L)
            ts_cpp, d_cpp = load_cpp(U, L)
            d_py_interp = np.interp(ts_cpp, ts_py, d_py)
            diff = np.max(np.abs(d_py_interp - d_cpp))
            print(f"  U={U:g} Lbath={2*L}: max|d_cpp - d_python| = {diff:.4e}")


if __name__ == "__main__":
    main()
