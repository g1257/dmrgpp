#!/usr/bin/env python3
"""
Compare the production C++ GBEK solver's new docc/energy observables
(ImpuritySolverNeqGBEK::dumpDoccAndEnergy, ported 2026-07-15) against the
Python reference's Fig. 9 reproduction (run_fig9_scan.py /
gbek_selfconsistency.py::compute_energy_observables), for the true atomic
limit (NeqAtomicLimit=1, nBath=0) at U=2,4 and Lbath=4,6,8 (L=2,3,4).

L=4 (Lbath=8) originally hit a separate C++-only bug: lanczosGS()
(used only when nsites_ext_>8, i.e. L>=4) seeded Lanczos with a naive
uniform vector, converging to the WRONG extremal state (reported ground
energy -2807.96 vs. the true product-state energy -4000 for L=4) instead
of the true pre-quench ground state, giving a wrong t=0 state
(docc(0)~0.24, Ekin(0)~0.47 instead of the exact atomic-limit 0/0).
FIXED 2026-07-16: lanczosGS() now seeds Lanczos with the known analytic
atomic-limit product state instead (exact for nBath=0, since the
pre-quench Hamiltonian then has no hopping at all) -- see project memory
project_gbek_cpp_lanczosGS_bug for the full diagnosis and fix. L=2,3 use
full diagonalization (nsites_ext_<=8) and were never affected.

C++ data source: cincuenta/TestSuite/gbek_reference/cpp_docc_energy/
  atomic-limit-gbek-L{L}[-U{U}]-docc-energy, format "t docc Ekin Eint Etot"
  (U=2 files have no -U{U} suffix, matching the original single-U inputs).
Python data source: fig9_energy_U{U}_L{L}.npz (run_fig9_scan.py's output).

Usage:
    uv run --with numpy --with matplotlib python3 plot_cpp_vs_python_fig9.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from provenance import write_provenance

U_VALUES = (2.0, 4.0)     # the paper's actual Fig. 9 curves
L_VALUES_CPP = (2, 3, 4)
L_VALUES_PY = (2, 3, 4)
QUANTITIES = [("Ekin_t", "Ekin", "green"), ("Eint_t", "Eint", "blue"), ("Etot_t", "Etot", "red")]
U_STYLES = {2.0: "-", 4.0: "--"}
L_WIDTHS = {2: 1.2, 3: 2.0, 4: 2.8}


def load_cpp(U, L):
    suffix = "" if U == 2.0 else f"-U{U:g}"
    fname = f"cpp_docc_energy/atomic-limit-gbek-L{L}{suffix}-docc-energy"
    raw = np.loadtxt(fname)
    return {"ts": raw[:, 0], "docc_t": raw[:, 1], "Ekin_t": raw[:, 2],
            "Eint_t": raw[:, 3], "Etot_t": raw[:, 4]}


def load_py(U, L):
    return np.load(f"fig9_energy_U{U:g}_L{L}.npz")


def main():
    py_data = {(U, L): load_py(U, L) for U in U_VALUES for L in L_VALUES_PY}
    cpp_data = {(U, L): load_cpp(U, L) for U in U_VALUES for L in L_VALUES_CPP}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    ax_py, ax_cpp = axes

    for ax, dataset, L_list, title in [
        (ax_py, py_data, L_VALUES_PY, "Python reference (gbek_selfconsistency.py)"),
        (ax_cpp, cpp_data, L_VALUES_CPP, "C++ production solver (ImpuritySolverNeqGBEK)"),
    ]:
        tq = None
        for U in U_VALUES:
            for L in L_list:
                d = dataset[(U, L)]
                ts = d["ts"]
                if "tq" in getattr(d, "files", []):
                    tq = float(d["tq"])
                for key, label, color in QUANTITIES:
                    ax.plot(ts, d[key], color=color, linestyle=U_STYLES[U],
                            linewidth=L_WIDTHS[L],
                            label=f"{label} U={U:g} Lbath={2*L}")
        if tq is not None:
            ax.axvline(tq, color="gray", linestyle=":", linewidth=1.0)
        ax.set_xlabel("t")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=6, ncol=2, loc="lower right")

    ax_py.set_ylabel("energy (units of v0)")
    fig.suptitle("cf. GBEK Fig. 9: Python reference vs. C++ production solver")
    fig.tight_layout()
    out = "cpp_vs_python_fig9.png"
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")

    extra_files = (["plot_cpp_vs_python_fig9.py"]
                   + [f"fig9_energy_U{U:g}_L{L}.npz" for U in U_VALUES for L in L_VALUES_PY]
                   + [f"cpp_docc_energy/atomic-limit-gbek-L{L}{'' if U == 2.0 else f'-U{U:g}'}-docc-energy"
                      for U in U_VALUES for L in L_VALUES_CPP])
    prov = write_provenance(out, extra_files=extra_files,
                             notes="figure=9 C++-vs-Python, all three Lbath")
    print(f"Wrote {prov}")

    print("\n--- C++ vs Python: Etot(t) agreement ---")
    for U in U_VALUES:
        for L in L_VALUES_CPP:
            py = py_data[(U, L)]
            cpp = cpp_data[(U, L)]
            py_etot = np.interp(cpp["ts"], py["ts"], py["Etot_t"])
            diff = np.max(np.abs(py_etot - cpp["Etot_t"]))
            print(f"  U={U:g} Lbath={2*L}: max|Etot_cpp - Etot_python| = {diff:.4e}")

    print("\n--- C++ t=0 sanity check: should be exactly docc=0, Ekin=0, "
          "Etot=-U/4 for all Lbath ---")
    for U in U_VALUES:
        for L in L_VALUES_CPP:
            raw = np.loadtxt(f"cpp_docc_energy/atomic-limit-gbek-L{L}"
                              f"{'' if U == 2.0 else f'-U{U:g}'}-docc-energy")
            print(f"  U={U:g} Lbath={2*L}: t=0 docc={raw[0,1]:.6f} Ekin={raw[0,2]:.6f} "
                  f"Etot={raw[0,4]:.6f} (expected 0, 0, {-U/4:.6f})")


if __name__ == "__main__":
    main()
