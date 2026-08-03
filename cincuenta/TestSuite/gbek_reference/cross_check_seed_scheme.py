"""
Decisive, from-scratch cross-check of the "seed-once-forward" vs
"reseed-at-each-time" two-time Green's function constructions. This is the
script that found and confirmed a real bug in cincuenta's C++ solver
(ImpuritySolverNeqGBEK.h), fixed the same day (2026-07-13) -- see
propagateOneStep's and gLesserRowGBEKSector's doc comments in that file for
the fixed implementation. This script is kept as the PERMANENT source of the
brute-force ground-truth values hardcoded into the C++ regression test
cincuenta/src/tests/test_ImpuritySolverNeqGBEK.cpp's "G^< two-time
construction matches an independent Python brute-force ground truth"
TEST_CASE (the `expected` array there is G_bruteforce's off-diagonal
entries, printed by this script). If the seeding or propagation scheme is
ever intentionally changed, rerun this script and update that TEST_CASE's
expected values from its printed G_bruteforce table.

Background (the bug, now fixed)
--------------------------------
C++ used to do (WRONG, fixed 2026-07-13):
    PsiHist[0] = c_{0,up} |GS_pre-quench>        (seed ONCE, at t=0)
    PsiHist[n] = U_{N-1}(t_n, 0) PsiHist[0]      (propagate FORWARD, N-1 sector
                                                   Hamiltonian only, for all n)
    G_cpp(n,j) = i <PsiHist[j] | PsiHist[n]>
This is "G_cpp-style" below -- kept only as a negative control, to document
what the bug looked like and confirm this script still detects it.

Python (gbek_dynamics.py::compute_g_lesser) does the CORRECT construction,
matching the fixed C++ code:
    phi(t_n)   = U_N(t_n, 0) |GS_pre-quench>     (N-sector reference trajectory)
    b(t_n)     = c_{0,up} phi(t_n)               (reseed AT EACH time n)
    G_py(n,k)  = i <b(t_k) | U_{N-1}(t_k, t_n) | b(t_n)>   (k <= n)

These two formulas are identical only if c_{0,up} commutes with the full
propagator (equivalently, if phi(t_n) never needs to be re-consulted after
t=0).  That requires [H, c] = 0, which is false whenever the Hubbard U term
is active (U n_up n_dn does not commute with c_up).  So the two constructions
provably disagree once U != 0, growing with |n-j| (more time for the error
to accumulate) while leaving G(n,n) -- which only involves phi(t_n) via its
own norm -- comparatively unaffected at early times. This is exactly the
signature the original bug showed in production runs.

This script builds BOTH constructions from scratch, plus a fully independent
"brute-force" dense-matrix-exponential ground truth (does not call
propagate()/compute_g_lesser() at all, so it can't inherit a shared bug),
and compares all three on the same small L=1-style star-geometry system
prior cross-check scripts already validated (nsites=3: impurity=site0, plus
sites 1,2 standing in for the L=1 second bath's "empty"/"occupied" legs --
the SAME toy setup as cross_check_gbek_hamiltonian.py /
cross_check_gbek_propagation.py).

psi0 (the N-sector pre-quench seed) is NOT picked arbitrarily -- it is
derived from and verified against the already-established N-1-sector seed
used in cross_check_gbek_propagation.py; see the derivation and the
assertion immediately after psi0 is constructed below.
"""
import numpy as np
from scipy.linalg import expm

from gbek_ed import build_basis, c_operator
from gbek_dynamics import build_hamiltonian, propagate, compute_g_lesser

# ---- physical setup (matches cross_check_gbek_hamiltonian.py conventions) ----
nsites = 3
U = 2.0
tq = 0.25
tstar_f = 1.0
dt = 0.05
nsteps = 6  # -> t up to 0.30, well past tq=0.25 so hop(t) has saturated to 1 for some steps


def v_ramp(t, tq):
    if t <= 0:
        return 0.0
    if t >= tq:
        return 1.0
    return 0.5 * (1 - np.cos(np.pi * t / tq))


def hop(t):
    return tstar_f * v_ramp(t, tq)


# N sector (2,1); N-1 sector (1,1): c_{0,up} removes the impurity's up
# electron.
basis_N, index_N = build_basis(nsites, 2, 1)
basis_Nm1, index_Nm1 = build_basis(nsites, 1, 1)

c_op = c_operator(nsites, basis_N, index_N, basis_Nm1, index_Nm1, site=0, spin=0)

# psi0 = the ACTUAL pre-quench N-sector ground state cincuenta seeds from
# (ImpuritySolverNeqGBEK.h::buildSector, eigvecsN's column 0), NOT an
# arbitrary hop=0 diagonalization pick.
#
# Why not just diagonalize build_hamiltonian(..., hop=[0,0])?  That toy
# Hamiltonian only has the Hubbard-U diagonal term (U*(n0up-1/2)(n0dn-1/2))
# -- it omits the pre-quench extended-Fock-space site potentials
# (potPre[1]=+bigEps forcing the "empty" second-bath site empty,
# potPre[2]=-bigEps forcing the "occupied" second-bath site occupied) that
# cincuenta actually uses to fix the initial second-bath occupation. Without
# those potentials sites 1 and 2 are degenerate at hop=0, so
# np.linalg.eigh's column-0 pick for the "ground state" is an ARBITRARY one
# of several degenerate basis states -- verified directly: it returns
# (up=3,dn=2) (site0 up, site1 up+dn occupied), not the physical seed.
#
# The physical seed is instead derived from the ALREADY-VALIDATED N-1-sector
# seed used in cross_check_gbek_propagation.py and the "Time propagation..."
# C++ unit test: that seed is the pure basis state (up=4,dn=4) in the (1,1)
# sector (site 2 doubly occupied, impurity and site 1 empty -- matching the
# production pre-quench potentials: site 2 is the "occupied second-bath
# site", forced occupied by potPre=-bigEps).  Since that (1,1)-sector state
# is itself c_{0,up}|psi0>, psi0 must be the (2,1)-sector state obtained by
# adding the impurity's up electron back: up_word = 4 | (1<<0) = 5, dn_word
# unchanged = 4.  Verified directly below (c_op applied to (5,4) must land
# exactly on (4,4) with coefficient 1, matching the established seed).
psi0 = np.zeros(len(basis_N), dtype=complex)
psi0[index_N[(5, 4)]] = 1.0

_check = c_op.dot(psi0)
_nz = np.nonzero(_check)[0]
assert len(_nz) == 1 and basis_Nm1[_nz[0]] == (4, 4) and abs(_check[_nz[0]] - 1.0) < 1e-12, \
    "psi0=(5,4) must map to the established (4,4) N-1 seed under c_{0,up} -- if this fails, " \
    "the seed derivation above is stale and cross_check_gbek_propagation.py's own seed should " \
    "be re-verified against the current production code first."
print("psi0 = pure basis state (up=5,dn=4) in the (2,1) N-sector "
      "(verified: c_{0,up} psi0 = established (4,4) N-1 seed)")

# ---- Hamiltonian sequences (piecewise constant, one per dt step) ----
ts = [(n + 0.5) * dt for n in range(nsteps)]  # midpoint of each step, simplest convention
H_seq_N = [build_hamiltonian(nsites, basis_N, index_N, U, hop=[hop(t), hop(t)]) for t in ts]
H_seq_Nm1 = [build_hamiltonian(nsites, basis_Nm1, index_Nm1, U, hop=[hop(t), hop(t)]) for t in ts]

# ==================== 1) Python's construction ====================
phi_states = propagate(psi0, H_seq_N, dt)
G_py = compute_g_lesser(phi_states, c_op, H_seq_Nm1, dt)

# ==================== 2) C++-style construction ====================
# Seed ONCE at t=0, propagate the (N-1)-sector state FORWARD using only
# the (N-1)-sector Hamiltonian -- mirrors PsiHist[n] = U_{N-1}(t_n,0) PsiHist[0].
b0 = c_op.dot(psi0)
PsiHist = propagate(b0, H_seq_Nm1, dt)
Ntot = nsteps
G_cpp = np.zeros((Ntot + 1, Ntot + 1), dtype=complex)
for n in range(Ntot + 1):
    for k in range(n + 1):
        G_cpp[n, k] = 1j * np.vdot(PsiHist[k], PsiHist[n])

# ==================== 3) Brute-force ground truth ====================
# Fully independent: dense scipy.linalg.expm, no shared code with propagate()
# or compute_g_lesser().
H_N_dense = [H.toarray() for H in H_seq_N]
H_Nm1_dense = [H.toarray() for H in H_seq_Nm1]

phi_bf = [psi0]
for H in H_N_dense:
    phi_bf.append(expm(-1j * H * dt) @ phi_bf[-1])

b_bf = [c_op.dot(phi) for phi in phi_bf]

G_bf = np.zeros((Ntot + 1, Ntot + 1), dtype=complex)
for n in range(Ntot + 1):
    psi = b_bf[n]
    G_bf[n, n] = 1j * np.vdot(b_bf[n], psi)
    for k in range(n - 1, -1, -1):
        # backward step: undo H_Nm1_dense[k] (the propagator that advanced
        # index k -> k+1 going forward), i.e. apply exp(+i H dt)
        psi = expm(1j * H_Nm1_dense[k] * dt) @ psi
        G_bf[n, k] = 1j * np.vdot(b_bf[k], psi)

# ==================== Compare ====================
print()
print(f"{'n':>2} {'k':>2}  {'G_bruteforce':>28}  {'G_python':>28}  {'G_cpp-style':>28}  "
      f"{'|bf-py|':>10}  {'|bf-cpp|':>10}")
for n in range(1, Ntot + 1):
    for k in range(n + 1):
        bf, py, cpp = G_bf[n, k], G_py[n, k], G_cpp[n, k]
        print(f"{n:>2} {k:>2}  {bf.real:+.6f}{bf.imag:+.6f}j  "
              f"{py.real:+.6f}{py.imag:+.6f}j  "
              f"{cpp.real:+.6f}{cpp.imag:+.6f}j  "
              f"{abs(bf-py):>10.3e}  {abs(bf-cpp):>10.3e}")

max_err_py = np.max(np.abs(G_bf - G_py))
max_err_cpp = np.max(np.abs(G_bf - G_cpp))
print()
print(f"max|G_bruteforce - G_python|    = {max_err_py:.3e}  (expect: ~0, Python matches ground truth)")
print(f"max|G_bruteforce - G_cpp-style| = {max_err_cpp:.3e}  (expect: large, growing off-diagonal)")

# Diagonal-only comparison (should match well for BOTH, per the observed C++ behavior)
diag_err_cpp = np.max(np.abs(np.diag(G_bf) - np.diag(G_cpp)))
print(f"max diagonal |G_bf - G_cpp|     = {diag_err_cpp:.3e}  (expect: small, diagonal is spared)")
