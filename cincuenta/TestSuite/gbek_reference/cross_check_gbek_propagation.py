"""
Independent, from-scratch cross-check of cincuenta's GBEK time propagation
(Krylov exponential + seeded initial state), used to derive the expected
G^< values hardcoded in test_ImpuritySolverNeqGBEK.cpp's
"time propagation matches an independent Python reconstruction" test.

cross_check_gbek_hamiltonian.py (see that script's docstring) already
confirmed the STATIC Hamiltonian, for one fixed V, is built correctly. This
script extends the same from-scratch construction to a two-step propagation
with two DIFFERENT, externally-fixed vMid values (never touching
NeqBathDecomposition or the self-consistency loop), and to the seeded
initial state itself -- the two remaining pieces that script didn't cover.

Initial state: cincuenta's seedState() builds the N-sector (nupExt=2,
ndownExt=1) ground state of a "pre-quench" Hamiltonian where the second
bath sites are pinned by a huge (+-500) potential -- empty site forced
unoccupied, occupied site forced doubly occupied -- so the N-sector ground
state is, to floating-point precision, the single product state with the
impurity carrying the one remaining up-electron. Removing that electron
(c_{imp,up}) lands in the N-1 sector's single basis state with BOTH second
bath sites in their initial configuration and the impurity/empty site both
empty. Confirmed directly against cincuenta's own PsiHist[0]: it is a pure
basis state (amplitude 1, all others 0) at the state with up_word=dn_word=4
(binary 100 -- site 2, "occ", bit set) in the nup=1/ndown=1, nsites=3 basis
built the same way as cross_check_gbek_hamiltonian.py's Hamiltonian check.

Run standalone to (re)generate the values the C++ test hardcodes:
    uv run --with numpy --with scipy python3 cross_check_gbek_propagation.py
"""
import numpy as np
from scipy.linalg import expm
import sys
sys.path.insert(0, '.')
from gbek_ed import build_basis, number_operator
from gbek_dynamics import hop_operator

U = 2.0
potPost = [-U / 2, 0.0, 0.0]   # site 0 (impurity), site 1 (empty), site 2 (occ)
dt = 0.05
vMid_step0 = 0.3 + 0.2j         # propagates t=0 -> t=dt
vMid_step1 = 0.1 - 0.15j        # propagates t=dt -> t=2*dt

nsites = 3
nupExt, ndownExt = 1, 1        # N-1 sector, same as cross_check_gbek_hamiltonian.py
basis, index = build_basis(nsites, nupExt, ndownExt)
dim = len(basis)

n = [[number_operator(basis, site, spin).diagonal() for spin in (0, 1)]
     for site in range(nsites)]
diag = (U * n[0][0] * n[0][1]
        + sum(potPost[site] * (n[site][0] + n[site][1]) for site in range(nsites)))


def build_H(vMid):
    """Same construction as cross_check_gbek_hamiltonian.py, parameterized by vMid."""
    H = np.diag(diag).astype(complex)
    for spin in (0, 1):
        H += vMid * hop_operator(basis, index, 1, 0, spin).toarray()
        H += np.conj(vMid) * hop_operator(basis, index, 0, 1, spin).toarray()
        # Occupied orbital (site 2) has the OPPOSITE conjugation from the
        # empty orbital's (site 1) -- GBEK Eq. 49/52a: -iLambda^<_+ couples
        # the occupied orbital via V(t)V(t')^*, while +iLambda^>_+ =
        # (-iLambda^<_+)^* couples the empty orbital. See
        # gbek_selfconsistency.py::build_templates's docstring (validated
        # directly against the paper) and ImpuritySolverNeqGBEK.h's matching
        # fix. An earlier version of this script used the same (wrong,
        # symmetric) convention for both orbitals, matching a real bug that
        # was also present in the C++ at the time.
        H += np.conj(vMid) * hop_operator(basis, index, 2, 0, spin).toarray()
        H += vMid * hop_operator(basis, index, 0, 2, spin).toarray()
    return H


psi0 = np.zeros(dim, dtype=complex)
psi0[index[(4, 4)]] = 1.0

H_a = build_H(vMid_step0)
H_b = build_H(vMid_step1)

psi1 = expm(-1j * H_a * dt) @ psi0
psi2 = expm(-1j * H_b * dt) @ psi1
psis = [psi0, psi1, psi2]

print("norms:", [f"{np.vdot(p, p).real:.14f}" for p in psis])

print("G_lesser(n,j) = i * conj(psi_j) . psi_n:")
for nrow in range(3):
    for j in range(nrow + 1):
        g = 1j * np.vdot(psis[j], psis[nrow])
        print(f"  n={nrow} j={j}: ({g.real:.14f},{g.imag:.14f})")
