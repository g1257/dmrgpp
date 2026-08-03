"""
Independent, from-scratch reconstruction of a small GBEK extended-Fock-space
Hamiltonian, used to derive the expected eigenvalue spectrum hardcoded in
cincuenta/src/tests/test_ImpuritySolverNeqGBEK.cpp's
"Hamiltonian eigenvalue spectrum matches an independent Python reconstruction"
test. That C++ test exists specifically to catch Hamiltonian-construction
bugs (wrong Jordan-Wigner sign, wrong site layout, wrong occ/empty coupling
polarity, wrong diagonal/potential terms) that row-by-row V comparisons or
end-to-end self-consistency comparisons cannot cleanly localize, because V
and the converged Weiss field both depend on the self-consistency loop's own
history and are not themselves gauge/basis invariant. Eigenvalues of a
FIXED, externally-supplied-V Hamiltonian are: this script and the C++ test
build the identical physical system independently (this script never calls
into cincuenta or LanczosPlusPlus) and compare spectra, which are
basis-ordering-independent -- so a mismatch can only come from the physics
of the Hamiltonian itself, not from a basis-labeling difference between the
two languages.

System: nBath=0 (true atomic limit), L=1 second-bath pair, nup=1/ndown=0
spin-imbalanced seed (this session's production setup). Reconstructs the N-1
sector used for G^< at the impurity: if the N sector is (nupExt, ndownExt) =
(nup+L, ndown+L) = (2, 1), the N-1 sector (one impurity up-electron removed)
is (nupExt-1, ndownExt) = (1, 1) -- NOT (nupExt, ndownExt) itself. (An
earlier draft of this script used (2, 1) for the N-1 sector by mistake and
got a different, wrong trace; re-derive this from buildSector's own
`model.createBasis(nupExt-1, ndownExt)` call if it's ever unclear again.)

Site layout matches cincuenta's own convention exactly (site 0 = impurity,
site 1 = "empty" second-bath partner, site 2 = "occupied" second-bath
partner -- see ImpuritySolverNeqGBEK.h's buildSector: emptyS = 1+nBath+p,
occS = 1+nBath+L+p, both counted from 0).

Coupling convention: site 1 (empty) and site 2 (occ) use OPPOSITE V/conj(V)
polarity, matching gbek_selfconsistency.py's build_templates() and GBEK
Eq. 49/52a (-iLambda^<_+ couples the occupied orbital via V(t)V(t')^*, while
+iLambda^>_+ = (-iLambda^<_+)^* couples the empty orbital). This was, for a
while, a live, not-yet-resolved question -- an earlier version of this
script instead matched cincuenta's then-current (symmetric, WRONG)
buildHextCSR convention, which gave both orbitals the SAME polarity. That
was resolved 2026-07-14: an independently-propagated real beta trajectory
(cross_check_beta_trajectory.py) proved the asymmetric convention is the
one that restores the exact particle-hole-symmetric sum
Galpha_up(t,t)+Gbeta_up(t,t)=1 at every timestep, for any bath V(t), while
the symmetric convention (both this script's and cincuenta's old code) does
not. ImpuritySolverNeqGBEK.h's buildHextCSR/applyHext were fixed to match.

Run standalone to (re)generate the eigenvalues the C++ test hardcodes:
    uv run --with numpy --with scipy python3 cross_check_gbek_hamiltonian.py
"""
import numpy as np
import sys
sys.path.insert(0, '.')
from gbek_ed import build_basis, number_operator
from gbek_dynamics import hop_operator

U = 2.0
potPost = [-U / 2, 0.0, 0.0]   # site 0 (impurity), site 1 (empty), site 2 (occ)
vMid = 0.3 + 0.2j              # matches the C++ test's fixed vMid

nsites = 3
nupExt, ndownExt = 1, 1        # N-1 sector: (nupExt_N - 1, ndownExt_N) = (2-1, 1)
basis, index = build_basis(nsites, nupExt, ndownExt)
dim = len(basis)
print(f"dim={dim}")

n = [[number_operator(basis, site, spin).diagonal() for spin in (0, 1)]
     for site in range(nsites)]

diag = (U * n[0][0] * n[0][1]
        + sum(potPost[site] * (n[site][0] + n[site][1]) for site in range(nsites)))

H = np.diag(diag).astype(complex)

# Second bath (L=1, p=0): impurity <-> empty (site 1), impurity <-> occ (site 2).
# hop_operator(basis, index, site_a, site_b, spin) is c^dagger_{site_a,spin}
# c_{site_b,spin} (hop FROM site_b TO site_a), matching cincuenta's
# tryHop(from=site_b, to=site_a, ...) convention exactly.
for spin in (0, 1):
    # impurity -> empty (site 1): plain vMid.  empty -> impurity: conj(vMid).
    H += vMid * hop_operator(basis, index, 1, 0, spin).toarray()
    H += np.conj(vMid) * hop_operator(basis, index, 0, 1, spin).toarray()
    # occ (site 2) has the OPPOSITE polarity from empty's:
    # impurity -> occ: conj(vMid).  occ -> impurity: plain vMid.
    H += np.conj(vMid) * hop_operator(basis, index, 2, 0, spin).toarray()
    H += vMid * hop_operator(basis, index, 0, 2, spin).toarray()

herm_err = np.max(np.abs(H - H.conj().T))
print(f"max|H - H^dagger| = {herm_err:.3e}")

eigs = np.sort(np.linalg.eigvalsh(H))
print("eigenvalues:", " ".join(f"{e:.14f}" for e in eigs))
