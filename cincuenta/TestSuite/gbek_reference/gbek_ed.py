"""
Standalone exact-diagonalization building blocks for a star-geometry SIAM
(1 impurity + Lbath bath sites), independent of cincuenta/LanczosPlusPlus.

Sites are numbered 0..nsites-1 with site 0 the impurity.  Basis states for a
fixed (nup, ndown) sector are represented as (up_word, dn_word) bitmasks.
"""
import itertools
import numpy as np
import scipy.sparse as sp


def combos(nsites, n):
    """All bitmasks over nsites bits with exactly n bits set."""
    words = []
    for bits in itertools.combinations(range(nsites), n):
        w = 0
        for b in bits:
            w |= (1 << b)
        words.append(w)
    return words


def build_basis(nsites, nup, ndown):
    up_words = combos(nsites, nup)
    dn_words = combos(nsites, ndown)
    basis = [(u, d) for u in up_words for d in dn_words]
    index = {state: i for i, state in enumerate(basis)}
    return basis, index


def jw_sign(word, site):
    """Fermion sign from anticommuting c_site past occupied sites < site."""
    mask = (1 << site) - 1
    return -1 if bin(word & mask).count("1") % 2 else 1


def annihilate(word, site):
    """Return (new_word, sign) for c_site acting on `word`, or None if empty."""
    if not (word >> site) & 1:
        return None
    sign = jw_sign(word, site)
    return word & ~(1 << site), sign


def create(word, site):
    """Return (new_word, sign) for c_site^dagger acting on `word`, or None if occupied."""
    if (word >> site) & 1:
        return None
    sign = jw_sign(word, site)
    return word | (1 << site), sign


def c_operator(nsites, basis_from, index_from, basis_to, index_to, site, spin):
    """
    Sparse matrix for c_{site,spin}: basis_from (dim Nf) -> basis_to (dim Nt).
    spin: 0=up, 1=down.
    """
    rows, cols, vals = [], [], []
    for j, (u, d) in enumerate(basis_from):
        word = u if spin == 0 else d
        res = annihilate(word, site)
        if res is None:
            continue
        new_word, sign = res
        new_state = (new_word, d) if spin == 0 else (u, new_word)
        i = index_to.get(new_state)
        if i is None:
            continue
        rows.append(i)
        cols.append(j)
        vals.append(sign)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(basis_to), len(basis_from)))


def cdag_operator(nsites, basis_from, index_from, basis_to, index_to, site, spin):
    rows, cols, vals = [], [], []
    for j, (u, d) in enumerate(basis_from):
        word = u if spin == 0 else d
        res = create(word, site)
        if res is None:
            continue
        new_word, sign = res
        new_state = (new_word, d) if spin == 0 else (u, new_word)
        i = index_to.get(new_state)
        if i is None:
            continue
        rows.append(i)
        cols.append(j)
        vals.append(sign)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(basis_to), len(basis_from)))


def number_operator(basis, site, spin):
    n = np.zeros(len(basis))
    for j, (u, d) in enumerate(basis):
        word = u if spin == 0 else d
        n[j] = (word >> site) & 1
    return sp.diags(n)


def double_occupation_operator(basis, site=0):
    """Diagonal sparse operator n_{site,up} n_{site,down}."""
    n_up = number_operator(basis, site, 0)
    n_dn = number_operator(basis, site, 1)
    return sp.diags(n_up.diagonal() * n_dn.diagonal())


if __name__ == "__main__":
    # Smoke test: isolated Hubbard atom (nsites=1), particle-hole symmetric,
    # mu = -U/2.  Ground states at half filling (nup=1,ndown=0 or 0,1) should
    # have energy -U/2 (single occupied orbital, no double-occupancy cost).
    U = 2.0
    nsites = 1
    basis, index = build_basis(nsites, 1, 0)
    print("dim (nup=1,ndown=0):", len(basis))
    n0up = number_operator(basis, 0, 0).toarray()
    n0dn = number_operator(basis, 0, 1).toarray()
    H = U * (n0up - 0.5).dot(n0dn - 0.5)
    print("H =", H, " (expect -U/2 =", -U / 2, ")")
