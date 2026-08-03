"""
Hamiltonian construction and real-time propagation for the star-geometry
SIAM, plus two-time Green's function extraction.  Built on gbek_ed.py.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from gbek_ed import annihilate, create, build_basis, c_operator, number_operator


def hop_operator(basis, index, site_a, site_b, spin):
    """
    Sparse matrix (within a single fixed-(nup,ndown) sector) for
    c^dag_{site_a,spin} c_{site_b,spin}: annihilate at site_b, then create at
    site_a, tracking the combined Jordan-Wigner sign, staying in-sector since
    particle number is conserved overall.
    """
    rows, cols, vals = [], [], []
    for j, (u, d) in enumerate(basis):
        word = u if spin == 0 else d
        res = annihilate(word, site_b)
        if res is None:
            continue
        mid_word, sign1 = res
        res2 = create(mid_word, site_a)
        if res2 is None:
            continue
        new_word, sign2 = res2
        new_state = (new_word, d) if spin == 0 else (u, new_word)
        i = index.get(new_state)
        if i is None:
            continue
        rows.append(i)
        cols.append(j)
        vals.append(sign1 * sign2)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)))


def build_hamiltonian(nsites, basis, index, U, hop):
    """
    H = U (n0up-1/2)(n0dn-1/2) + sum_{p=1..nsites-1, spin} hop[p-1] * (c0^dag c_p + h.c.)
    hop: array of length nsites-1, real hopping amplitude impurity<->bath site p.
    """
    n0up_diag = number_operator(basis, 0, 0).diagonal()
    n0dn_diag = number_operator(basis, 0, 1).diagonal()
    H = sp.diags(U * (n0up_diag - 0.5) * (n0dn_diag - 0.5)).tocsr().astype(complex)
    for p in range(1, nsites):
        v = hop[p - 1]
        if v == 0:
            continue
        for spin in (0, 1):
            term = hop_operator(basis, index, 0, p, spin)
            H = H + v * (term + term.conj().T)
    return H


def expm_multiply_taylor(H, psi, coeff, tol=1e-13, max_terms=60):
    """
    Compute exp(coeff*H) @ psi via a truncated Taylor series.  For our use
    case (small dt * ||H||, genuinely sparse H) this converges in ~8-15
    terms and is much cheaper than scipy's general-purpose expm_multiply,
    which pays for adaptive norm estimation and scaling/squaring it doesn't
    need here.

    Occasionally a Cholesky "optimal update" step (see gbek_cholesky.py, the
    n > L branch) produces an ill-conditioned/large bath coupling, making
    dt*||H|| too large for a low-order Taylor series -- the same failure
    mode documented in NeqBathDecomposition.h's history (spikes from a
    near-zero pivot).  Detect that (term growing instead of shrinking, or
    not converged within max_terms) and fall back to scipy's robust
    (scaling-and-squaring) expm_multiply for just that one call.
    """
    term = psi.astype(complex, copy=True)
    result = term.copy()
    prev_norm = np.linalg.norm(term)
    for k in range(1, max_terms + 1):
        term = (coeff / k) * (H @ term)
        term_norm = np.linalg.norm(term)
        if not np.isfinite(term_norm) or (k > 5 and term_norm > prev_norm):
            return spla.expm_multiply(coeff * H, psi)
        result = result + term
        if term_norm < tol * (np.linalg.norm(result) + 1e-300):
            return result
        prev_norm = term_norm
    return spla.expm_multiply(coeff * H, psi)


def propagate(state, H_seq, dt):
    """
    Propagate `state` forward through a sequence of (piecewise-constant)
    Hamiltonians H_seq, each held for time dt.  Returns list of states
    [state(t0), state(t0+dt), ...] of length len(H_seq)+1.
    """
    states = [state]
    psi = state
    for H in H_seq:
        psi = expm_multiply_taylor(H, psi, -1j * dt)
        states.append(psi)
    return states


def propagate_backward(state, H_seq, dt):
    """
    Propagate `state` backward through H_seq (same convention as propagate,
    but using +i dt): if `state` lives at time index k (after H_seq[:k] applied
    going forward), this undoes H_seq[k-1], H_seq[k-2], ... in reverse.
    Returns list [state(t_k), state(t_{k-1}), ...].
    """
    states = [state]
    psi = state
    for H in reversed(H_seq):
        psi = expm_multiply_taylor(H, psi, 1j * dt)
        states.append(psi)
    return states


def compute_g_lesser(phi_states, c_op, H_seq_minus, dt):
    """
    phi_states: [|phi(t_0)>, ..., |phi(t_N)>] forward-propagated reference
      states in the N-particle sector.
    c_op: sparse matrix for c_{0,spin}, mapping N-sector -> (N-1)-sector.
    H_seq_minus: Hamiltonian sequence for the (N-1)-sector, same dt grid.

    Returns G[n, j] = i <c^dag_0(t_j) c_0(t_n)> for j <= n (else 0).
    """
    N = len(phi_states) - 1
    b_states = [c_op.dot(phi) for phi in phi_states]
    G = np.zeros((N + 1, N + 1), dtype=complex)
    for n in range(N + 1):
        psi = b_states[n]
        G[n, n] = 1j * np.vdot(b_states[n], psi)
        for k in range(n - 1, -1, -1):
            psi = expm_multiply_taylor(H_seq_minus[k], psi, 1j * dt)
            G[n, k] = 1j * np.vdot(b_states[k], psi)
    return G


def compute_diag_observable(states, diag_op):
    """states: list of state vectors (e.g. phi_states from propagate()).
    diag_op: sparse diagonal operator (e.g. from double_occupation_operator).
    Returns real array, one expectation value per state."""
    diag = diag_op.diagonal()
    return np.array([np.vdot(psi, diag * psi).real for psi in states])


if __name__ == "__main__":
    # Validation: U=0, static hopping V between 2 sites (impurity + 1 bath),
    # single spin-up electron.  Exact single-particle result: with H_1p =
    # [[0, V], [V, 0]] in the {imp, bath} basis, starting in |imp>, the
    # single-particle amplitude is c_imp(t) = cos(V t), c_bath(t) = -i sin(V t).
    # G^<(t,t') = i <c^dag(t') c(t)> for a single fixed particle (no dn spin)
    # reduces to the single-particle propagator overlap.
    nsites = 2
    V = 0.7
    dt = 0.02
    nsteps = 50

    basis, index = build_basis(nsites, 1, 0)
    H = build_hamiltonian(nsites, basis, index, U=0.0, hop=[V])
    H_seq = [H] * nsteps

    # initial state: |imp> occupied, i.e. up_word = 0b01 (site 0 set)
    psi0 = np.zeros(len(basis), dtype=complex)
    psi0[index[(1, 0)]] = 1.0

    states = propagate(psi0, H_seq, dt)

    # amplitude to be back on impurity at time t_n = n*dt
    imp_state_index = index[(1, 0)]
    amp_numeric = np.array([s[imp_state_index] for s in states])
    ts = np.arange(nsteps + 1) * dt
    amp_exact = np.cos(V * ts)

    err = np.max(np.abs(amp_numeric - amp_exact))
    print("max |numeric - exact| amplitude error:", err)
    assert err < 1e-8, "propagation does not match analytic single-particle result"
    print("OK: propagation matches analytic single-particle oscillation")

    # Validate compute_g_lesser: for a single particle in a pure state,
    # G^<(t,t') = i * amp_imp(t) * conj(amp_imp(t')) exactly.
    basis_m1, index_m1 = build_basis(nsites, 0, 0)
    c0up = c_operator(nsites, basis, index, basis_m1, index_m1, site=0, spin=0)
    H_seq_minus = [build_hamiltonian(nsites, basis_m1, index_m1, U=0.0, hop=[V])] * nsteps

    G = compute_g_lesser(states, c0up, H_seq_minus, dt)
    G_exact = 1j * np.outer(amp_numeric, np.conj(amp_numeric))
    mask = np.tril(np.ones_like(G, dtype=bool))  # only j<=n populated
    err_g = np.max(np.abs(G[mask] - G_exact[mask]))
    print("max |G_numeric - G_exact| (lower triangle):", err_g)
    assert err_g < 1e-8, "two-time G^< does not match analytic single-particle result"
    print("OK: two-time G^< matches analytic single-particle formula")

    # Validate compute_diag_observable / double_occupation_operator: this
    # sector has no down-spin electron at all (nup=1, ndown=0), so d(t) must
    # be exactly 0 for every t -- a minimal smoke test before wiring the
    # observable into the full self-consistency loop (see also
    # check_alpha_beta_docc_symmetry.py for the alpha/beta symmetry check,
    # which needs the self-consistency machinery and lives in its own script).
    from gbek_ed import double_occupation_operator
    docc_op = double_occupation_operator(basis, site=0)
    d_t = compute_diag_observable(states, docc_op)
    print("max |d(t)| (expect exactly 0, no down-spin electron in this sector):",
          np.max(np.abs(d_t)))
    assert np.max(np.abs(d_t)) == 0.0, "double occupation must be exactly 0 with no dn electron"
    print("OK: double_occupation_operator / compute_diag_observable smoke test passed")
