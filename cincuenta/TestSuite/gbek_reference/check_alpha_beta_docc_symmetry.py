#!/usr/bin/env python3
"""
Empirically verify the simplifying claim used by the Fig. 7/8 double-occupation
plan: that d_alpha(t) = d_beta(t) identically, so the double occupation
d(t) = <n0up(t) n0dn(t)> can be computed directly from the existing alpha
trajectory (impurity seeded spin-up) without ever propagating a separate beta
trajectory (impurity seeded spin-down).

Why this should hold: gbek_selfconsistency.py's H_seq is built with
build_templates(), whose occ_term/empty_term sum hop_operator(...) over BOTH
spins already -- so H is symmetric under swapping the up and down channels.
The beta initial condition (impurity spin-down occupied, same doubly-occupied
bath pairs) lives in the mirror sector (nup=L, ndown=1+L instead of
nup=1+L, ndown=L) and is obtained from the alpha initial condition by relabeling
every basis state's (up_word, dn_word) as (dn_word, up_word). Since H is
spin-symmetric, propagating the relabeled state under the same H_seq gives the
relabeled trajectory, and n0up(t)*n0dn(t) is invariant under that relabeling
(it maps to n0dn(t)*n0up(t), the same operator) -- so d_beta(t) = d_alpha(t)
follows from H's spin symmetry, not from anything specific to U=2 or small L.
This script checks it directly rather than just asserting it, per the
project's "flag weak/circular tests" convention.

Method: run a few self-consistency iterations at small L, N to get a
converged V (bath coupling). Build H_seq for BOTH the alpha basis
(nup=1+L, ndown=L) and the beta basis (nup=L, ndown=1+L) from that same V,
propagate each sector's own initial state, and compare
d_alpha(t) = <n0up n0dn>_alpha(t) against d_beta(t) = <n0up n0dn>_beta(t).
"""
import numpy as np

from gbek_ed import build_basis, double_occupation_operator
from gbek_dynamics import propagate, compute_diag_observable
from gbek_cholesky import cholesky_causal
from gbek_selfconsistency import (
    build_templates, make_h_seq, initial_state, run_self_consistency,
)


def beta_initial_state(nsites, L, basis, index):
    """Mirror of initial_state(): impurity spin-down occupied instead of
    spin-up, same doubly-occupied bath pairs. Lives in the (nup=L, ndown=1+L)
    sector."""
    up_word = 0
    dn_word = 1  # impurity spin-down occupied
    for p in range(L):
        occ_site = 1 + 2 * p
        up_word |= (1 << occ_site)
        dn_word |= (1 << occ_site)
    psi0 = np.zeros(len(basis), dtype=complex)
    psi0[index[(up_word, dn_word)]] = 1.0
    return psi0


def main():
    L, N, dt, U, tstar_f, tq = 1, 10, 0.04, 2.0, 1.0, 0.25
    nsites = 1 + 2 * L

    print(f"Getting a converged V from run_self_consistency(L={L}, N={N}, U={U})...")
    Lambda, V, hop_t, ts, history, _ = run_self_consistency(
        L, N, dt, U, tstar_f, tq, n_iterations=10, tol=1e-8, verbose=False)
    print("convergence history:", [f"{h:.2e}" for h in history])

    # Alpha sector (nup=1+L, ndown=L) -- same as run_self_consistency's basis_N.
    basis_alpha, index_alpha = build_basis(nsites, 1 + L, L)
    Uterm_a, occ_a, empty_a = build_templates(nsites, basis_alpha, index_alpha, U, L)
    psi0_alpha = initial_state(nsites, L, basis_alpha, index_alpha)
    H_seq_alpha = make_h_seq(Uterm_a, occ_a, empty_a, V[:N, :])
    states_alpha = propagate(psi0_alpha, H_seq_alpha, dt)
    docc_op_alpha = double_occupation_operator(basis_alpha, site=0)
    d_alpha = compute_diag_observable(states_alpha, docc_op_alpha)

    # Beta sector (nup=L, ndown=1+L) -- the up<->down mirror sector.
    basis_beta, index_beta = build_basis(nsites, L, 1 + L)
    Uterm_b, occ_b, empty_b = build_templates(nsites, basis_beta, index_beta, U, L)
    psi0_beta = beta_initial_state(nsites, L, basis_beta, index_beta)
    H_seq_beta = make_h_seq(Uterm_b, occ_b, empty_b, V[:N, :])
    states_beta = propagate(psi0_beta, H_seq_beta, dt)
    docc_op_beta = double_occupation_operator(basis_beta, site=0)
    d_beta = compute_diag_observable(states_beta, docc_op_beta)

    err = np.max(np.abs(d_alpha - d_beta))
    print(f"\nmax|d_alpha(t) - d_beta(t)| over t in [0,{N*dt:.2f}]: {err:.3e}")
    print(f"d_alpha(t=0) = {d_alpha[0]:.6f} (expect 0.0, impurity singly occupied)")
    print(f"d_alpha(t=tmax) = {d_alpha[-1]:.6f}")

    assert d_alpha[0] == 0.0, "impurity is singly occupied at t=0, d(0) must be exactly 0"
    assert err < 1e-10, (
        "d_alpha(t) != d_beta(t): the plan's no-averaging-needed simplification is WRONG, "
        "an explicit beta-trajectory propagation is required after all")
    print("\nOK: d_alpha(t) == d_beta(t) to numerical precision -- "
          "no beta-trajectory propagation is needed for double occupation.")


if __name__ == "__main__":
    main()
