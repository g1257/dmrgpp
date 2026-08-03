#!/usr/bin/env python3
"""
Test whether gbek_selfconsistency.py's single-trajectory shortcut
(Gbeta_0sigma = Galpha_0,-sigma, asserted, not computed) is actually
equivalent to running a REAL, independently-propagated beta trajectory, as
GBEK's paper itself literally prescribes (Sec. VI.A, Eq. 70):

    "we average over two Green functions Ga and Gb, where the impurity of
    system alpha (beta) is populated initially by a single up-spin
    (down-spin) electron... Gsigma(t,t') = 1/2(Ga_0sigma + Gb_0sigma).
    Taking the average restores particle-hole symmetry, which is not given
    for Ga or Gb alone."

This is the paper's own architecture: TWO independently-seeded systems,
BOTH measured with the SAME spin operator sigma=up, then averaged -- not
one trajectory measured with two different operators (up and down). This
script builds the real beta system (impurity seeded spin-down instead of
spin-up, same bath V(t)) and propagates it independently, so Gbeta_0up can
be directly compared against:
  (a) gbek_selfconsistency.py's shortcut proxy G_dn (from the alpha
      trajectory alone, using the down-spin operator -- what the module's
      docstring assumes equals Gbeta_0up by symmetry, but never checks), and
  (b) cincuenta's C++ sectorBeta_ diagonal (see ImpuritySolverNeqGBEK.h),
      which is a real, independently-propagated second sector exactly like
      this script's beta trajectory.

Both trajectories share the SAME rank-L Cholesky bath coupling V(t) --
obtained from a converged run_self_consistency() call using the standard
alpha-only shortcut -- since V(t) is a property of the (converged) input
Weiss field, not of which trajectory is being used to probe it.

Usage:
    uv run --with numpy --with scipy python3 cross_check_beta_trajectory.py \\
        --L 3 --N 100 --dt 0.04 --U 2.0 --tq 0.25 --iterations 30 --tol 1e-8
"""
import argparse

import numpy as np

from gbek_ed import build_basis, c_operator, double_occupation_operator
from gbek_dynamics import propagate, compute_g_lesser
from gbek_selfconsistency import (
    run_self_consistency, build_templates, make_h_seq, v_ramp,
)


def initial_state_beta(nsites, L, basis, index):
    """Mirror of gbek_selfconsistency.initial_state(): impurity's single
    extra electron is spin-DOWN instead of spin-up. Bath occ/empty
    structure (doubly-occupied vs. empty pairs) is unchanged -- it's
    spin-symmetric by construction, so both up_word and dn_word get the
    same occ_site bits regardless of which system (alpha/beta) this is."""
    up_word = 0
    dn_word = 1  # impurity spin-down occupied (mirror of alpha's up_word=1)
    for p in range(L):
        occ_site = 1 + 2 * p
        up_word |= (1 << occ_site)
        dn_word |= (1 << occ_site)
    psi0 = np.zeros(len(basis), dtype=complex)
    psi0[index[(up_word, dn_word)]] = 1.0
    return psi0


def compute_beta_g_up(L, N, dt, U, V):
    """Build and propagate the REAL beta system (nup_ext=L, ndown_ext=1+L,
    impurity seeded spin-down) and measure G_up(t,t') on it -- i.e. Gbeta_0up
    of GBEK Eq. 70, computed exactly as prescribed (independent trajectory,
    same operator), not asserted via a symmetry shortcut."""
    nsites = 1 + 2 * L
    nup_ext, ndown_ext = L, 1 + L  # mirror of alpha's (1+L, L)
    basis_N, index_N = build_basis(nsites, nup_ext, ndown_ext)
    basis_up, index_up = build_basis(nsites, nup_ext - 1, ndown_ext)

    Uterm_N, occ_N, empty_N = build_templates(nsites, basis_N, index_N, U, L)
    Uterm_up, occ_up, empty_up = build_templates(nsites, basis_up, index_up, U, L)

    c0up = c_operator(nsites, basis_N, index_N, basis_up, index_up, site=0, spin=0)

    psi0 = initial_state_beta(nsites, L, basis_N, index_N)

    H_seq_N = make_h_seq(Uterm_N, occ_N, empty_N, V[:N, :])
    H_seq_up = make_h_seq(Uterm_up, occ_up, empty_up, V[:N, :])

    phi_states = propagate(psi0, H_seq_N, dt)
    G_up_beta = compute_g_lesser(phi_states, c0up, H_seq_up, dt)
    return G_up_beta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--L", type=int, default=3)
    ap.add_argument("--N", type=int, default=100)
    ap.add_argument("--dt", type=float, default=0.04)
    ap.add_argument("--U", type=float, default=2.0)
    ap.add_argument("--tstar-f", type=float, default=1.0)
    ap.add_argument("--tq", type=float, default=0.25)
    ap.add_argument("--iterations", type=int, default=30)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--check-times", default="1,2,3,4",
                     help="comma-separated t values to report (must be multiples of dt)")
    args = ap.parse_args()

    L, N, dt = args.L, args.N, args.dt

    print("Step 1: converge the alpha-only shortcut self-consistency loop "
          "(same as gbek-atomic-limit-exact-lesser) to get V(t)...")
    Lambda, V, hop_t, ts, history, _ = run_self_consistency(
        L, N, dt, args.U, args.tstar_f, args.tq,
        n_iterations=args.iterations, tol=args.tol, verbose=False)
    print(f"  converged after {len(history)} iterations, final diff={history[-1]:.2e}")

    # Recompute Galpha_0up on this converged V (mirrors what
    # run_self_consistency does internally on its last iteration).
    nsites = 1 + 2 * L
    nup_ext, ndown_ext = 1 + L, L
    basis_N, index_N = build_basis(nsites, nup_ext, ndown_ext)
    basis_up, index_up = build_basis(nsites, nup_ext - 1, ndown_ext)
    basis_dn, index_dn = build_basis(nsites, nup_ext, ndown_ext - 1)
    Uterm_N, occ_N, empty_N = build_templates(nsites, basis_N, index_N, args.U, L)
    Uterm_up, occ_up, empty_up = build_templates(nsites, basis_up, index_up, args.U, L)
    Uterm_dn, occ_dn, empty_dn = build_templates(nsites, basis_dn, index_dn, args.U, L)
    c0up = c_operator(nsites, basis_N, index_N, basis_up, index_up, site=0, spin=0)
    c0dn = c_operator(nsites, basis_N, index_N, basis_dn, index_dn, site=0, spin=1)
    from gbek_selfconsistency import initial_state
    psi0 = initial_state(nsites, L, basis_N, index_N)
    H_seq_N = make_h_seq(Uterm_N, occ_N, empty_N, V[:N, :])
    H_seq_up = make_h_seq(Uterm_up, occ_up, empty_up, V[:N, :])
    H_seq_dn = make_h_seq(Uterm_dn, occ_dn, empty_dn, V[:N, :])
    phi_states = propagate(psi0, H_seq_N, dt)
    G_up_alpha = compute_g_lesser(phi_states, c0up, H_seq_up, dt)
    G_dn_alpha = compute_g_lesser(phi_states, c0dn, H_seq_dn, dt)  # the shortcut proxy

    print("\nStep 2: build and propagate the REAL independent beta trajectory "
          "(impurity seeded spin-down), same V(t), measure G_up on it...")
    G_up_beta = compute_beta_g_up(L, N, dt, args.U, V)

    check_ns = [int(round(float(t) / dt)) for t in args.check_times.split(",")]
    print(f"\n{'t':>5} {'Galpha_up(t,t)':>16} {'shortcut Gdn_alpha':>20} "
          f"{'REAL Gbeta_up(t,t)':>20} {'shortcut vs real diff':>22} "
          f"{'0.5(Ga+Gb) [real]':>20} {'0.5(Ga+Gb) [shortcut]':>22}")
    for n in check_ns:
        t = n * dt
        ga = G_up_alpha[n, n].imag
        gd_shortcut = G_dn_alpha[n, n].imag
        gb_real = G_up_beta[n, n].imag
        diff = gd_shortcut - gb_real
        avg_real = 0.5 * (ga + gb_real)
        avg_shortcut = 0.5 * (ga + gd_shortcut)
        print(f"{t:5.2f} {ga:16.6f} {gd_shortcut:20.6f} {gb_real:20.6f} "
              f"{diff:22.6f} {avg_real:20.6f} {avg_shortcut:22.6f}")

    print("\nInterpretation:")
    print("  'shortcut vs real diff' ~ 0  => the module's symmetry assumption "
          "(Gbeta_0sigma = Galpha_0,-sigma) holds numerically; the shortcut is "
          "valid and Python's exact-0.5 diagonal is real physics, not an artifact.")
    print("  '0.5(Ga+Gb) [real]' close to 0.5  => the TRUE paper-Eq.-70 average, "
          "built the way GBEK's own reference implementation and cincuenta's C++ "
          "both do it (independent trajectories, same operator), also restores "
          "particle-hole symmetry -- meaning C++'s drift is a bug in its beta-sector "
          "construction, not an inherent property of the independent-trajectory method.")


if __name__ == "__main__":
    main()
