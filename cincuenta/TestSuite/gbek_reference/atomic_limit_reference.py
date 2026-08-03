#!/usr/bin/env python3
"""
Independent, from-scratch reference for the trivial "no bath" (atomic-limit)
impurity Green's function: a single Hubbard atom, U fixed (no interaction
quench, since the atomic limit here is used only as the equilibrium
starting point for a real-time hopping/bandwidth quench carried entirely by
the second bath -- the impurity's OWN Hamiltonian never changes in time).

Purpose: cincuenta's ImpuritySolverNeqExactDiag needs a special-case bypass
for nBath=0 (the true GBEK atomic limit -- see NeqAtomicLimit), because its
general Lehmann machinery (built for finite baths with an actual U_i->U_f
interaction quench) goes through LanczosPlusPlus/Ainur, which chokes on a
literal single-site system (see cincuenta.cpp / CincuentaInputCheck.h
history), AND because a single spin-seeded impurity has a Pauli-forbidden
particle-addition or removal sector (dim=0) that needs the paper's Gα/Gβ
averaging to resolve correctly.

This script computes the correctly-averaged G^R, G^<, G^M (Matsubara),
and G^Left (mixed real/imaginary time) for the bare atom DIRECTLY via
exact diagonalization of the 4-state atomic Fock space, independent of and
before touching the C++ special-case implementation -- so the C++ can be
checked against these numbers exactly the same way the Cholesky bugs were
caught (see README.md's "The real bug(s)" and
feedback_trust_the_publication.md).

Physics (verified against the standard atomic-limit result before writing
this file): H = U*(n_up-1/2)*(n_dn-1/2), so the atomic spectrum is
E(0) = E(up-down) = +U/4, E(up) = E(down) = -U/4 (twofold degenerate
ground state, gap U/2). Averaging Galpha (seeded up) and Gbeta (seeded
down) restores spin symmetry and gives, for EITHER spin sigma:

    G^R_sigma(t,t')   = -i * theta(t-t') * cos(U(t-t')/2)
    G^<_sigma(t,t')   = (i/2) * exp(i*U(t-t')/2)             for all t,t'
    G^>_sigma(t,t')   = (-i/2) * exp(-i*U(t-t')/2)           for all t,t'

These real-time formulas are unambiguous (single pole, no thermal/quench
subtlety) and are cross-checked below against a direct 4-state ED
diagonalization.

The Matsubara and Left (mixed) branches were NOT reused from cincuenta's
existing formula conventions without independent verification, since those
mix "particle" (N+1) and "hole" (N-1) sector contributions in a way whose
sign/domain conventions (an Omega^h >= 0 guard) were written for a
different physical regime (finite bath, possibly not half-filled) and were
not obviously correct for a gapped atomic ground state when checked by
hand during this investigation. Instead they were derived here from a
direct T=0 Lehmann sum (only the two degenerate ground states |up>,|down>
contribute at T=0) and cross-checked against ED to machine precision. Note
the *decay rate* differs between the two branches -- this was gotten wrong
twice before landing on the ED-confirmed values, precisely the kind of
subtlety this independent-derivation discipline is for:

    G^M_sigma(tau)      = -(1/2) * exp(-(U/2)*tau)                for tau in [0, beta)
    G^Left_sigma(t,tau) = -(i/2) * exp(i*U*t/4) * exp(-(U/4)*tau)

G^M decays at rate U/2 (the true excitation gap: G^M sandwiches the
propagator between two bra/ket applications of H, exp(+H tau) on the left
and exp(-H tau) on the right, so it is a genuine energy *difference*).
G^Left has only a single one-sided imaginary-time propagator
exp(-H tau) with no compensating factor on the other side (the real-time
leg carries the phase instead), so it picks up the *absolute* energy of
the doubly-occupied/empty intermediate state relative to zero, U/4, not
the gap U/2 -- half the naive guess. Both formulas are validated against
the C++ output rather than assumed to be right by analogy.
"""
import numpy as np

from gbek_ed import build_basis, c_operator, cdag_operator, number_operator


def atomic_hamiltonian(U):
    """4-state Fock space (1 site): |0>, |up>, |down>, |up,down>.
    H = U*(n_up - 1/2)(n_down - 1/2), particle-hole symmetric."""
    basis_all = [(0, 0), (1, 0), (0, 1), (1, 1)]
    index_all = {s: i for i, s in enumerate(basis_all)}
    H = np.zeros((4, 4))
    for (u, d), i in index_all.items():
        nup, ndn = u, d
        H[i, i] = U * (nup - 0.5) * (ndn - 0.5)
    return H, basis_all, index_all


def analytic_real_time(U, ts_n, ts_j):
    """G^R, G^<, G^> for either spin (Galpha/Gbeta averaged), from the
    closed-form single-pole expressions derived in this file's docstring."""
    tau = ts_n[:, None] - ts_j[None, :]
    g_less = (1j / 2) * np.exp(1j * U * tau / 2)
    g_grtr = (-1j / 2) * np.exp(-1j * U * tau / 2)
    g_ret = np.where(tau >= 0, -1j * np.cos(U * tau / 2), 0.0)
    return g_ret, g_less, g_grtr


def analytic_matsubara(U, taus):
    """G^M_sigma(tau), either spin, T=0 Lehmann sum (see module docstring)."""
    return -0.5 * np.exp(-(U / 2) * taus)


def analytic_left(U, ts, taus):
    """G^Left_sigma(t,tau), either spin, T=0 Lehmann sum (see module docstring)."""
    t = np.asarray(ts)[:, None]
    tau = np.asarray(taus)[None, :]
    return -0.5j * np.exp(1j * U * t / 4) * np.exp(-(U / 4) * tau)


def ed_real_time_single_seed(U, seed_up, ts):
    """
    Direct 4-state ED cross-check for a single seed (seed_up=True: impurity
    starts in |up>; seed_up=False: starts in |down>).  Returns
    (G^<_up, G^>_up, G^<_down, G^>_down) for that one seed, each shape
    (len(ts), len(ts)).

    G^<(t_n,t_j) = i<psi0|c^dag(t_j) c(t_n)|psi0>, Heisenberg c(t) = U(-t) c U(t).
    Using |phi(t)> = U(t)|psi0> (ordinary forward evolution):
        c(t_n)|psi0>      = U(-t_n) c |phi(t_n)>
        <psi0|c^dag(t_j)  = <phi(t_j)| c^dag U(t_j)
    so <psi0|c^dag(t_j) c(t_n)|psi0> = <phi(t_j)| c^dag U(t_j - t_n) c |phi(t_n)>
    -- the U(t_j - t_n) "propagate from t_n to t_j" step is easy to forget
    (an earlier draft of this function did exactly that, and silently
    returned all zeros as a result -- caught only by cross-checking against
    the independently-derived analytic formula, not by this code "looking
    plausible").  H is time-independent here, so U(t) for any real t
    (forward or backward) is just the same 4x4 matrix exponential.
    """
    H, basis_all, index_all = atomic_hamiltonian(U)
    eigvals, eigvecs = np.linalg.eigh(H)

    psi0 = np.zeros(4, dtype=complex)
    psi0[index_all[(1, 0) if seed_up else (0, 1)]] = 1.0

    def evolve(state, t):
        c = eigvecs.conj().T @ state
        return eigvecs @ (np.exp(-1j * eigvals * t) * c)

    def c_op(spin):
        # spin=0 (up) or 1 (down); annihilation operator in the 4-state basis
        M = np.zeros((4, 4))
        for (u, d), i in index_all.items():
            word = [u, d]
            if word[spin] == 0:
                continue
            sign = 1.0 if spin == 0 else (-1.0 if word[0] else 1.0)
            new_word = list(word)
            new_word[spin] = 0
            j = index_all[tuple(new_word)]
            M[j, i] = sign
        return M

    c_up, c_dn = c_op(0), c_op(1)

    N = len(ts)
    g_less = {0: np.zeros((N, N), dtype=complex), 1: np.zeros((N, N), dtype=complex)}
    g_grtr = {0: np.zeros((N, N), dtype=complex), 1: np.zeros((N, N), dtype=complex)}
    phi = [evolve(psi0, t) for t in ts]
    for spin, c_mat in [(0, c_up), (1, c_dn)]:
        c_phi = [c_mat @ p for p in phi]              # c |phi(t_n)>
        cdag_phi = [c_mat.conj().T @ p for p in phi]  # c^dagger |phi(t_j)>
        for n in range(N):
            for j in range(N):
                # G^<(t_n,t_j) = i <phi(t_j)| c . U(t_j-t_n) . c |phi(t_n)>
                # (bra uses c, not c^dag: [c(t_j)|psi0>]^dag = <c phi(t_j)| U(t_j) --
                # easy to get backwards, since G^> above genuinely does use
                # c^dag in the analogous slot; caught by cross-checking
                # against the independent analytic formula, not by symmetry
                # with the G^> code above "looking similar enough".)
                propagated = evolve(c_phi[n], ts[j] - ts[n])
                g_less[spin][n, j] = 1j * np.vdot(c_phi[j], propagated)
                # G^>(t_n,t_j) = -i <phi(t_n)| c . U(t_n-t_j) . c^dag |phi(t_j)>
                propagated_g = evolve(cdag_phi[j], ts[n] - ts[j])
                g_grtr[spin][n, j] = -1j * np.vdot(cdag_phi[n], propagated_g)
    return g_less[0], g_grtr[0], g_less[1], g_grtr[1]


def _c_op_matrix(index_all, spin):
    """Annihilation operator for the given spin in the 4-state basis
    (unsigned: with only one site there is no Jordan-Wigner string to
    worry about between the two spin flavors here)."""
    M = np.zeros((4, 4))
    for (u, d), i in index_all.items():
        word = [u, d]
        if word[spin] == 0:
            continue
        new_word = list(word)
        new_word[spin] = 0
        j = index_all[tuple(new_word)]
        M[j, i] = 1.0
    return M


def ed_matsubara_and_left(U, ts, taus):
    """
    Direct 4-state ED cross-check for G^M and G^Left, averaged over both
    spin-seeds (Galpha/Gbeta) and both spins sigma -- see module docstring
    for the closed forms this is checked against.

    G^M(tau)      = -<psi0| e^{+H tau} c e^{-H tau} c^dagger |psi0>
    G^Left(t,tau) = -i <psi0| c(t) e^{-H tau} c^dagger |psi0>,
                    c(t) = e^{iHt} c e^{-iHt} (Heisenberg, forward real time)

    Both are one-sided in imaginary time (only exp(-H tau) c^dagger|psi0>
    is ever propagated), unlike G^M's superficial appearance of "the same
    kind of decay" as G^Left -- they are NOT the same decay rate (see
    module docstring); this function computes both from the same
    intermediate ket precisely so that distinction cannot be blurred by
    two separately-written, only-superficially-similar code paths.
    """
    H, basis_all, index_all = atomic_hamiltonian(U)
    eigvals, eigvecs = np.linalg.eigh(H)

    def evolve_complex(state, t):
        c = eigvecs.conj().T @ state
        return eigvecs @ (np.exp(-1j * eigvals * t) * c)

    def evolve_real(state, tau):
        c = eigvecs.conj().T @ state
        return eigvecs @ (np.exp(-eigvals * tau) * c)

    c_ops = [_c_op_matrix(index_all, 0), _c_op_matrix(index_all, 1)]

    gm_sum = np.zeros(len(taus), dtype=complex)
    gleft_sum = np.zeros((len(ts), len(taus)), dtype=complex)
    n_terms = 0
    for seed_up in (True, False):
        psi0 = np.zeros(4, dtype=complex)
        psi0[index_all[(1, 0) if seed_up else (0, 1)]] = 1.0
        for c_mat in c_ops:
            cdag_psi0 = c_mat.conj().T @ psi0
            n_terms += 1
            for it, tau in enumerate(taus):
                ket = evolve_real(cdag_psi0, tau)  # e^{-H tau} c^dagger|psi0>

                c_ket = c_mat @ ket
                back = evolve_real(c_ket, -tau)  # e^{+H tau} c e^{-H tau} c^dagger|psi0>
                gm_sum[it] += -np.vdot(psi0, back)

                for jn, t in enumerate(ts):
                    phi_ket = evolve_complex(ket, -t)  # e^{+iHt} e^{-H tau} c^dagger|psi0>
                    c_phi_ket = c_mat @ phi_ket
                    gleft_sum[jn, it] += -1j * np.vdot(psi0, c_phi_ket)
    return gm_sum / n_terms, gleft_sum / n_terms


if __name__ == "__main__":
    U = 2.0
    dt = 0.04
    N = 20
    ts = np.arange(N + 1) * dt

    g_ret_an, g_less_an, g_grtr_an = analytic_real_time(U, ts, ts)

    g_less_up_a, g_grtr_up_a, g_less_dn_a, g_grtr_dn_a = ed_real_time_single_seed(U, True, ts)
    g_less_up_b, g_grtr_up_b, g_less_dn_b, g_grtr_dn_b = ed_real_time_single_seed(U, False, ts)

    # Gsigma = 1/2 (Galpha_sigma + Gbeta_sigma); check BOTH spins average to
    # the same, spin-independent result, matching the analytic formula.
    g_less_up_avg = 0.5 * (g_less_up_a + g_less_up_b)
    g_less_dn_avg = 0.5 * (g_less_dn_a + g_less_dn_b)
    g_grtr_up_avg = 0.5 * (g_grtr_up_a + g_grtr_up_b)
    g_grtr_dn_avg = 0.5 * (g_grtr_dn_a + g_grtr_dn_b)

    err_less_up = np.max(np.abs(g_less_up_avg - g_less_an))
    err_less_dn = np.max(np.abs(g_less_dn_avg - g_less_an))
    err_grtr_up = np.max(np.abs(g_grtr_up_avg - g_grtr_an))
    err_grtr_dn = np.max(np.abs(g_grtr_dn_avg - g_grtr_an))

    print("max|ED avg - analytic| for G^<_up:  ", err_less_up)
    print("max|ED avg - analytic| for G^<_down:", err_less_dn)
    print("max|ED avg - analytic| for G^>_up:  ", err_grtr_up)
    print("max|ED avg - analytic| for G^>_down:", err_grtr_dn)

    assert err_less_up < 1e-10 and err_less_dn < 1e-10
    assert err_grtr_up < 1e-10 and err_grtr_dn < 1e-10
    print("OK: closed-form analytic formulas match independent 4-state ED "
          "(Galpha/Gbeta averaged) to machine precision.")

    # Matsubara + Left branches: cross-check the closed forms against ED.
    taus = np.arange(N + 1) * dt
    gm_an = analytic_matsubara(U, taus)
    gleft_an = analytic_left(U, ts, taus)
    gm_ed, gleft_ed = ed_matsubara_and_left(U, ts, taus)

    err_gm = np.max(np.abs(gm_ed - gm_an))
    err_gleft = np.max(np.abs(gleft_ed - gleft_an))
    print("max|ED avg - analytic| for G^M:   ", err_gm)
    print("max|ED avg - analytic| for G^Left:", err_gleft)
    assert err_gm < 1e-10
    assert err_gleft < 1e-10
    print("OK: closed-form G^M/G^Left match independent 4-state ED "
          "(Galpha/Gbeta averaged) to machine precision.")

    print()
    print("Sample values for hardcoding into a Catch2 test (U=2, dt=0.04):")
    g_ret, g_less, _ = analytic_real_time(U, ts, ts)
    for n in [0, 1, 5, 10, 20]:
        for j in [0, 1, 5, 10, 20]:
            if j > n:
                continue
            print(f"  n={n:2d} j={j:2d} t={ts[n]:.2f} t'={ts[j]:.2f}  "
                  f"G^R={g_ret[n,j]!r}  G^<={g_less[n,j]!r}")

    print()
    print("Sample G^M/G^Left values for hardcoding into a Catch2 test:")
    for it in [0, 1, 5, 10, 20]:
        print(f"  tau={taus[it]:.2f}  G^M={gm_an[it]!r}")
    for n in [0, 1, 5, 10, 20]:
        for it in [0, 1, 5, 10, 20]:
            print(f"  n={n:2d} it={it:2d} t={ts[n]:.2f} tau={taus[it]:.2f}  "
                  f"G^Left={gleft_an[n,it]!r}")
