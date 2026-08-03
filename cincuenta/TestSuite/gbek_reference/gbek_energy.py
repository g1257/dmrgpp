"""
Kinetic/interaction/total energy from the GBEK atomic-limit self-consistency
loop (see gbek_selfconsistency.py), needed for paper Fig. 9 ("test of energy
conservation") and Fig. 10's Eint columns.

Physics (unnumbered equations in the paper, Sec. VI / Fig. 9-10 discussion):
    <Ekin(t)>  = -i Sum_sigma [Lambda_sigma * G_sigma]^<(t,t)   (Keldysh
                 contour convolution, evaluated on the diagonal)
    <Eint(t)>  = U * (<d(t)> - 1/4)
    <Etot(t)>  = <Ekin(t)> + <Eint(t)>

Array conventions (must match gbek_selfconsistency.py / gbek_dynamics.py
exactly -- this is where three of the six prior Cholesky/GBEK bugs in this
codebase came from, all sign/conjugation mistakes in exactly this kind of
two-time bookkeeping):
    - Lambda_less/Lambda_great are hop(t)*hop(t')*(-i)*G^{<,>}(t,t') (see
      gbek_selfconsistency.py's Lambda_new construction), stored
      lower-triangular (array[n,j] = value at (t_n,t_j), only j<=n
      filled). This quantity is HERMITIAN: X(t,t') = conj(X(t',t)) -- the
      extra -i flips the parity of G^{<,>}'s own anti-Hermitian relation
      (see gbek_selfconsistency.py's dump_lesser() docstring).
    - G_less/G_great are G^{<,>}(t,t') themselves (compute_g_lesser's and
      compute_g_greater's direct return values, no extra prefactor), also
      stored lower-triangular. This quantity is ANTI-Hermitian:
      X(t,t') = -conj(X(t',t)) -- the standard relation for a single-particle
      lesser/greater Green's function.
    - G_less/G_great passed to compute_kinetic_energy should be the
      SPIN-AVERAGED G_avg = 0.5*(G_up+G_dn): Lambda does not depend on spin in
      this SU(2)-symmetric problem, so Sum_sigma Lambda*G_sigma =
      Lambda*(G_up+G_dn) = 2*Lambda*G_avg -- the factor of 2 is applied
      explicitly inside compute_kinetic_energy.
"""
import numpy as np


def _full_hermitian(X_lower):
    """
    Reconstruct the full (N+1,N+1) matrix of a Hermitian-type two-time
    quantity (X(t,t') = conj(X(t',t))) from its lower-triangular storage
    (X_lower[n,j] = X(t_n,t_j), j<=n). Matches the convention already used
    for Lambda in gbek_selfconsistency.py's dump_lesser()/_probe().
    """
    N = X_lower.shape[0] - 1
    full = np.array(X_lower, dtype=complex, copy=True)
    for n in range(N + 1):
        for j in range(n + 1, N + 1):
            full[n, j] = np.conj(X_lower[j, n])
    return full


def _full_antihermitian(X_lower):
    """
    Same as _full_hermitian but for anti-Hermitian two-time quantities
    (X(t,t') = -conj(X(t',t))) -- the standard relation for G^{<,>} itself,
    as opposed to Lambda (Hermitian, see _full_hermitian).
    """
    N = X_lower.shape[0] - 1
    full = np.array(X_lower, dtype=complex, copy=True)
    for n in range(N + 1):
        for j in range(n + 1, N + 1):
            full[n, j] = -np.conj(X_lower[j, n])
    return full


def _retarded_full(X_less_full, X_great_full):
    """X^R(t,t') = theta(t-t') * [X^>(t,t') - X^<(t,t')], full (N+1,N+1) grid."""
    N = X_less_full.shape[0] - 1
    theta = np.tril(np.ones((N + 1, N + 1)))
    return theta * (X_great_full - X_less_full)


def compute_kinetic_energy(Lambda_less, Lambda_great, G_less, G_great, dt):
    """
    <Ekin(t_n)> via the Langreth rule for [Lambda*G]^<(t,t):
        [A*B]^<(t,t) = int_0^t ds [A^R(t,s) B^<(s,t) + A^<(t,s) B^A(s,t)]

    Lambda_less/Lambda_great, G_less/G_great: (N+1,N+1) complex arrays in the
    exact storage conventions documented in this module's docstring.

    Returns a length-(N+1) complex array. Physical energies are real -- the
    imaginary part should vanish to numerical precision; callers should
    assert this explicitly rather than silently discarding it (e.g. via
    np.imag(Ekin_t)).
    """
    N = Lambda_less.shape[0] - 1
    Lam_less_full = _full_hermitian(Lambda_less)
    Lam_great_full = _full_hermitian(Lambda_great)
    G_less_full = _full_antihermitian(G_less)
    G_great_full = _full_antihermitian(G_great)

    Lam_R_full = _retarded_full(Lam_less_full, Lam_great_full)
    G_R_full = _retarded_full(G_less_full, G_great_full)
    G_A_full = np.conj(G_R_full.T)  # G^A(t,t') = [G^R(t',t)]^*

    Ekin = np.zeros(N + 1, dtype=complex)
    for n in range(N + 1):
        integrand = (Lam_R_full[n, :n + 1] * G_less_full[:n + 1, n]
                     + Lam_less_full[n, :n + 1] * G_A_full[:n + 1, n])
        Ekin[n] = -1j * 2.0 * np.trapezoid(integrand, dx=dt)
    return Ekin


def compute_interaction_energy(d_t, U):
    """<Eint(t)> = U*(<d(t)> - 1/4)."""
    return U * (np.asarray(d_t) - 0.25)


def compute_total_energy(Ekin_t, Eint_t):
    """<Etot(t)> = <Ekin(t)> + <Eint(t)>."""
    return np.asarray(Ekin_t) + np.asarray(Eint_t)
