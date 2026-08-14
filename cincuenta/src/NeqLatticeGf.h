#ifndef NEQ_LATTICE_GF_H
#define NEQ_LATTICE_GF_H

#include "KadanoffBaym.h"
#include "ParamsNeqDmftSolver.h"
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>
#include <cassert>
#include <cmath>
#include <complex>

namespace Dmft {

/*!
 * \brief NeqLatticeGf
 * Non-equilibrium Weiss field G_0(t,t') for the Bethe lattice.
 *
 * Self-consistency (Bethe lattice):
 *   Δ(t,t') = t*² G_imp(t,t')    where t* = D = W/2 (half-bandwidth)
 *
 * Dyson equation (Volterra integro-differential):
 *   [i d/dt - μ] G_0(t,t') = δ_C(t,t') + (Δ ⊛ G_0)(t,t')
 *
 * Usage per time step n:
 *   1. Call initialize(gimp) once after the equilibrium run (sets t=0 BCs).
 *   2. updateDelta(n, gimp) — copy t*² G_imp → Δ for the n-th row.
 *   3. advance(n) — solve the Volterra equation for G_0(n, j), j ≤ n.
 */
template <typename ComplexOrRealType> class NeqLatticeGf {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using VectorStringType  = PsimagLite::Vector<PsimagLite::String>::Type;
	using KBType            = KadanoffBaym<ComplexOrRealType>;
	using KBDerivType       = KBDerivative<ComplexOrRealType>;
	using ParamsNeqType     = ParamsNeqDmftSolver<ComplexOrRealType>;

	explicit NeqLatticeGf(const ParamsNeqType& params)
	    : params_(params)
	    , nTau_(params.eqParams.nMatsubaras)
	    , dtau_(params.eqParams.ficticiousBeta
	            / static_cast<RealType>(params.eqParams.nMatsubaras))
	    , tStar_(params.neqAtomicLimit ? RealType(0) : parseTstar(params.eqParams.latticeGf))
	    , tStarSq_(tStar_ * tStar_)
	    , tStarFinal_(params.bandwidthFinal > RealType(0)
	                      ? RealType(0.25) * params.bandwidthFinal
	                      : tStar_)
	    , g0_(params.nT,
	          params.eqParams.nMatsubaras,
	          params.dt,
	          params.eqParams.ficticiousBeta
	              / static_cast<RealType>(params.eqParams.nMatsubaras))
	    , delta_(params.nT,
	             params.eqParams.nMatsubaras,
	             params.dt,
	             params.eqParams.ficticiousBeta
	                 / static_cast<RealType>(params.eqParams.nMatsubaras))
	    , g0_der_(params.nT, params.eqParams.nMatsubaras)
	    , g0_der_new_(params.nT, params.eqParams.nMatsubaras)
	    , h_(params.nT + 1, ComplexType(params.eqParams.mu, 0))
	{ }

	/*!
	 * \brief initialize
	 * Set t=0 boundary conditions and Matsubara components from equilibrium G_imp.
	 *
	 * G_0^M(iω_k) = 1 / (iω_k + μ - t*² G_imp^M(iω_k))
	 * G_0^M(τ_j) via inverse Matsubara sum
	 * G_0^{Left}(0,j) = -i G_0^M(β - τ_j)
	 * G_0^R(0,0) = -i
	 * G_0^<(0,0) = G_0^{Left}(0,0) = -i G_0^M(β)
	 *
	 * \param[in] gimp Equilibrium impurity Green's function.
	 */
	void initialize(const KBType& gimp)
	{
		const int         Ntau = static_cast<int>(nTau_);
		const RealType    beta = params_.eqParams.ficticiousBeta;
		const RealType    mu   = params_.eqParams.mu;
		const ComplexType I(0, 1);

		// Δ^M and G_0^M in Matsubara frequency
		for (SizeType k = 0; k < nTau_; ++k)
			delta_.matsubara_w[k] = tStarSq_ * gimp.matsubara_w[k];

		for (SizeType k = 0; k < nTau_; ++k) {
			const RealType omk = matsubaraFreq(k, nTau_, beta);
			g0_.matsubara_w[k]
			    = ComplexType(1) / (I * omk + mu - delta_.matsubara_w[k]);
		}

		// G_0^M(τ_j) = (1/β) Σ_k G_0^M(iω_k) exp(-iω_k τ_j)
		for (SizeType j = 0; j <= nTau_; ++j) {
			const RealType tau = j * dtau_;
			ComplexType    gm  = 0;
			for (SizeType k = 0; k < nTau_; ++k) {
				const RealType omk   = matsubaraFreq(k, nTau_, beta);
				const RealType phase = -omk * tau;
				gm += g0_.matsubara_w[k]
				    * ComplexType(std::cos(phase), std::sin(phase));
			}
			g0_.matsubara_t[j] = gm / beta;
		}

		// Δ^M(τ_j) = t*² G_imp^M(τ_j)
		for (SizeType j = 0; j <= nTau_; ++j)
			delta_.matsubara_t[j] = tStarSq_ * gimp.matsubara_t[j];

		// t=0 imaginary-time slice
		// G_0^{Left}(0,j) = -i G_0^M(β - τ_j) = -i matsubara_t[nTau - j]
		for (int j = 0; j <= Ntau; ++j)
			g0_.left_mixing(0, j) = -I * g0_.matsubara_t[Ntau - j];

		// Δ^{Left}(0,j) = t*² G_imp^{Left}(0,j)
		for (int j = 0; j <= Ntau; ++j)
			delta_.left_mixing(0, j) = tStarSq_ * gimp.left_mixing(0, j);

		// t=0 retarded boundary condition
		g0_.retarded(0, 0)    = ComplexType(0, -1);
		delta_.retarded(0, 0) = tStarSq_ * gimp.retarded(0, 0);

		// t=0 lesser: G_0^<(0,0) = i * n  where n is the equilibrium occupancy.
		// The Matsubara sum at tau=beta suffers high-frequency tail truncation and
		// gives the wrong value; use the positive-frequency sum instead, which
		// converges correctly: n = 1/2 + (1/beta) Re[sum_{omega>0} G_0^M(i*omega)].
		{
			const RealType betaVal = params_.eqParams.ficticiousBeta;
			const SizeType halfN   = nTau_ / 2; // first positive-frequency index
			ComplexType    posSum  = 0;
			for (SizeType k = halfN; k < nTau_; ++k)
				posSum += g0_.matsubara_w[k];
			const RealType n_weiss = RealType(0.5) + PsimagLite::real(posSum) / betaVal;
			g0_.lesser(0, 0)       = ComplexType(0, 1) * n_weiss;
			g0_.left_mixing(0, 0)
			    = g0_.lesser(0, 0); // left_mixing(0, tau=0) = lesser(0,0)
		}
		delta_.lesser(0, 0) = tStarSq_ * gimp.lesser(0, 0);

		// Initial RK derivatives d/dt G_0(0, ·) needed for the n=1 predictor
		computeDerivativesAt0();
	}

	/*!
	 * \brief updateDelta
	 * Compute Delta(t_n, t_j) = t*(t_n) t*(t_j) G_imp for all retarded, lesser,
	 * and left-mixing components at row n.  t*(t) follows the ramp shape
	 * specified by params_.quenchShape / params_.quenchDuration.
	 *
	 * \param[in] n    Time row index.
	 * \param[in] gimp Impurity Green's function at time step n.
	 */
	void updateDelta(int n, const KBType& gimp)
	{
		const RealType tsn = tStarAt(n);
		for (int j = 0; j <= n; ++j) {
			const RealType tsj    = tStarAt(j);
			delta_.retarded(n, j) = tsn * tsj * gimp.retarded(n, j);
			delta_.lesser(n, j)   = tsn * tsj * gimp.lesser(n, j);
			delta_.lesser(j, n)   = tsj * tsn * gimp.lesser(j, n);
		}
		// Left-mixing: real-time hopping x imaginary-time (equilibrium) hopping.
		for (SizeType j = 0; j <= nTau_; ++j)
			delta_.left_mixing(n, j) = tsn * tStar_ * gimp.left_mixing(n, j);
	}

	/*!
	 * \brief advance
	 * Advance G_0 to time step n via volterra_intdiff.
	 * Precondition: updateDelta(n, gimp) called for n and all n' < n.
	 *
	 * \param[in] n Time row index to advance to.
	 */
	void advance(int n)
	{
		assert(n >= 1);
		g0_.volterra_intdiff(n, h_, delta_, g0_der_, g0_der_new_);
		g0_der_.update(static_cast<SizeType>(n), nTau_, g0_der_new_);
	}

	const KBType& g0() const { return g0_; }
	const KBType& delta() const { return delta_; }

private:

	// Evaluate t*(t_n = n*dt) according to the quench ramp shape.
	// "cosine" matches GBEK PRB 88, 235106 (recommended t_q=0.25).
	// "tanh"   sigmoid centered at t_q/2 with characteristic width t_q/6.
	// "step" (or quenchDuration=0): t*_i at n=0, t*_f for n>=1.
	RealType tStarAt(int n) const
	{
		const RealType tq = params_.quenchDuration;
		if (tq <= RealType(0) || params_.quenchShape == "step")
			return (n == 0) ? tStar_ : tStarFinal_;

		const RealType t = n * params_.dt;
		if (t >= tq)
			return tStarFinal_;

		RealType shape;
		if (params_.quenchShape == "tanh") {
			const RealType x = (t / tq - RealType(0.5)) * RealType(6);
			shape            = RealType(0.5) * (RealType(1) + std::tanh(x));
		} else { // "cosine" (default for any unrecognised value)
			shape = RealType(0.5) * (RealType(1) - std::cos(M_PI * t / tq));
		}
		return tStar_ + (tStarFinal_ - tStar_) * shape;
	}

	/*!
	 * \brief parseTstar
	 * Extract the Bethe lattice hopping t* = D/2 = W/4 from "energy,semicircular,W".
	 * D = W/2 is the half-bandwidth; t* = D/2 satisfies <epsilon^2> = t*^2
	 * for the semicircular DOS, which is the coefficient in Delta = t*^2 G_latt.
	 *
	 * \param[in] latticeGf Lattice GF descriptor string, e.g. "energy,semicircular,W".
	 * \return t* = W/4.
	 */
	static RealType parseTstar(const PsimagLite::String& latticeGf)
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, latticeGf, ",");
		if (tokens.size() < 3)
			err("NeqLatticeGf: LatticeGf must be 'energy,semicircular,W'; got: "
			    + latticeGf + "\n");
		return RealType(0.25) * PsimagLite::atof(tokens[2]); // t* = W/4
	}

	/*!
	 * \brief matsubaraFreq
	 * Fermionic Matsubara frequency ω_k = (2k - N + 1) π / β, k = 0..N-1.
	 *
	 * \param[in] k    Frequency index.
	 * \param[in] N    Total number of frequencies.
	 * \param[in] beta Inverse temperature.
	 * \return ω_k.
	 */
	static RealType matsubaraFreq(SizeType k, SizeType N, RealType beta)
	{
		return RealType(2 * static_cast<int>(k) - static_cast<int>(N) + 1) * M_PI / beta;
	}

	/*!
	 * \brief trapz
	 * Trapezoidal quadrature: half-weight at endpoints.  Returns 0 for i == j.
	 *
	 * \param[in] f Array of values to integrate.
	 * \param[in] i Start index (inclusive).
	 * \param[in] j End index (inclusive).
	 * \return Trapezoidal sum over f[i..j].
	 */
	static ComplexType trapz(const VectorComplexType& f, int i, int j)
	{
		if (i == j)
			return ComplexType(0);
		ComplexType s = RealType(0.5) * (f[i] + f[j]);
		for (int k = i + 1; k < j; ++k)
			s += f[k];
		return s;
	}

	/*!
	 * \brief computeDerivativesAt0
	 * Compute d/dt G_0(t=0, ·) by evaluating the ODE at t=0.
	 *
	 * At t=0 the real-time Volterra integral is empty, leaving only:
	 *   d/dt G_0^R(0,0)     = -i μ G_0^R(0,0)
	 *   d/dt G_0^{Left}(0,j) = +I dtau ∫_0^{τ_j} Δ^L(0,l) G_0^M(β+τ_l-τ_j) dτ
	 *                         - I dtau ∫_{τ_j}^β  Δ^L(0,l) G_0^M(τ_l-τ_j)   dτ
	 *                         - i μ G_0^{Left}(0,j)
	 *   d/dt G_0^<(0,0)     = -dtau ∫_0^β Δ^L(0,l) [G_0^{Left}(0,β-τ_l)]^* dτ
	 *                         - i μ G_0^<(0,0)
	 *
	 * Signs follow the volterra_intdiff convention (match factor -I*(-I) for lesser,
	 * +I and -I for the two left-mixing Matsubara halves).
	 */
	void computeDerivativesAt0()
	{
		const ComplexType I(0, 1);
		const int         Ntau = static_cast<int>(nTau_);
		const RealType    mu   = params_.eqParams.mu;
		VectorComplexType tmp(Ntau + 1);

		// Retarded diagonal
		g0_der_.retarded[0] = -I * mu * g0_.retarded(0, 0);

		// Left-mixing
		for (int j = 0; j <= Ntau; ++j) {
			for (int l = 0; l <= j; ++l)
				tmp[l] = delta_.left_mixing(0, l) * g0_.matsubara_t[Ntau + l - j];
			g0_der_.left_mixing[j] = I * dtau_ * trapz(tmp, 0, j);
			for (int l = j; l <= Ntau; ++l)
				tmp[l] = delta_.left_mixing(0, l) * g0_.matsubara_t[l - j];
			g0_der_.left_mixing[j] -= I * dtau_ * trapz(tmp, j, Ntau);
			g0_der_.left_mixing[j] -= I * mu * g0_.left_mixing(0, j);
		}

		// Lesser diagonal — coefficient matches -I*(-I) = -1 from volterra_intdiff
		for (int l = 0; l <= Ntau; ++l)
			tmp[l] = delta_.left_mixing(0, l)
			    * PsimagLite::conj(g0_.left_mixing(0, Ntau - l));
		g0_der_.lesser[0] = -dtau_ * trapz(tmp, 0, Ntau) - I * mu * g0_.lesser(0, 0);
	}

	const ParamsNeqType& params_;
	SizeType             nTau_;
	RealType             dtau_;
	RealType             tStar_;
	RealType             tStarSq_;
	RealType             tStarFinal_; ///< t*_f (post-quench); equals tStar_ if BandwidthFinal=0
	KBType               g0_;
	KBType               delta_;
	KBDerivType          g0_der_;
	KBDerivType          g0_der_new_;
	VectorComplexType    h_; ///< h[n] = μ (constant single-particle term)
};

} // namespace Dmft
#endif // NEQ_LATTICE_GF_H
