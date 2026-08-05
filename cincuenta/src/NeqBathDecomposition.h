#ifndef NEQ_BATH_DECOMPOSITION_H
#define NEQ_BATH_DECOMPOSITION_H
#include "KadanoffBaym.h"
#include <PsimagLite/Svd.h>
#include <PsimagLite/Vector.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>

// Bath decomposition for neq-DMFT following Gramsch, Balzer, Eckstein, Kollar
// PRB 88, 235106 (2013).
//
// Splits the hybridization Δ = Δ⁻ + Δ⁺:
//   Δ⁻ : first bath — memory of the initial equilibrium state,
//         computed analytically from the free propagators of the equilibrium {V_α, ε_α}.
//   Δ⁺ : second bath — neq dynamics accumulated for t>0,
//         represented via a rank-L low-rank Cholesky of i·Δ⁺_<.
//
// At each time step n (in order 0,1,2,...) call update(n, delta) after
// delta has been filled up to row n.  The resulting second-bath hoppings
// are available via Vplus(n, p).
//
// bathParams layout: [V_0 .. V_{N-1},  ε_0 .. ε_{N-1}]
namespace Dmft {

template <typename ComplexOrRealType> class NeqBathDecomposition {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorRealType    = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using MatrixComplexType = PsimagLite::Matrix<ComplexType>;
	using KBType            = KadanoffBaym<ComplexOrRealType>;

	NeqBathDecomposition(SizeType              rank,
	                     RealType              beta,
	                     RealType              mu,
	                     const VectorRealType& bathParams,
	                     SizeType              nT,
	                     SizeType              nTau,
	                     RealType              dt,
	                     RealType              dtau)
	    : rank_(rank)
	    , beta_(beta)
	    , mu_(mu)
	    , dt_(dt)
	    , dtau_(dtau)
	    , maxDiagSeen_(0)
	    , nT_(nT)
	    , nTau_(nTau)
	    , V_(nT + 1, std::max(rank, SizeType(1)), ComplexType(0))
	{
		const SizeType nBath = bathParams.size() / 2;
		hoppings_.resize(nBath);
		bathEps_.resize(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings_[i] = bathParams[i];
			bathEps_[i]  = bathParams[nBath + i];
		}
	}

	// Incremental Cholesky update for time step n.
	// delta.lesser(n, j) for j <= n must be filled before calling.
	//
	// n=0 is skipped: at t=0 the system is in equilibrium so Δ⁺<(0,0)=0.
	// Skipping ensures V_(0,p)=0 exactly. The standard Cholesky phase then
	// seeds column p from pivot row p+1 (n=1..L), giving true rank-L behavior.
	// Without the skip col 0 of V stays zero (L_eff=1) because V_(0,0)≈0 from
	// the bath-fit residual causes V_(1,0)=Δ⁺<(1,0)/V_(0,0) to blow up.
	void update(int n, const KBType& delta)
	{
		if (rank_ == 0)
			return;
		if (n == 0)
			return;

		const int L = static_cast<int>(rank_);

		maxDiagSeen_
		    = std::max(maxDiagSeen_, std::abs(-std::real(iDeltaPlusLesser(n, n, delta))));

		if (n <= L) {
			// Standard Cholesky: n=1 seeds col 0, n=2 seeds col 1, ..., n=L seeds col
			// L-1. Pivot for column p lives at row p+1 (row 0 is degenerate).
			const int col = n - 1;
			for (int p = 0; p < col; ++p)
				V_(n, p) = offDiagElement(n, p, delta);

			RealType d = -std::real(iDeltaPlusLesser(n, n, delta));
			for (int k = 0; k < col; ++k)
				d -= std::norm(V_(n, k));
			// Floor instead of clamping to exactly 0 -- see
			// gbek_cholesky.py::cholesky_causal's docstring (2026-07-09) for
			// the full story: the last seeding column's residual `d` is often
			// the tiny difference of O(1)-scale terms (a genuine property of
			// Lambda's available precision, verified via an independent
			// 50-digit mpmath replay, not a linear-algebra artifact), so its
			// sign is essentially undetermined noise. Clamping to EXACTLY 0
			// there is catastrophic: it makes that column's entire row/column
			// in every future Gram matrix identically zero, which then
			// structurally forces every later optimal-update solve to keep
			// that component at exactly zero forever (confirmed: this is
			// exactly what collapsed column 4 to zero for the whole run in a
			// real L=5 atomic-limit reconstruction). A small positive floor
			// keeps the column numerically alive so later rows -- seeing
			// genuine, non-noise data -- remain free to grow it into
			// something meaningful.
			const RealType floor
			    = RealType(1e-10) * std::max(maxDiagSeen_, RealType(1e-300));
			V_(n, col) = std::sqrt(std::max(d, floor));
		} else {
			// Optimal update: minimize ||Q q - a||² where Q_{kp} = conj(V_(k,p)),
			// k = 1..L (reference rows 1..L; row 0 is degenerate and excluded).
			choleskyOptimalUpdate(n, delta);
		}
	}

	// V^+(t_n, p): second-bath hopping for orbital p = 0..rank_-1
	ComplexType Vplus(int n, int p) const
	{
		assert(p >= 0 && static_cast<SizeType>(p) < rank_);
		return V_(n, p);
	}

	// V^-(t_n, α): first-bath hopping V_α exp(−i ε_α t_n)
	ComplexType Vminus(int n, int alpha) const
	{
		const ComplexType I(0, 1);
		return static_cast<ComplexType>(hoppings_[alpha])
		    * std::exp(-I * static_cast<ComplexType>(bathEps_[alpha] * n * dt_));
	}

	// Δ⁻_<(t_n, t_j): first-bath lesser hybridization
	ComplexType deltaMinusLesser(int n, int j) const
	{
		const ComplexType I(0, 1);
		ComplexType       result(0);
		for (SizeType a = 0; a < hoppings_.size(); ++a) {
			const RealType V2    = hoppings_[a] * hoppings_[a];
			const RealType eps   = bathEps_[a];
			const RealType fermi = fermiFunc(eps);
			result += V2 * I * fermi
			    * std::exp(-I * static_cast<ComplexType>(eps * (n - j) * dt_));
		}
		return result;
	}

	// Δ⁻_R(t_n, t_j): first-bath retarded hybridization  (n >= j)
	ComplexType deltaMinusRetarded(int n, int j) const
	{
		assert(n >= j);
		const ComplexType I(0, 1);
		ComplexType       result(0);
		for (SizeType a = 0; a < hoppings_.size(); ++a) {
			const RealType V2  = hoppings_[a] * hoppings_[a];
			const RealType eps = bathEps_[a];
			result += V2 * (-I)
			    * std::exp(-I * static_cast<ComplexType>(eps * (n - j) * dt_));
		}
		return result;
	}

	// Δ⁻_⌐(t_n, τ_j): first-bath left-mixing hybridization
	// Uses g_α^M(τ) = −f(ε_α) exp(−ε_α τ) for the free bath-site Matsubara GF
	ComplexType deltaMinusLeftMixing(int n, int j) const
	{
		const ComplexType I(0, 1);
		ComplexType       result(0);
		for (SizeType a = 0; a < hoppings_.size(); ++a) {
			const RealType V2    = hoppings_[a] * hoppings_[a];
			const RealType eps   = bathEps_[a];
			const RealType fermi = fermiFunc(eps);
			const RealType tau   = j * dtau_;
			result += V2 * std::exp(-I * static_cast<ComplexType>(eps * n * dt_))
			    * (-fermi) * static_cast<ComplexType>(std::exp(-eps * tau));
		}
		return result;
	}

	SizeType rank() const { return rank_; }
	SizeType nFirst() const { return hoppings_.size(); }

	// Reconstruct -iΔ⁺_<(t_n, t_j) = Σ_p V(n,p) conj(V(j,p)) and write to file.
	// File columns: n  j  Re(-iΔ⁺_<)  Im(-iΔ⁺_<)
	void dumpPlusBath(const std::string& filename) const
	{
		std::ofstream f(filename);
		for (SizeType n = 0; n <= nT_; ++n) {
			for (SizeType j = 0; j <= nT_; ++j) {
				ComplexType val(0);
				for (SizeType p = 0; p < rank_; ++p)
					val += V_(n, p) * std::conj(V_(j, p));
				f << n * dt_ << " " << j * dt_ << " " << val.real() << " "
				  << val.imag() << "\n";
			}
		}
	}

	// Dump the raw Cholesky factor V_(n,p) itself (not the reconstructed
	// product), for row-by-row comparison against an independent offline
	// trace (see cincuenta/TestSuite/gbek_reference -- the online-vs-offline
	// discrepancy diagnosis). File columns: n  p  Re(V)  Im(V)
	void dumpV(const std::string& filename) const
	{
		std::ofstream f(filename);
		// Full double precision: this comparison is used to diagnose
		// near-singular seeding rows where the default 6-sig-fig ostream
		// precision would itself be the dominant source of error.
		f.precision(17);
		for (SizeType n = 0; n <= nT_; ++n) {
			for (SizeType p = 0; p < rank_; ++p) {
				f << n << " " << p << " " << V_(n, p).real() << " "
				  << V_(n, p).imag() << "\n";
			}
		}
	}

private:

	// Fermi function f(ε) = 1 / (1 + exp(β(ε − μ)))
	RealType fermiFunc(RealType eps) const
	{
		const RealType x = beta_ * (eps - mu_);
		if (x > 500)
			return 0;
		if (x < -500)
			return 1;
		return 1 / (1 + std::exp(x));
	}

	// i Δ⁺_<(t_n, t_j): helper with antisymmetry for j > n
	ComplexType iDeltaPlusLesser(int n, int j, const KBType& delta) const
	{
		const ComplexType I(0, 1);
		ComplexType       dless;
		if (j <= n)
			dless = static_cast<ComplexType>(delta.lesser(n, j));
		else
			dless = -std::conj(static_cast<ComplexType>(delta.lesser(j, n)));
		return I * (dless - deltaMinusLesser(n, j));
	}

	// Standard Cholesky off-diagonal element V_(n, p) in the standard phase.
	// Pivot for column p is at row p+1 (row 0 is degenerate and skipped).
	ComplexType offDiagElement(int n, int p, const KBType& delta) const
	{
		ComplexType val = -iDeltaPlusLesser(n, p + 1, delta);
		for (int k = 0; k < p; ++k)
			val -= V_(n, k) * std::conj(V_(p + 1, k));
		const ComplexType pivot = V_(p + 1, p);
		if (std::abs(pivot) < 1e-14)
			return ComplexType(0);
		return val / std::conj(pivot);
	}

	// Optimal update for n > rank_: solve the JOINT minimization from GBEK
	// Eq. 63,
	//   F(q) = 2||Q_s q - a||^2 + |q^H q - d|^2,
	// where Q_s is the s x L matrix of ALL previously-determined rows 1..s
	// (s = n-1), growing by one row every step, and d = -iΔ⁺_<(t_n,t_n) is
	// the target diagonal value.
	//
	// Three bugs, fixed in sequence (see
	// cincuenta/TestSuite/gbek_reference/gbek_cholesky.py for the
	// independent Python implementation all three were caught against):
	//
	//  1. Fixed reference window: an earlier version of this function
	//     hardcoded the reference set to the first L seeding rows forever,
	//     instead of growing it to all s = n-1 rows. This discards almost
	//     all previously-built structure once n grows large relative to L,
	//     causing far-too-fast (near-total) collapse of the reconstructed
	//     hybridization.
	//  2. Missing diagonal constraint: after fixing (1), this function
	//     still solved only the linear part of Eq. 63 (Q^H Q q = Q^H a),
	//     ignoring the q^H q ~ d term -- NOT separable from the linear fit,
	//     since F(q) couples them. This systematically undershoots the
	//     diagonal by a smaller but still real amount (verified: ~5-15% on
	//     realistic full-rank targets, growing worse with n).
	//  3. Conjugation direction (found 2026-07-10, then REVERTED here on
	//     2026-07-10 after further investigation -- see below). The
	//     factorization this class targets is
	//     Lambda(n,k) = Sum_p V(n,p) * conj(V(k,p)), so solving a new row
	//     q=V(n,:) against established rows Q_raw(k,p)=V(k,p) means
	//     a = conj(Q_raw) @ q, i.e. the correct design matrix in the
	//     standard complex least-squares sense is conj(Q_raw), giving
	//     QtQ = Q_raw^T @ conj(Q_raw) = Sum_k V(k,p)*conj(V(k,pp)) and
	//     Qta = Q_raw^T @ a = Sum_k V(k,p)*a[k-1] (NO conjugate on V here).
	//     That is exactly what this function computed originally. A
	//     2026-07-08 session (commit 5a058f4d) found this function's
	//     output disagreeing with gbek_cholesky.py's Python
	//     `_solve_optimal_update` on a genuinely complex-phased target and
	//     "fixed" C++ to match Python -- but Python had the mirror-image
	//     bug at the time (its call sites passed the raw, unconjugated
	//     Q_raw into a solver that internally computes the standard
	//     Q^H Q / Q^H a form, which needs Q_raw pre-conjugated to be
	//     correct). Cross-validating two implementations that share the
	//     same bug confirms nothing (see gbek_cholesky.py's own module
	//     docstring for the general version of this lesson) -- in this
	//     case it actively made things worse, since the ORIGINAL C++ here
	//     was right and got swapped to match the wrong one. Re-derived
	//     2026-07-10 while chasing a ~50%-then-confirmed amplitude
	//     discrepancy against GBEK Fig. 8's double occupation, verified via
	//     (a) an exactly-solvable atomic-limit rank-1 target reconstructing
	//     to machine precision only with THIS convention, (b) a
	//     lambda-scaling test converging to the analytic O(v^2) benchmark
	//     only with THIS convention, and (c) the actual paper comparison
	//     landing within a few percent of Fig. 8's quoted peak only with
	//     THIS convention. The 2026-07-08 Catch2 regression test's
	//     hardcoded reference values were generated from the
	//     then-also-buggy Python and have been regenerated to match.
	//
	// a[k] = -iΔ⁺_<(t_n, t_{k+1}) for k = 0..s-1 (row k+1), one entry per
	// row of Q_s.
	void choleskyOptimalUpdate(int n, const KBType& delta)
	{
		const int L = static_cast<int>(rank_);
		const int s = n - 1; // number of previously-determined rows (1..n-1)

		VectorComplexType a(s);
		for (int k = 0; k < s; ++k)
			a[k] = -iDeltaPlusLesser(n, k + 1, delta);
		// Gram matrix QtQ[p][p'] = Sum_{k=1}^{s} V_(k,p) * conj(V_(k,p')).
		// See the bug-3 discussion above -- this is NOT literally "Q^H Q"
		// with Q_{kp}=V_(k,p); it's Q_raw^T @ conj(Q_raw), which is what
		// the underlying Hermitian factorization actually requires.
		MatrixComplexType QtQ(L, L, ComplexType(0));
		for (int p = 0; p < L; ++p) {
			for (int pp = 0; pp < L; ++pp) {
				ComplexType sum(0);
				for (int k = 1; k <= s; ++k)
					sum += V_(k, p) * std::conj(V_(k, pp));
				QtQ(p, pp) = sum;
			}
		}

		// Rhs Qta[p] = Sum_{k=1}^{s} V_(k,p) * a[k-1]  (i.e. Q_raw^T @ a,
		// NO conjugate on V -- see bug-3 discussion above).
		VectorComplexType Qta(L, ComplexType(0));
		for (int p = 0; p < L; ++p)
			for (int k = 1; k <= s; ++k)
				Qta[p] += V_(k, p) * a[k - 1];

		const RealType    d = -std::real(iDeltaPlusLesser(n, n, delta));
		VectorComplexType q = solveOptimalUpdateJoint(QtQ, Qta, d, L);
		for (int p = 0; p < L; ++p)
			V_(n, p) = q[p];
	}

	// Solve GBEK Eq. 63's joint objective F(q) = 2||Qq-a||^2 + |q^H q-d|^2
	// exactly. Setting the Wirtinger gradient dF/dq* to zero gives (the
	// factors of 2 cancel completely -- easy to get wrong; a draft of this
	// fix kept a stray factor of 2 and got answers close to, but measurably
	// different from, the true minimum found by numerically minimizing F(q)
	// directly with BFGS in Python):
	//
	//   [QtQ + (mu - d) I] q(mu) = Qta,   mu = q(mu)^H q(mu)
	//
	// This is ONE nonlinear SCALAR equation in mu (not a joint multivariate
	// optimization): for fixed mu, q(mu) solves a small L x L linear
	// system; g(mu) = ||q(mu)||^2 - mu = 0 is then a 1D root-find, solved
	// here by bisection after bracketing -- no nonlinear-optimization
	// library needed.
	//
	// At mu = d, this reduces exactly to the old (bug 2) linear-only
	// solution: [QtQ] q = Qta. ||q(mu)||^2 is monotonically decreasing in mu
	// (standard secular-equation property: A(mu) = QtQ + (mu-d)I only gains
	// positive-definiteness as mu grows, so its inverse only shrinks), hence
	// g(mu) = ||q(mu)||^2 - mu is strictly decreasing too. Since g(0) >= 0
	// generically and g(mu) -> -infinity as mu -> infinity, exactly one root
	// exists somewhere in (0, infinity) -- but NOT always in [0, d]. When
	// g(d) >= 0 (the diagonal-unconstrained fit q(d) already has norm >= d),
	// monotonicity forces the true root above d, not at d: bracket upward
	// instead of returning q(d), which in that regime is not even a
	// stationary point of Eq. 63. (Found 2026-07-09: the previous version of
	// this function only ever searched [0, d] and silently accepted q(d) --
	// including its norm overshooting d -- whenever the linear-only fit's
	// norm already met or exceeded the diagonal target. This under-corrects
	// the diagonal on targets whose off-diagonal magnitude does not decay,
	// e.g. GBEK's true atomic-limit self-consistency target, producing a
	// spurious diagonal overshoot distinct from the paper's own
	// documented decay-only low-rank artifact.)
	//
	// Bug 4 (found 2026-07-10, ported from gbek_cholesky.py's
	// _solve_optimal_update -- see that function's docstring for the full
	// derivation and the brute-force-optimization verification): qOf(mu)
	// used to solve [QtQ+(mu-d)I]q=Qta via a FRESH SVD-based solveLinear()
	// call at every mu tried during bisection. Fine when QtQ is
	// well-conditioned, but at low rank it's common for a just-established
	// column to carry almost no signal yet (observed on real GBEK data:
	// QtQ eigenvalues [1.1e-8, 7.4], condition number ~7e8). As mu varies,
	// the shifted near-zero eigenvalue crosses solveLinear's rcond
	// threshold inconsistently, making the numerically computed g(mu)
	// behave non-monotonically instead of the smooth, strictly-decreasing
	// function bisection assumes -- confirmed by comparing against a
	// brute-force numerical minimization of F(q): the old solver's answer
	// had ||q||^2 substantially undershooting d (bug 2's symptom recurring
	// via a different mechanism) and a measurably higher F(q) than the true
	// minimum. Fixed by diagonalizing QtQ ONCE (Hermitian, real eigenvalues)
	// and evaluating g(mu) via the exact secular-equation form
	// g(mu) = sum_i |b_i|^2/(lambda_i+mu-d)^2 - mu, b = W^H Qta in QtQ's
	// eigenbasis -- an explicit algebraic function of mu with no linear
	// solve (hence no rcond-driven rank flips) per trial mu. Also handles
	// the classic trust-region-subproblem "hard case" explicitly: when the
	// smallest eigenvalue AND its b-component are both negligible, that
	// direction is a near-free parameter for the linear fit but still
	// counts fully toward ||q||^2=d, and the generic secular equation has a
	// genuine pole there that bisection cannot resolve reliably (the true
	// root can sit within ~1e-12 of it). Handled by building q from the
	// well-conditioned directions alone, then setting the near-null
	// direction's magnitude to whatever hits ||q||^2=d exactly, using b's
	// own (tiny but nonzero) phase as the natural choice.
	static VectorComplexType solveOptimalUpdateJoint(const MatrixComplexType& QtQ,
	                                                 const VectorComplexType& Qta,
	                                                 RealType                 d,
	                                                 int                      L)
	{
		MatrixComplexType W(QtQ); // zheev overwrites W with eigenvectors (columns)
		VectorRealType    lam; // ascending eigenvalues
		PsimagLite::diag(W, lam, 'V');

		VectorComplexType b(L, ComplexType(0));
		for (int i = 0; i < L; ++i)
			for (int k = 0; k < L; ++k)
				b[i] += std::conj(W(k, i)) * Qta[k];

		auto qOf = [&](RealType mu)
		{
			VectorComplexType q(L, ComplexType(0));
			for (int i = 0; i < L; ++i) {
				const ComplexType c = b[i] / ComplexType(lam[i] + mu - d, 0);
				for (int p = 0; p < L; ++p)
					q[p] += W(p, i) * c;
			}
			return q;
		};
		auto gOf = [&](RealType mu)
		{
			RealType sum = 0;
			for (int i = 0; i < L; ++i) {
				const RealType denom = lam[i] + mu - d;
				sum += std::norm(b[i]) / (denom * denom);
			}
			return sum - mu;
		};

		// Hard-case check (see docstring above).
		{
			const RealType    lamMax = lam[L - 1];
			std::vector<bool> hard(L, false);
			RealType          hardNorm2 = 0;
			int               nHard     = 0;
			for (int i = 0; i < L; ++i) {
				if (lam[i] < 1e-6 * lamMax) {
					hard[i] = true;
					hardNorm2 += std::norm(b[i]);
					++nHard;
				}
			}
			if (nHard > 0 && hardNorm2 < 1e-6 * std::max(d, RealType(1.0))) {
				RealType softNorm2 = 0;
				for (int i = 0; i < L; ++i)
					if (!hard[i])
						softNorm2 += std::norm(b[i]) / (lam[i] * lam[i]);
				const RealType remainder = d - softNorm2;
				if (remainder > 0) {
					const RealType    tau = std::sqrt(remainder / nHard);
					VectorComplexType q(L, ComplexType(0));
					for (int i = 0; i < L; ++i) {
						ComplexType c;
						if (hard[i]) {
							const RealType    absB  = std::abs(b[i]);
							const ComplexType phase = (absB > 0)
							    ? (b[i] / absB)
							    : ComplexType(1, 0);
							c = ComplexType(tau, 0) * phase;
						} else {
							c = b[i] / ComplexType(lam[i], 0);
						}
						for (int p = 0; p < L; ++p)
							q[p] += W(p, i) * c;
					}
					RealType n2Hard    = 0;
					bool     allFinite = true;
					for (int p = 0; p < L; ++p) {
						if (!std::isfinite(q[p].real())
						    || !std::isfinite(q[p].imag()))
							allFinite = false;
						n2Hard += std::norm(q[p]);
					}
					if (allFinite
					    && std::abs(n2Hard - d)
					        < 1e-6 * std::max(d, RealType(1.0)))
						return q;
				}
			}
		}

		RealType muLo;
		RealType muHi;
		if (gOf(d) >= 0) {
			// True root is above d: bracket upward.
			muLo           = d;
			muHi           = d;
			RealType step  = std::max(d, RealType(1.0)) * 0.5;
			bool     found = false;
			for (int iter = 0; iter < 60; ++iter) {
				muHi += step;
				if (gOf(muHi) <= 0) {
					found = true;
					break;
				}
				step *= 1.5;
			}
			if (!found)
				return qOf(d);
		} else {
			// True root is below d, as in the original implementation.
			muHi          = d;
			muLo          = d;
			RealType step = std::max(d, RealType(1.0)) * 0.5;
			for (int iter = 0; iter < 60; ++iter) {
				muLo = std::max(muLo - step, RealType(0.0));
				if (gOf(muLo) >= 0 || muLo <= 0.0)
					break;
				step *= 1.5;
			}
		}

		for (int iter = 0; iter < 100; ++iter) {
			const RealType muMid = 0.5 * (muLo + muHi);
			if (gOf(muMid) >= 0)
				muLo = muMid;
			else
				muHi = muMid;
			if (muHi - muLo < 1e-14 * std::max(muHi, RealType(1.0)))
				break;
		}
		return qOf(muLo);
	}

	// Least-squares solve of A q = b (A is L×L, Hermitian in practice since
	// it is QtQ + scalar*I, but solved generally) via SVD-based pseudo-
	// inverse, truncating singular values below rcond*max(s) -- matching
	// numpy.linalg.lstsq's default behavior.
	//
	// This replaced a plain Gaussian-elimination-with-partial-pivoting
	// solve (skip-if-|pivot|<1e-14) that was found, via independent
	// cross-check against cincuenta/TestSuite/gbek_reference/gbek_cholesky.py
	// on the true atomic-limit target (Delta^- === 0 exactly, an extreme,
	// genuinely rank-deficient seed), to silently return large-but-finite
	// garbage on NEAR-singular (not exactly singular) systems -- exactly
	// the "near-singular Gram matrix blowup" failure mode documented in
	// gbek_cholesky.py's module docstring, which was fixed on the Python
	// side (switched to np.linalg.lstsq) but never ported here. The two
	// implementations agreed on the near-atomic (W_i=0.1) target because
	// that target's early rows are not exactly singular, only degenerate;
	// the true atomic limit's Delta_<(0,0)=0 exactly is a harsher seed and
	// exposed the gap.
	static VectorComplexType
	solveLinear(const MatrixComplexType& A, const VectorComplexType& b, int L)
	{
		MatrixComplexType            u(A);
		VectorRealType               s;
		MatrixComplexType            vt;
		PsimagLite::Svd<ComplexType> svd;
		svd('S', u, s, vt);

		// rcond=1e-10 matches gbek_cholesky.py::_solve_optimal_update's
		// np.linalg.lstsq(A, Qta, rcond=1e-10) exactly -- deliberately far
		// looser than machine-epsilon truncation, since the near-singular
		// directions here are not noise but genuinely poorly-determined
		// combinations that must be dropped, not merely rounding error.
		const RealType sMax
		    = (s.empty()) ? RealType(0) : *std::max_element(s.begin(), s.end());
		const RealType rcond  = RealType(1e-10);
		const RealType thresh = rcond * sMax;

		// x = V * diag(1/s_i, truncated) * U^H * b
		VectorComplexType uhB(L, ComplexType(0));
		for (int i = 0; i < L; ++i)
			for (int k = 0; k < L; ++k)
				uhB[i] += std::conj(u(k, i)) * b[k];

		VectorComplexType x(L, ComplexType(0));
		for (int i = 0; i < L; ++i) {
			if (s[i] <= thresh)
				continue;
			const ComplexType coeff = uhB[i] / s[i];
			for (int c = 0; c < L; ++c)
				x[c] += std::conj(vt(i, c)) * coeff;
		}
		return x;
	}

	SizeType          rank_;
	RealType          beta_, mu_, dt_, dtau_;
	RealType          maxDiagSeen_; // running max |diag| for the seeding floor, see update()
	SizeType          nT_, nTau_;
	VectorRealType    hoppings_; // V_α from equilibrium bath
	VectorRealType    bathEps_; // ε_α from equilibrium bath
	MatrixComplexType V_; // second-bath hoppings V_[n][p], (nT+1) × rank_
};

} // namespace Dmft
#endif // NEQ_BATH_DECOMPOSITION_H
