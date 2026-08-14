#ifndef IMPURITYSOLVER_NEQ_BASE_H
#define IMPURITYSOLVER_NEQ_BASE_H
#include "CincuentaInputCheck.h"
#include "KadanoffBaym.h"
#include "ModelParams.h"
#include <PsimagLite/InputNg.h>
#include <PsimagLite/Matsubaras.h>
#include <PsimagLite/Vector.h>

namespace Dmft {

/*!
 * \brief ImpuritySolverNeqBase
 * Abstract interface for non-equilibrium impurity solvers.
 * Separates the two-time G_imp computation from the equilibrium ImpuritySolverBase
 * hierarchy (which returns a 1D Matsubara vector) to avoid polluting that interface.
 */
template <typename ComplexOrRealType> class ImpuritySolverNeqBase {

public:

	using RealType        = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType     = std::complex<RealType>;
	using VectorRealType  = typename PsimagLite::Vector<RealType>::Type;
	using KBType          = KadanoffBaym<ComplexOrRealType>;
	using InputNgType     = PsimagLite::InputNg<CincuentaInputCheck>;
	using ModelParamsType = ModelParams<RealType>;

	virtual ~ImpuritySolverNeqBase() { }

	/*!
	 * \brief solve
	 * Run the impurity solver for the given bath parameters.
	 * For fixed-bath solvers (ExactDiag): diagonalizes H(U_i) and H(U_f) so that
	 * subsequent computeGimp(n) calls can fill time slices on demand.
	 *
	 * \param[in] bathParams Bath parameters passed to the solver.
	 */
	virtual void solve(const VectorRealType& bathParams) = 0;

	/*!
	 * \brief computeGimp
	 * Fill G_imp KB components for all (t_n, t_j) with j <= n, and G^{Left}(t_n, tau_j).
	 * Called once per time step inside the neq self-consistency loop.
	 *
	 * \param[in/out] gimp Two-time Green's function container to fill.
	 * \param[in] n Current time-step index.
	 */
	virtual void computeGimp(KBType& gimp, int n) const = 0;

	// Called after updateDelta fills row n of Delta, before the corrector computeGimp(n).
	// Default: no-op (ExactDiag has a fixed bath that doesn't update per step).
	// GBEK overrides this to update the Cholesky bath decomposition for step n.
	virtual void prepareTimeStep(int /*n*/, const KBType& /*delta*/) { }

	// Write -i*Delta^+_<(t_n,t_j) = sum_p V(n,p)*conj(V(j,p)) to file.
	// No-op for solvers that have no second bath.
	virtual void dumpPlusBath(const std::string& /*filename*/) const { }

	// Write the raw Cholesky factor V(n,p) itself (not the reconstructed
	// product) to file, for row-by-row comparison against an independent
	// offline trace. No-op for solvers that have no second bath.
	virtual void dumpV(const std::string& /*filename*/) const { }

	// Write one line per time step, "t docc Ekin Eint Etot", for double
	// occupation and energy-conservation observables. No-op for solvers
	// that don't compute these (currently only GBEK does).
	virtual void dumpDoccAndEnergy(const std::string& /*filename*/) const { }

	/*!
	 * \brief gimp
	 * Access the last-computed two-time impurity GF.
	 *
	 * \return Const reference to the KadanoffBaym Green's function.
	 */
	virtual const KBType& gimp() const = 0;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_BASE_H
