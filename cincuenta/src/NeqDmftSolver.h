#ifndef NEQ_DMFT_SOLVER_H
#define NEQ_DMFT_SOLVER_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqGBEK.h"
#include "NeqHybridizationTarget.h"
#include "NeqRealTimeGf.h"
#include "ParamsNeqDmftSolver.h"
#include <iostream>

namespace Dmft {

/*!
 * \brief NeqDmftSolver
 * Non-equilibrium DMFT self-consistency loop for an interaction quench U_i -> U_f.
 *
 * At each step the GBEK progressive scheme predicts G_imp, makes the Bethe
 * hybridization target Λ(n,j) = v(n)v(j)G_imp(n,j), updates the bath
 * decomposition, then corrects G_imp with that bath. There is no propagated
 * hybridization-target Green's function in this path.
 *
 * ImpSolverTemplate selects the impurity solver.  Production uses
 * ImpuritySolverNeqGBEK, the two-bath Cholesky scheme (GBEK PRB 88, 235106).
 */
template <typename ComplexOrRealType,
          template <typename> class ImpSolverTemplate = ImpuritySolverNeqGBEK>
class NeqDmftSolver {

public:

	using RealType                = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType             = std::complex<RealType>;
	using VectorRealType          = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType       = typename PsimagLite::Vector<ComplexType>::Type;
	using RealTimeGfType          = NeqRealTimeGf<ComplexOrRealType>;
	using ParamsNeqType           = ParamsNeqDmftSolver<ComplexOrRealType>;
	using InputNgType             = PsimagLite::InputNg<CincuentaInputCheck>;
	using ImpSolverType           = ImpSolverTemplate<ComplexOrRealType>;
	using HybridizationTargetType = NeqHybridizationTarget<ComplexOrRealType>;

	NeqDmftSolver(const ParamsNeqType& params, typename InputNgType::Readable& io)
	    : params_(params)
	    , impSolver_(params, io)
	    , hybridizationTarget_(params)
	    , gimp_(params.nT, params.dt)
	{ }

	/*!
	 * \brief Run the full neq-DMFT calculation.
	 *
	 * \param[in] bathParams Equilibrium bath parameters {V_0..V_{nBath-1}, ε_0..ε_{nBath-1}}
	 *                       obtained from the preceding equilibrium DMFT run.
	 */
	void solve(const VectorRealType& bathParams)
	{
		std::cout << "NeqDmftSolver: running impurity solver setup\n";
		impSolver_.solve(bathParams);

		// Populate t=0 and seed the Cholesky decomposition.  The target is
		// established before prepareTimeStep(0); it is never propagated as a
		// Green's function.
		impSolver_.computeGimp(gimp_, 0);
		hybridizationTarget_.updateTimeRow(0, gimp_);
		impSolver_.prepareTimeStep(0, hybridizationTarget_.lambda());

		std::cout << "NeqDmftSolver: starting time propagation to t_max=" << params_.tMax
		          << " with nT=" << params_.nT << " steps\n";

		for (int n = 1; n <= static_cast<int>(params_.nT); ++n) {
			timeStep(n);
			if (n % 10 == 0 || n == static_cast<int>(params_.nT))
				std::cout << "  step " << n << " / " << params_.nT << "\n";
		}

		std::cout << "NeqDmftSolver: done\n";
	}

	/*!
	 * \brief Access the impurity GF (populated after solve()).
	 */
	const RealTimeGfType& gimp() const { return gimp_; }

	/*!
	 * \brief Write retained real-time Green's functions to files.
	 *
	 * When params_.neqOutputPrefix is set, filenames are "{prefix}-green-retarded"
	 * etc.; otherwise "green-retarded" etc.
	 */
	void dumpGreenFunctions() const
	{
		const std::string& p = params_.neqOutputPrefix;
		gimp_.dump(p.empty() ? "green" : p + "-green");
		hybridizationTarget_.lambda().dump(p.empty() ? "lambda" : p + "-lambda");
		impSolver_.dumpPlusBath(p.empty() ? "plus-bath-lesser" : p + "-plus-bath-lesser");
		impSolver_.dumpV(p.empty() ? "cholesky-V" : p + "-cholesky-V");
		impSolver_.dumpDoccAndEnergy(p.empty() ? "docc-energy" : p + "-docc-energy");
	}

private:

	// Advance all real-time components by one time step n.
	//
	// Self-consistency following GBEK Fig. 2(b) progressive scheme (PRB 88, 235106):
	//   Predictor: V[n] is a zero-row cold start until prepareTimeStep updates it.
	//   Corrector iterations: updateTimeRow fills Λ(n,j), prepareTimeStep updates
	//   the Cholesky bath V[n] from that complete target row, then computeGimp
	//   re-evaluates G_imp. No post-corrector target propagation is performed.
	void timeStep(int n)
	{
		// Predictor: G_imp(n,j) with V[n]'s zero-row cold start.
		impSolver_.computeGimp(gimp_, n);

		for (SizeType iter = 0; iter < params_.neqDmftIter; ++iter) {
			// Λ(n,j) = v(n)v(j)G_imp(n,j): fill the complete target row.
			hybridizationTarget_.updateTimeRow(n, gimp_);

			// Update bath from that complete target row.
			impSolver_.prepareTimeStep(n, hybridizationTarget_.lambda());

			// Corrector: re-evaluate G_imp with the updated bath.
			// This final evaluation makes G^< consistent with V[n], rather than
			// with the penultimate bath update.
			impSolver_.computeGimp(gimp_, n);
		}
	}

	const ParamsNeqType&    params_;
	ImpSolverType           impSolver_;
	HybridizationTargetType hybridizationTarget_;
	RealTimeGfType          gimp_; ///< local copy filled step by step
};

} // namespace Dmft
#endif // NEQ_DMFT_SOLVER_H
