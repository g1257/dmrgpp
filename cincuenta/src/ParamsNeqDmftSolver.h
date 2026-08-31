#ifndef PARAMS_NEQ_DMFT_SOLVER_H
#define PARAMS_NEQ_DMFT_SOLVER_H
#include "CincuentaInputCheck.h"
#include "ParamsMatsubaraGrid.h"

namespace Dmft {

/*!
 * \brief ParamsNeqDmftSolver
 * Parameters for the non-equilibrium DMFT extension: the interaction-quench
 * and real-time grid parameters, plus (via composition) the Matsubara-grid
 * parameters shared with the equilibrium stage.
 *
 * Deliberately does NOT compose the full equilibrium ParamsDmftSolver: under
 * NeqAtomicLimit=1 there is no equilibrium DMFT stage at all (see
 * cincuenta.cpp), so equilibrium-bath-fit-only parameters (bath size,
 * impurity solver choice, fit options, minimizer settings) have no meaning
 * here and are never parsed by this class.
 */
template <typename ComplexOrRealType> struct ParamsNeqDmftSolver {

	using RealType       = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType    = PsimagLite::InputNg<CincuentaInputCheck>;
	using GridParamsType = ParamsMatsubaraGrid<ComplexOrRealType>;

	ParamsNeqDmftSolver(typename InputNgType::Readable& io)
	    : neqAtomicLimit(peekNeqAtomicLimit(io))
	    , grid(io, neqAtomicLimit)
	{
		io.readline(uInitial, "HubbardU=");
		io.readline(uFinal, "HubbardUFinal=");
		io.readline(tMax, "TmaxNeq=");
		io.readline(nT, "NtNeq=");
		if (nT == 0)
			throw PsimagLite::RuntimeError("Non-equilibrium GBEK requires NtNeq>0; "
			                               "at least one real-time step is required\n");
		dt = tMax / static_cast<RealType>(nT);

		try {
			io.readline(neqDmftIter, "NeqDmftIter=");
		} catch (std::exception&) { }

		try {
			io.readline(neqBathRank, "NeqBathRank=");
		} catch (std::exception&) { }

		if (neqDmftIter == 0)
			throw PsimagLite::RuntimeError(
			    "Non-equilibrium GBEK requires NeqDmftIter>0; "
			    "at least one corrector iteration is required at each time step\n");
		if (neqBathRank == 0)
			throw PsimagLite::RuntimeError(
			    "Non-equilibrium GBEK requires NeqBathRank>0; "
			    "the frozen rank-zero bath path is no longer supported\n");
		if (neqBathRank > nT)
			throw PsimagLite::RuntimeError(
			    "Non-equilibrium GBEK requires NeqBathRank<=NtNeq; "
			    "every bath column must have a real-time seed\n");

		try {
			io.readline(bandwidthFinal, "BandwidthFinal=");
		} catch (std::exception&) { }

		try {
			io.readline(quenchShape, "QuenchShape=");
		} catch (std::exception&) { }

		try {
			io.readline(quenchDuration, "QuenchDuration=");
		} catch (std::exception&) { }

		try {
			io.readline(neqOutputPrefix, "NeqOutputPrefix=");
		} catch (std::exception&) { }
	}

	/// True atomic limit start (GBEK PRB 88, 235106 (2013), Sec. VI): no
	/// equilibrium bath at all, so the lattice hopping ramp starts from exactly
	/// 0 and the impurity solver starts from an empty bath (see cincuenta.cpp).
	/// Read before `grid` (declaration order) so `grid` knows whether
	/// LatticeGf= is even meaningful.
	const bool neqAtomicLimit;

	GridParamsType grid; ///< Matsubara-grid parameters shared with the equilibrium stage.

	RealType uInitial = 0; ///< Interaction quench: U_i -> U_f at t = 0
	RealType uFinal   = 0; ///<

	RealType tMax = 0; ///< Real-time grid
	SizeType nT   = 0; ///<
	RealType dt   = 0; ///<

	SizeType neqDmftIter = 10; ///< Fixed corrector count at each time step

	/// Rank L of the low-rank Cholesky second bath (GBEK scheme).
	/// Production neq-DMFT requires L>0; zero remains the omitted-input
	/// sentinel and is rejected by the production dispatcher.
	/// L>0 creates a second bath with 2L orbitals (L empty + L occupied at t=0).
	SizeType neqBathRank = 0;

	/// Hopping quench: Bethe lattice bandwidth for t > 0.
	/// 0 (default) means no quench — use the equilibrium bandwidth from LatticeGf.
	RealType bandwidthFinal = 0;

	/// Shape of the bandwidth ramp: "step" (default), "cosine", "tanh".
	std::string quenchShape = "step";

	/// Duration of the ramp in real time. 0 means instantaneous step at t=0.
	RealType quenchDuration = 0;

	/// Optional prefix for output Green's function files.
	/// Empty (default) → "green-retarded" etc.  Non-empty → "{prefix}-green-retarded" etc.
	std::string neqOutputPrefix = "";

private:

	static bool peekNeqAtomicLimit(typename InputNgType::Readable& io)
	{
		int tmp = 0;
		try {
			io.readline(tmp, "NeqAtomicLimit=");
		} catch (std::exception&) {
			return false;
		}
		return tmp > 0;
	}
};

} // namespace Dmft
#endif // PARAMS_NEQ_DMFT_SOLVER_H
