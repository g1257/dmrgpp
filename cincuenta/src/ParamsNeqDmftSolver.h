#ifndef PARAMS_NEQ_DMFT_SOLVER_H
#define PARAMS_NEQ_DMFT_SOLVER_H
#include "CincuentaInputCheck.h"
#include "ParamsDmftSolver.h"

namespace Dmft {

/*!
 * \brief ParamsNeqDmftSolver
 * Parameters for the non-equilibrium DMFT extension.
 * Contains all equilibrium parameters (via composition) plus the
 * interaction-quench and real-time grid parameters.
 */
template <typename ComplexOrRealType> struct ParamsNeqDmftSolver {

	using RealType     = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType  = PsimagLite::InputNg<CincuentaInputCheck>;
	using EqParamsType = ParamsDmftSolver<ComplexOrRealType>;

	ParamsNeqDmftSolver(typename InputNgType::Readable& io)
	    : eqParams(io)
	{
		io.readline(uInitial, "HubbardU=");
		io.readline(uFinal, "HubbardUFinal=");
		io.readline(tMax, "TmaxNeq=");
		io.readline(nT, "NtNeq=");
		dt = tMax / static_cast<RealType>(nT);

		try {
			io.readline(neqDmftIter, "NeqDmftIter=");
		} catch (std::exception&) {
			neqDmftIter = eqParams.dmftIter;
		}

		try {
			io.readline(neqDmftError, "NeqDmftTolerance=");
		} catch (std::exception&) {
			neqDmftError = eqParams.dmftError;
		}

		try {
			io.readline(neqBathRank, "NeqBathRank=");
		} catch (std::exception&) { }

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

		try {
			int tmp = 0;
			io.readline(tmp, "NeqAtomicLimit=");
			neqAtomicLimit = (tmp > 0);
		} catch (std::exception&) { }
	}

	EqParamsType eqParams; ///< Equilibrium DMFT parameters (beta, mu, nBath, etc.)

	RealType uInitial = 0; ///< Interaction quench: U_i -> U_f at t = 0
	RealType uFinal   = 0; ///<

	RealType tMax = 0; ///< Real-time grid
	SizeType nT   = 0; ///<
	RealType dt   = 0; ///<

	SizeType neqDmftIter  = 10; ///< Inner DMFT convergence at each time step
	RealType neqDmftError = 1e-4; ///<

	/// Rank L of the low-rank Cholesky second bath (GBEK scheme).
	/// L=0: first bath only (no GBEK second bath; equivalent to single-shot ExactDiag).
	/// L>0: second bath with 2L orbitals (L empty + L occupied at t=0).
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

	/// True atomic limit start (GBEK PRB 88, 235106 (2013), Sec. VI): no
	/// hopping at all at t=0, so the lattice hopping t*(t) ramp must start
	/// from exactly 0, not from whatever small-but-nonzero bandwidth the
	/// LatticeGf= equilibrium spec happens to use (that spec is only kept
	/// non-degenerate to keep the discarded equilibrium bath-fit minimizer
	/// well-posed; it should not leak into the real-time quench ramp's
	/// baseline -- see NeqLatticeGf's tStar_ construction).
	bool neqAtomicLimit = false;
};

} // namespace Dmft
#endif // PARAMS_NEQ_DMFT_SOLVER_H
