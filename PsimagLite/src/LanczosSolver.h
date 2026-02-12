#ifndef PSI_LANCZOS_SOLVER_H
#define PSI_LANCZOS_SOLVER_H
#include "LanczosCore.h"
#include "MatrixSolverBase.hh"
#include "Profiling.h"
#include "Vector.h"

namespace PsimagLite {

template <typename MatrixType_> class LanczosSolver : public MatrixSolverBase<MatrixType_> {

public:

	using MatrixType            = MatrixType_;
	using BaseType              = MatrixSolverBase<MatrixType_>;
	using ComplexOrRealType     = typename BaseType::ComplexOrRealType;
	using VectorType            = typename BaseType::VectorType;
	using RealType              = typename BaseType::RealType;
	using VectorRealType        = typename BaseType::VectorRealType;
	using VectorVectorType      = typename BaseType::VectorVectorType;
	using ParametersSolverType  = ParametersForSolver<RealType>;
	using LanczosCoreType       = LanczosCore<ParametersSolverType, MatrixType_, VectorType>;
	using TridiagonalMatrixType = typename LanczosCoreType::TridiagonalMatrixType;
	using PostProcType          = ContinuedFraction<TridiagonalMatrixType>;

	LanczosSolver(const MatrixType& mat, const ParametersSolverType& params)
	    : ls_(mat, params, BaseType::isReorthoEnabled(params.options, params.lotaMemory))
	{ }

	void computeOneState(RealType&         energy,
	                     VectorType&       z,
	                     const VectorType& initialVector,
	                     SizeType          excited)
	{
		Profiling profiling("LanczosSolver", std::cout);

		TridiagonalMatrixType ab;
		ls_.decomposition(initialVector, ab, excited);

		VectorRealType                            eigs(ab.size());
		typename LanczosCoreType::DenseMatrixType ritz;
		ab.buildDenseMatrix(ritz);
		diag(ritz, eigs, 'V');

		energy = eigs[excited];
		ls_.excitedVector(z, ritz, excited);

		String str = "LanczosSolver: computeOneState: ";
		if (norm(z) < 1e-6)
			throw RuntimeError(str + " norm is zero\n");

		const RealType norma = norm(initialVector);
		const SizeType iter  = ls_.steps();

		if (norma < 1e-5 || norma > 100)
			std::cerr << "norma=" << norma << "\n";

		OstringStream msg(std::cout.precision());
		String        what = "lowest";
		if (excited > 0)
			what = ttos(excited) + " excited";
		msg() << "Found " << what << " eigenvalue= " << energy << " after " << iter;
		msg() << " iterations, "
		      << " orig. norm=" << norma << " excited=" << excited;
		profiling.end(msg().str());
	}

	void computeAllStatesBelow(VectorRealType&   eigs,
	                           VectorVectorType& z,
	                           const VectorType& initialVector,
	                           SizeType          excited)
	{
		TridiagonalMatrixType ab;
		ls_.decomposition(initialVector, ab, excited);

		if (excited > ab.size())
			throw RuntimeError("Excited too big\n");

		typename LanczosCoreType::DenseMatrixType ritz;
		ab.buildDenseMatrix(ritz);
		diag(ritz, eigs, 'V');

		SizeType n = z.size();
		if (n > excited + 1)
			n = excited + 1;
		for (SizeType i = 0; i < n; ++i)
			ls_.excitedVector(z[i], ritz, i);
	}

	void decomposition(const VectorType& initVector, TridiagonalMatrixType& ab)
	{
		return ls_.decomposition(initVector, ab, ls_.params().eigsForStop);
	}

	void lanczosVectorsSwap(typename LanczosCoreType::DenseMatrixType& V)
	{
		ls_.lanczosVectorsSwap(V);
	}

	SizeType steps() const { return ls_.steps(); }

private:

	LanczosCoreType ls_;
};
} // namespace PsimagLite
#endif // PSI_LANCZOS_SOLVER_H
