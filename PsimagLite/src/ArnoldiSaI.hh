#ifndef ARNOLDISAI_HH
#define ARNOLDISAI_HH
#include "ArnoldiIteration.hh"
#include "CrsMatrix.h"
#include "MatrixSolverBase.hh"
#include "ParametersForSolver.h"
#include <vector>

namespace PsimagLite {

/* The template parameter MatrixType_ allows the user to pass an object
 * so that the matrix to be decomposed doesn't have to be in memory,
 * only the A*v operation is needed.
 * But in practice, this is limited by the need for an inversion
 */
template <typename MatrixType_> class ArnoldiSaI : public MatrixSolverBase<MatrixType_> {

public:

	using MatrixType           = MatrixType_;
	using BaseType             = MatrixSolverBase<MatrixType_>;
	using ComplexOrRealType    = typename BaseType::ComplexOrRealType;
	using VectorType           = typename BaseType::VectorType;
	using RealType             = typename BaseType::RealType;
	using VectorRealType       = typename BaseType::VectorRealType;
	using VectorVectorType     = typename BaseType::VectorVectorType;
	using ArnoldiIterationType = ArnoldiIteration<CrsMatrix<ComplexOrRealType>>;
	using ParametersSolverType = ParametersForSolver<RealType>;

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Constructor
	 *
	 * \param[in]  a       The matrix for the Arnolid shift-and-invert method
	 * \param[in]  params  The parameters object; see ParametersForSolver.h
	 * \param[in]  sigma   The value for the small shift, must be positive
	 */
	ArnoldiSaI(const MatrixType& a, const ParametersSolverType& params, const RealType& sigma)
	    : arnoldi_iteration_(params)
	    , a_(a)
	    , sigma_(sigma)
	{
		if (a_.rows() != a_.cols()) {
			throw RuntimeError("computeOneState::ctor(): Matrix not square\n");
		}

		if (sigma_ <= 0.) {
			throw RuntimeError("computeOneState::ctor(): sigma must be positive\n");
		}
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Compute the lowest real positive eigenvector of A
	 *
	 * Note that the matrix A is passed in the ctor.
	 * This uses the shift-and-invert method followed by Arnoldi
	 *
	 * \param[out]  eigenvalue  The eigenvalue (because it's not zero actually)
	 * \param[out]  eigenvector The eigenvector
	 * \param[in]   init_vector The initial vector for Arnoldi
	 * \param[in]   index       Must be zero
	 */
	void computeOneState(RealType&         eigenvalue,
	                     VectorType&       eigenvector,
	                     const VectorType& init_vector,
	                     SizeType) final
	{
		std::vector<VectorType>   q;
		Matrix<ComplexOrRealType> h;

		// Here lies the problem:
		// I only know how to invert dense matrices
		// Needs more study: FIXME TODO

		const CrsMatrix<ComplexOrRealType>& crs    = a_.toCRS();
		Matrix<ComplexOrRealType>           a_copy = crs.toDense();

		// shift first
		SizeType nrows = a_.rows();

		// we check this at runtime in the ctor:
		assert(nrows == a_.cols());

		for (SizeType i = 0; i < nrows; ++i) {
			a_copy(i, i) -= sigma_;
		}

		// invert now
		inverse(a_copy);

		CrsMatrix<ComplexOrRealType> a_sai(a_copy);
		arnoldi_iteration_.decompose(q, h, a_sai, init_vector);

		// eigs_of_h are always complex even if h is real
		std::vector<std::complex<RealType>> eigs_of_h;

		Matrix<ComplexOrRealType> h_right_eigenvectors;
		arnoldi_iteration_.diagHessenberg(eigs_of_h, h_right_eigenvectors, h);

		SizeType index_wanted = 0; // we want the first
		if (eigs_of_h.size() < index_wanted) {
			throw RuntimeError("Wanted eigenvalue index " + ttos(index_wanted)
			                   + " but there are only " + ttos(eigs_of_h.size())
			                   + " values.\n");
		}

		// Eigenvalues is always complex even if h is real.
		std::complex<RealType> eigs_of_sai = eigs_of_h[index_wanted];
		if (eigs_of_sai == 0.) {
			// We shifted by sigma so that A - sigmaI be invertible
			// If this isn't the case, we need to throw.
			eigenvalue            = 0.;
			std::string warn_zero = "Eigenvalue of the shifted-and-inverted is 0\n";
			std::cerr << warn_zero;
			std::cout << warn_zero;
		} else {
			// Recover the original eigenvalue
			eigenvalue = std::real(sigma_ + 1. / eigs_of_sai);
		}

		// The eigenvector is the same
		eigenvector
		    = arnoldi_iteration_.getEigenvector(h_right_eigenvectors, q, index_wanted);
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Compute lowest eigenvectors and eigenvalues below nexcited
	 *
	 * Defers to computeOneState if nexcited = 1, else throws.
	 *
	 * Note that the matrix A is passed in the ctor.
	 *
	 * \param[out]  eigs      The eigenvalues
	 * \param[out]  zs        The eigenvectors
	 * \param[in]   init_v    The initial vector for Arnoldi
	 * \param[in]   nexcited  Must be one or throws
	 */
	void computeAllStatesBelow(VectorRealType&   eigs,
	                           VectorVectorType& zs,
	                           const VectorType& init_v,
	                           SizeType          nexcited) final
	{
		if (nexcited != 1) {
			throw RuntimeError("ArnoldiSaI::computeAllStatesBelow() only for g.s.\n");
		} else {
			RealType   eigenvalue = 0;
			VectorType z;
			computeOneState(eigenvalue, z, init_v, nexcited);
			eigs.clear();
			zs.clear();
			eigs.push_back(eigenvalue);
			zs.push_back(z);
		}
	}

private:

	ArnoldiIterationType arnoldi_iteration_;
	const MatrixType&    a_;
	RealType             sigma_;
};

}
#endif // ARNOLDISAI_HH
