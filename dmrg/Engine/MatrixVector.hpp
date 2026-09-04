//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   dmrg/Engine/MatrixVector.hpp
 * \brief  Runtime-selectable matrix-vector facade
 */
//---------------------------------------------------------------------------//

#ifndef DMRG_MATRIX_VECTOR_HPP
#define DMRG_MATRIX_VECTOR_HPP

#include "MatrixVectorKron/MatrixVectorKron.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include <memory>

namespace Dmrg {

//===========================================================================//
/*!
 * \class MatrixVector
 *
 * \brief Own and expose the configured matrix-vector implementation
 *
 * This facade preserves the matrix interface required by the eigensolvers
 * while selecting MatrixVectorStored, MatrixVectorOnTheFly, or
 * MatrixVectorKron at runtime. The selected implementation is owned by the
 * facade and is valid for the facade's lifetime.
 *
 * \tparam ComplexOrRealType_ Scalar type used by vectors and matrices
 */
//===========================================================================//
template <typename ComplexOrRealType_> class MatrixVector {

public:

	using BaseType                  = MatrixVectorBase<ComplexOrRealType_>;
	using TypesType                 = typename BaseType::TypesType;
	using ModelType                 = typename BaseType::ModelType;
	using ParametersType            = typename ModelType::ParametersType;
	using ModelHelperType           = typename BaseType::ModelHelperType;
	using RealType                  = typename BaseType::RealType;
	using SparseMatrixType          = typename BaseType::SparseMatrixType;
	using ComplexOrRealType         = typename BaseType::ComplexOrRealType;
	using value_type                = typename BaseType::value_type;
	using VectorRealType            = typename BaseType::VectorRealType;
	using VectorType                = typename BaseType::VectorType;
	using FullMatrixType            = typename BaseType::FullMatrixType;
	using HamiltonianConnectionType = typename ModelType::HamiltonianConnectionType;
	using AuxType                   = typename ModelHelperType::Aux;

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Construct the configured matrix-vector implementation
	 *
	 * MatrixVectorStored takes precedence over MatrixVectorOnTheFly; when
	 * neither option is set, MatrixVectorKron is selected.
	 *
	 * \param[in] model Model providing parameters and matrix-vector operations
	 * \param[in] hc    Hamiltonian connection for the current sector
	 * \param[in] aux   Auxiliary data identifying the current sector
	 *
	 * 	hrows PsimagLite::RuntimeError if ArnoldiSaI is requested without
	 *         MatrixVectorStored
	 */
	MatrixVector(const ModelType&                 model,
	             const HamiltonianConnectionType& hc,
	             const AuxType&                   aux)
	    : ptr_(create(model, hc, aux))
	{ }

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Return the number of matrix rows
	 *
	 * \returns The dimension of the current sector
	 */
	SizeType rows() const { return ptr_->rows(); }

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Return the number of matrix columns
	 *
	 * \returns The dimension of the current sector
	 */
	SizeType cols() const { return ptr_->cols(); }

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Apply the represented matrix to a vector
	 *
	 * \param[out] x Result of the matrix-vector product
	 * \param[in]  y Input vector
	 */
	void matrixVectorProduct(VectorType& x, const VectorType& y) const
	{
		ptr_->matrixVectorProduct(x, y);
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Diagonalize the represented matrix as a dense matrix
	 *
	 * \param[out] eigs Eigenvalues in ascending order
	 * \param[out] fm   Eigenvectors stored as columns
	 */
	void fullDiag(VectorRealType& eigs, FullMatrixType& fm) const { ptr_->fullDiag(eigs, fm); }

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Return the explicitly stored sparse matrix
	 *
	 * \returns A reference owned by the selected implementation
	 *
	 * 	hrows PsimagLite::RuntimeError if the selected implementation does not
	 *         provide an explicit sparse matrix
	 */
	const SparseMatrixType& toCRS() const { return ptr_->toCRS(); }

private:

	static std::unique_ptr<BaseType>
	create(const ModelType& model, const HamiltonianConnectionType& hc, const AuxType& aux)
	{
		if (model.params().options.isSet("MatrixVectorStored"))
			return std::make_unique<MatrixVectorStored<ComplexOrRealType>>(
			    model, hc, aux);

		if (model.params().matrix_solver_enum
		    == ParametersType::MatrixSolverEnum::ARNOLDISAI) {
			throw PsimagLite::RuntimeError(
			    "MatrixSolver=ArnoldiSaI requires MatrixVectorStored\n");
		}

		if (model.params().options.isSet("MatrixVectorOnTheFly"))
			return std::make_unique<MatrixVectorOnTheFly<ComplexOrRealType>>(
			    model, hc, aux);

		return std::make_unique<MatrixVectorKron<ComplexOrRealType>>(model, hc, aux);
	}

	std::unique_ptr<BaseType> ptr_;
};

} // namespace Dmrg

#endif // DMRG_MATRIX_VECTOR_HPP
