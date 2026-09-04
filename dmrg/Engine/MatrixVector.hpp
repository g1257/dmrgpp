#ifndef DMRG_MATRIX_VECTOR_HPP
#define DMRG_MATRIX_VECTOR_HPP

#include "MatrixVectorKron/MatrixVectorKron.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include <memory>

namespace Dmrg {

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

	MatrixVector(const ModelType&                 model,
	             const HamiltonianConnectionType& hc,
	             const AuxType&                   aux)
	    : ptr_(create(model, hc, aux))
	{ }

	SizeType rows() const { return ptr_->rows(); }

	SizeType cols() const { return ptr_->cols(); }

	void matrixVectorProduct(VectorType& x, const VectorType& y) const
	{
		ptr_->matrixVectorProduct(x, y);
	}

	void fullDiag(VectorRealType& eigs, FullMatrixType& fm) const { ptr_->fullDiag(eigs, fm); }

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
