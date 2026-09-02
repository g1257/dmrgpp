#ifndef MATRIX_VECTOR_TYPES_HPP
#define MATRIX_VECTOR_TYPES_HPP

#include "Basis.h"
#include "BasisWithOperators.h"
#include "InputCheck.h"
#include "LeftRightSuper.h"
#include "ModelBase.h"
#include "ModelHelperLocal.h"
#include "ParametersDmrgSolver.h"
#include "Qn.h"
#include "SuperGeometry.h"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/InputNg.h>

namespace Dmrg {

template <typename ComplexOrRealType_> struct MatrixVectorTypes {

	using ComplexOrRealType         = ComplexOrRealType_;
	using InputNgType               = PsimagLite::InputNg<InputCheck>;
	using InputType                 = typename InputNgType::Readable;
	using SparseMatrixType          = PsimagLite::CrsMatrix<ComplexOrRealType>;
	using BasisType                 = Basis<SparseMatrixType>;
	using RealType                  = typename BasisType::RealType;
	using BasisWithOperatorsType    = BasisWithOperators<BasisType>;
	using LeftRightSuperType        = LeftRightSuper<BasisWithOperatorsType, BasisType>;
	using ModelHelperType           = ModelHelperLocal<LeftRightSuperType>;
	using ParametersType            = ParametersDmrgSolver<RealType, InputType, Qn>;
	using SuperGeometryType         = SuperGeometry<ComplexOrRealType, InputType, ProgramGlobals>;
	using ModelType = ModelBase<ModelHelperType,
	                            ParametersType,
	                            InputType,
	                            SuperGeometryType>;
	using HamiltonianConnectionType = typename ModelType::HamiltonianConnectionType;
};

} // namespace Dmrg

#endif // MATRIX_VECTOR_TYPES_HPP
