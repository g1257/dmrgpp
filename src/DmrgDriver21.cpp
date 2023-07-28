// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"

typedef PsimagLite::CrsMatrix<std::complex<RealType>> SparseMatrixInstance42Type;
typedef Dmrg::SuperGeometry<std::complex<RealType>, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance42Type;

typedef Dmrg::MatrixVectorOnTheFly<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance42Type>>, Dmrg::Basis<SparseMatrixInstance42Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance42Type>>
    MatrixVector42Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance42Type::RealType>,
    MatrixVector42Type,
    MatrixVector42Type::VectorType>
    LanczosSolver42Type;

template void mainLoop4<LanczosSolver42Type, Dmrg::VectorWithOffset<std::complex<RealType>, Dmrg::Qn>>(LanczosSolver42Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);

typedef PsimagLite::CrsMatrix<std::complex<RealType>> SparseMatrixInstance43Type;
typedef Dmrg::SuperGeometry<std::complex<RealType>, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance43Type;

typedef Dmrg::MatrixVectorStored<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance43Type>>, Dmrg::Basis<SparseMatrixInstance43Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance43Type>>
    MatrixVector43Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance43Type::RealType>,
    MatrixVector43Type,
    MatrixVector43Type::VectorType>
    LanczosSolver43Type;

template void mainLoop4<LanczosSolver43Type, Dmrg::VectorWithOffset<std::complex<RealType>, Dmrg::Qn>>(LanczosSolver43Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);
