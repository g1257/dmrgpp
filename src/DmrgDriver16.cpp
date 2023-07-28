// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"

typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance32Type;
typedef Dmrg::SuperGeometry<RealType, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance32Type;

typedef Dmrg::MatrixVectorKron<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance32Type>>, Dmrg::Basis<SparseMatrixInstance32Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance32Type>>
    MatrixVector32Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance32Type::RealType>,
    MatrixVector32Type,
    MatrixVector32Type::VectorType>
    LanczosSolver32Type;

template void mainLoop4<LanczosSolver32Type, Dmrg::VectorWithOffset<std::complex<RealType>, Dmrg::Qn>>(LanczosSolver32Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);

typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance33Type;
typedef Dmrg::SuperGeometry<RealType, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance33Type;

typedef Dmrg::MatrixVectorOnTheFly<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance33Type>>, Dmrg::Basis<SparseMatrixInstance33Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance33Type>>
    MatrixVector33Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance33Type::RealType>,
    MatrixVector33Type,
    MatrixVector33Type::VectorType>
    LanczosSolver33Type;

template void mainLoop4<LanczosSolver33Type, Dmrg::VectorWithOffsets<std::complex<RealType>, Dmrg::Qn>>(LanczosSolver33Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);
