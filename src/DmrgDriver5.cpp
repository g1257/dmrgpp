// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"

typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance10Type;
typedef Dmrg::SuperGeometry<RealType, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance10Type;

typedef Dmrg::MatrixVectorStored<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance10Type>>, Dmrg::Basis<SparseMatrixInstance10Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance10Type>>
    MatrixVector10Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance10Type::RealType>,
    MatrixVector10Type,
    MatrixVector10Type::VectorType>
    LanczosSolver10Type;

template void mainLoop4<LanczosSolver10Type, Dmrg::VectorWithOffsets<RealType, Dmrg::Qn>>(LanczosSolver10Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);

typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance11Type;
typedef Dmrg::SuperGeometry<RealType, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance11Type;

typedef Dmrg::MatrixVectorKron<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance11Type>>, Dmrg::Basis<SparseMatrixInstance11Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance11Type>>
    MatrixVector11Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance11Type::RealType>,
    MatrixVector11Type,
    MatrixVector11Type::VectorType>
    LanczosSolver11Type;

template void mainLoop4<LanczosSolver11Type, Dmrg::VectorWithOffsets<RealType, Dmrg::Qn>>(LanczosSolver11Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);
