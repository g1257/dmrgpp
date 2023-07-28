// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"

typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance26Type;
typedef Dmrg::SuperGeometry<RealType, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance26Type;

typedef Dmrg::MatrixVectorKron<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance26Type>>, Dmrg::Basis<SparseMatrixInstance26Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance26Type>>
    MatrixVector26Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance26Type::RealType>,
    MatrixVector26Type,
    MatrixVector26Type::VectorType>
    LanczosSolver26Type;

template void mainLoop4<LanczosSolver26Type, Dmrg::VectorWithOffset<std::complex<RealType>, Dmrg::Qn>>(LanczosSolver26Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);

typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance27Type;
typedef Dmrg::SuperGeometry<RealType, PsimagLite::InputNg<Dmrg::InputCheck>::Readable, Dmrg::ProgramGlobals> GeometryInstance27Type;

typedef Dmrg::MatrixVectorOnTheFly<
    Dmrg::ModelBase<
	Dmrg::ModelHelperLocal<
	    Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance27Type>>, Dmrg::Basis<SparseMatrixInstance27Type>>>,
	ParametersDmrgSolverType,
	InputNgType::Readable,
	GeometryInstance27Type>>
    MatrixVector27Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance27Type::RealType>,
    MatrixVector27Type,
    MatrixVector27Type::VectorType>
    LanczosSolver27Type;

template void mainLoop4<LanczosSolver27Type, Dmrg::VectorWithOffsets<std::complex<RealType>, Dmrg::Qn>>(LanczosSolver27Type::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);
