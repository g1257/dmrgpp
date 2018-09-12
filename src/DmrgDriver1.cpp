// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance2Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance2Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance2Type> >  >,Dmrg::Basis<SparseMatrixInstance2Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance2Type
 >
> MatrixVector2Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance2Type::RealType>,
	MatrixVector2Type, MatrixVector2Type::VectorType> LanczosSolver2Type;

template void mainLoop4<LanczosSolver2Type,Dmrg::VectorWithOffset<RealType, Dmrg::Qn> >
(LanczosSolver2Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance3Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance3Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance3Type> >  >,Dmrg::Basis<SparseMatrixInstance3Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance3Type
 >
> MatrixVector3Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance3Type::RealType>,
	MatrixVector3Type, MatrixVector3Type::VectorType> LanczosSolver3Type;

template void mainLoop4<LanczosSolver3Type,Dmrg::VectorWithOffsets<RealType, Dmrg::Qn> >
(LanczosSolver3Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

