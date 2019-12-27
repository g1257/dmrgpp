// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance4Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance4Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance4Type> >,Dmrg::Basis<SparseMatrixInstance4Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance4Type
 >
> MatrixVector4Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance4Type::RealType>,
	MatrixVector4Type, MatrixVector4Type::VectorType> LanczosSolver4Type;

template void mainLoop4<LanczosSolver4Type,Dmrg::VectorWithOffsets<RealType, Dmrg::Qn> >
(LanczosSolver4Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance5Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance5Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance5Type> >,Dmrg::Basis<SparseMatrixInstance5Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance5Type
 >
> MatrixVector5Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance5Type::RealType>,
	MatrixVector5Type, MatrixVector5Type::VectorType> LanczosSolver5Type;

template void mainLoop4<LanczosSolver5Type,Dmrg::VectorWithOffsets<RealType, Dmrg::Qn> >
(LanczosSolver5Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

