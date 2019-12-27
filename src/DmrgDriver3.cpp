// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance6Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance6Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance6Type> >,Dmrg::Basis<SparseMatrixInstance6Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance6Type
 >
> MatrixVector6Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance6Type::RealType>,
	MatrixVector6Type, MatrixVector6Type::VectorType> LanczosSolver6Type;

template void mainLoop4<LanczosSolver6Type,Dmrg::VectorWithOffset<RealType, Dmrg::Qn> >
(LanczosSolver6Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance7Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance7Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance7Type> >,Dmrg::Basis<SparseMatrixInstance7Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance7Type
 >
> MatrixVector7Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance7Type::RealType>,
	MatrixVector7Type, MatrixVector7Type::VectorType> LanczosSolver7Type;

template void mainLoop4<LanczosSolver7Type,Dmrg::VectorWithOffset<RealType, Dmrg::Qn> >
(LanczosSolver7Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

