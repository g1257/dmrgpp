// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance18Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance18Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance18Type> >  >,Dmrg::Basis<SparseMatrixInstance18Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance18Type
 >
> MatrixVector18Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance18Type::RealType>,
	MatrixVector18Type, MatrixVector18Type::VectorType> LanczosSolver18Type;

template void mainLoop4<LanczosSolver18Type,Dmrg::VectorWithOffset<RealType> >
(LanczosSolver18Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance19Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance19Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance19Type> >  >,Dmrg::Basis<SparseMatrixInstance19Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance19Type
 >
> MatrixVector19Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance19Type::RealType>,
	MatrixVector19Type, MatrixVector19Type::VectorType> LanczosSolver19Type;

template void mainLoop4<LanczosSolver19Type,Dmrg::VectorWithOffset<RealType> >
(LanczosSolver19Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

