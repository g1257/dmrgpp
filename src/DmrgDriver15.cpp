// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance30Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance30Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance30Type> >,Dmrg::Basis<SparseMatrixInstance30Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance30Type
 >
> MatrixVector30Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance30Type::RealType>,
	MatrixVector30Type, MatrixVector30Type::VectorType> LanczosSolver30Type;

template void mainLoop4<LanczosSolver30Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver30Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance31Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance31Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance31Type> >,Dmrg::Basis<SparseMatrixInstance31Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance31Type
 >
> MatrixVector31Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance31Type::RealType>,
	MatrixVector31Type, MatrixVector31Type::VectorType> LanczosSolver31Type;

template void mainLoop4<LanczosSolver31Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver31Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

