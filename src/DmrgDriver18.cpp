// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance36Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance36Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance36Type> >  >,Dmrg::Basis<SparseMatrixInstance36Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance36Type
 >
> MatrixVector36Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance36Type::RealType>,
	MatrixVector36Type, MatrixVector36Type::VectorType> LanczosSolver36Type;

template void mainLoop4<LanczosSolver36Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver36Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance37Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance37Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance37Type> >  >,Dmrg::Basis<SparseMatrixInstance37Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance37Type
 >
> MatrixVector37Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance37Type::RealType>,
	MatrixVector37Type, MatrixVector37Type::VectorType> LanczosSolver37Type;

template void mainLoop4<LanczosSolver37Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver37Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

