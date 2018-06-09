// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance24Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance24Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance24Type> >  >,Dmrg::Basis<SparseMatrixInstance24Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance24Type
 >
> MatrixVector24Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance24Type::RealType>,
	MatrixVector24Type, MatrixVector24Type::VectorType> LanczosSolver24Type;

template void mainLoop4<LanczosSolver24Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::EffectiveQuantumNumber<RealType> > >
(LanczosSolver24Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance25Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance25Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance25Type> >  >,Dmrg::Basis<SparseMatrixInstance25Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance25Type
 >
> MatrixVector25Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance25Type::RealType>,
	MatrixVector25Type, MatrixVector25Type::VectorType> LanczosSolver25Type;

template void mainLoop4<LanczosSolver25Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::EffectiveQuantumNumber<RealType> > >
(LanczosSolver25Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

